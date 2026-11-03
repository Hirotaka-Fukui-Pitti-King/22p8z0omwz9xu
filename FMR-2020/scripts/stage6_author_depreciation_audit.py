#!/usr/bin/env python3
"""Independent first-period economic audit of the ORIGINAL saved HANKSOME MATLAB path.

Reads the authors' .mat workspace directly; does not call fsolve, MATLAB, or
a Julia general-equilibrium solver. The paper/source model is the reference.
Do not reinterpret a conditional scalar root as a distinct full equilibrium.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.io import loadmat
from scipy.optimize import brentq


def vector(value):
    return np.asarray(value, dtype=float).reshape(-1)


def parameter(p, name):
    if not hasattr(p, name):
        raise ValueError("Original MATLAB params missing " + name)
    return getattr(p, name)


def scalar(p, name):
    a = vector(parameter(p, name))
    if a.size != 1 or not np.isfinite(a[0]):
        raise ValueError("Invalid MATLAB scalar params." + name)
    return float(a[0])


def read_original(path):
    names = (
        "params", "kl_exp", "resid", "tot_exp", "G0_shock",
        "In_dis_exp", "c_tran_fine_exp", "dec_tran_fine_exp",
        "r_t_exp", "w_nr_guess_exp", "Profits_exp", "q_agg_exp",
        "y_agg_exp", "c_t_exp", "c_t_agg_exp", "net_c_exp",
        "new_shock_b_agg", "z_t",
    )
    data = loadmat(str(path), variable_names=names,
                   squeeze_me=True, struct_as_record=False)
    absent = [n for n in names if n not in data]
    if absent:
        raise ValueError("Saved original workspace missing: " + ", ".join(absent))
    return data


def source_portfolio(agrid, mu, params):
    """Literal source backsolve_egm.m gross k and foreign b for pre-shock states."""
    k0 = scalar(params, "k_initial")
    q0 = scalar(params, "q_initial")
    khat = scalar(params, "k_const_hat")
    rich = float(np.sum(mu[agrid > 0, :] * agrid[agrid > 0, None]))
    if rich <= 0:
        raise ValueError("No positive-asset wealth for original portfolio mapping")
    ell = (1 - khat) * k0 * q0 / rich
    k = np.maximum(ell * agrid + khat * k0, 0)
    b = agrid - k
    return k, b


def scalar_goods(x, *, a_pre, a_target, wage_labor, profits, gross_return,
                 x_pre, B_pre, y, investment, adj_cost, d_export,
                 theta_export, omega, theta):
    """Literal flexible-FX source c_t_guess(x) - net_cons(x), date 1."""
    m = 1 + x ** (theta - 1) * omega / (1 - omega)
    c_budget = ((a_pre + B_pre * x_pre / x) * gross_return
                - a_target + wage_labor + profits) / m
    net_cons = y - investment - adj_cost - d_export * x ** (-theta_export)
    return c_budget - net_cons


def positive_roots(fn, low=1e-3, high=1e2, samples=1601):
    """Only bracketed positive zeros; no numerical minimization is a zero."""
    xx = np.geomspace(low, high, samples)
    ff = [fn(float(x)) for x in xx]
    roots = []
    for i in range(samples - 1):
        a, b = ff[i], ff[i + 1]
        if np.isfinite(a) and a == 0:
            roots.append(float(xx[i]))
        if np.isfinite(a) and np.isfinite(b) and a * b < 0:
            roots.append(brentq(fn, float(xx[i]), float(xx[i + 1]),
                                xtol=1e-13, rtol=1e-13))
    return [x for i, x in enumerate(roots)
            if all(abs(x - y) > 1e-8 for y in roots[:i])]


def first_period_euler(grid, consumption, savings, mu, params, returns):
    """Separate, independent t=1 fine-grid Euler/KKT diagnostic.

    Source backsolve_egm uses future return r_t(t+1), Markov transition
    piex(current,next), and the stationary terminal continuation. This
    evaluates the stored fine policies without replacing off-grid values.
    """
    if consumption.ndim != 3 or consumption.shape[2] < 2:
        return {"euler_max": math.nan, "euler_undefined_mass": 1.0}
    c0, c1 = consumption[:, :, 0], consumption[:, :, 1]
    prob = np.asarray(parameter(params, "piex"), dtype=float)
    gamma = scalar(params, "gamma")
    beta = scalar(params, "beta")
    bmin = scalar(params, "bmin")
    if prob.shape != (mu.shape[1], mu.shape[1]):
        raise ValueError("Invalid source productivity transition dimensions")
    valid = (mu > 0) & (c0 > 0) & np.isfinite(savings)
    valid &= (savings >= grid[0]) & (savings <= grid[-1])
    rhs = np.zeros_like(c0)
    for j in range(mu.shape[1]):
        cnext = PchipInterpolator(grid, c1[:, j], extrapolate=False)(savings)
        valid &= (cnext > 0) & np.isfinite(cnext) | (mu == 0)
        good = np.isfinite(cnext) & (cnext > 0)
        next_mu = np.full_like(cnext, np.nan)
        next_mu[good] = cnext[good] ** (-gamma)
        rhs += prob[None, :, j] * np.nan_to_num(next_mu, nan=0.0)
    rhs *= beta * returns[1]
    lhs = np.zeros_like(c0)
    good_c = np.isfinite(c0) & (c0 > 0)
    lhs[good_c] = c0[good_c] ** (-gamma)
    valid &= lhs > 0
    # For a'=bmin, the inequality is u'(c_t) >= beta R E u'(c_{t+1}).
    binding = savings <= bmin + 1e-9
    discrepancy = np.zeros_like(lhs)
    discrepancy[valid & ~binding] = np.abs(
        1 - rhs[valid & ~binding] / lhs[valid & ~binding])
    discrepancy[valid & binding] = np.maximum(
        rhs[valid & binding] / lhs[valid & binding] - 1, 0)
    return {
        "euler_max": float(np.max(discrepancy[valid])) if np.any(valid) else math.nan,
        "euler_undefined_mass": float(np.sum(mu[(mu > 0) & ~valid])),
        "euler_mass_above_1e4": float(np.sum(mu[valid & (discrepancy > 1e-4)])),
    }


def audit(data):
    p = data["params"]
    grid = vector(parameter(p, "Agrid_fine"))
    pre = np.asarray(data["G0_shock"], dtype=float)
    post = np.asarray(data["In_dis_exp"], dtype=float)
    fine_c = np.asarray(data["c_tran_fine_exp"], dtype=float)
    fine_s = np.asarray(data["dec_tran_fine_exp"], dtype=float)
    if pre.shape != post.shape or pre.shape[0] != len(grid):
        raise ValueError("Original fine-grid/distribution dimensions inconsistent")
    if pre.shape[1] != 7 or len(grid) != 769:
        raise ValueError("Unexpected original 769x7 reference dimensions")
    if fine_c.shape[:2] != pre.shape or fine_s.shape[:2] != pre.shape:
        raise ValueError("Original fine-grid household policy dimensions inconsistent")
    if not np.all(np.isfinite(pre)) or not np.all(np.isfinite(post)):
        raise ValueError("Nonfinite original household distributions")
    if np.any(pre < -1e-14) or np.any(post < -1e-14):
        raise ValueError("Negative original probability mass")
    if abs(pre.sum() - 1) > 1e-7 or abs(post.sum() - 1) > 1e-7:
        raise ValueError("Original distribution mass does not sum to one")
    T = len(vector(data["tot_exp"]))
    kl = vector(data["kl_exp"])
    if T != 200 or len(kl) != 2 * T - 1:
        raise ValueError("Expected source 399-dimensional, T=200 solution")
    k0 = scalar(p, "k_initial")
    q0 = scalar(p, "q_initial")
    B0 = scalar(p, "BG_shock")
    xpre = scalar(p, "tot_initial")
    x1 = vector(data["tot_exp"])[0]
    K = np.r_[k0, kl[:T - 1]]
    L = kl[T - 1:]
    r = vector(data["r_t_exp"])
    gross = 1 + r[0]
    w = vector(data["w_nr_guess_exp"])[0]
    prof = vector(data["Profits_exp"])[0]
    y = vector(data["y_agg_exp"])[0]
    aprime = vector(data["new_shock_b_agg"])[0] + vector(data["q_agg_exp"])[0] * K[1]
    alpha = scalar(p, "alpha")
    delta = scalar(p, "delta")
    phi = scalar(p, "phi_ac")
    epsf = scalar(p, "epsilon_f")
    omega = scalar(p, "omega")
    theta = scalar(p, "theta")
    dex = scalar(p, "d_s")
    theta_star = scalar(p, "theta_star")
    I1 = K[1] - (1 - delta) * K[0]
    adj1 = phi / 2 * K[0] * ((I1 - delta * K[0]) / K[0]) ** 2
    z = vector(data["z_t"])[0]
    y_check = z * K[0] ** alpha * L[0] ** (1 - alpha)
    if abs(y - y_check) > 1e-7:
        raise ValueError("Source date-1 output inconsistent with K/L")
    wage_check = (epsf - 1) / epsf * (1 - alpha) * y / L[0]
    if abs(w - wage_check) > 1e-7:
        raise ValueError("Source date-1 wage inconsistent with flexible markup")
    m = 1 + x1 ** (theta - 1) * omega / (1 - omega)
    ex = vector(parameter(p, "ex"))
    exmean = scalar(p, "ex_mean")
    if ex.shape[0] != pre.shape[1]:
        raise ValueError("Source labor productivity dimensions inconsistent")
    labor_income = ((epsf - 1) / epsf * (1 - alpha)
                    * y * ex / exmean)
    k, b = source_portfolio(grid, pre, p)
    wealth_post = k + b * xpre / x1
    resources = (gross * wealth_post[:, None] + labor_income[None, :]
                 + prof - scalar(p, "bmin"))
    occupied = pre > 0
    failing = occupied & (resources <= 0)
    offgrid = occupied & ((wealth_post[:, None] < grid[0]) |
                           (wealth_post[:, None] > grid[-1]))
    c1 = fine_c[:, :, 0]
    s1 = fine_s[:, :, 0]
    post_occupied = post > 0
    source_resource = (gross * grid[:, None] + labor_income[None, :]
                       + prof - scalar(p, "bmin"))
    budget = m * c1 + s1 - gross * grid[:, None] - labor_income[None, :] - prof
    market_cha = float(np.sum(post * c1))
    netc = y - I1 - adj1 - dex * x1 ** (-theta_star)
    f = lambda x: scalar_goods(
        x, a_pre=q0*k0, a_target=aprime, wage_labor=w*L[0],
        profits=prof, gross_return=gross, x_pre=xpre, B_pre=B0,
        y=y, investment=I1, adj_cost=adj1, d_export=dex,
        theta_export=theta_star, omega=omega, theta=theta)
    f1 = f(x1)
    h = max(1e-8, x1 * 1e-6)
    fderiv = (f(x1 + h) - f(x1 - h)) / (2 * h)
    saved_cbudget = vector(data["c_t_agg_exp"])[0]
    saved_netc = vector(data["net_c_exp"])[0]
    if abs(saved_netc - netc) > 1e-7:
        raise ValueError("Source saved first-period net-consumption identity fails")
    if abs(saved_cbudget - (netc + f1)) > 1e-7:
        raise ValueError("Source saved first-period scalar-consumption identity fails")
    asset_agg = float(np.sum(post * s1))
    saved_resid = vector(data["resid"])
    good = post_occupied & np.isfinite(budget)
    if not np.all(np.isfinite(budget[post_occupied])):
        raise ValueError("Nonfinite first-period budget residual at occupied post-state")
    if not np.all(np.isfinite(c1[post_occupied])):
        raise ValueError("Nonfinite consumption at occupied original source post-state")
    euler = first_period_euler(grid, fine_c, s1, post, p, 1 + r)
    metrics = {
        "source_shock_terms_x1": float(x1),
        "source_pre_shock_terms": float(xpre),
        "source_nominal_depreciation_100log": float(100 * math.log(xpre/x1)),
        "source_outer_max_absolute_residual": float(np.max(np.abs(saved_resid))),
        "source_scalar_goods_residual_t1": float(f1),
        "source_scalar_derivative_t1": float(fderiv),
        "source_scalar_positive_roots": "|".join(f"{q:.15g}" for q in positive_roots(f)),
        "source_scalar_saved_C_gap": float(saved_cbudget - saved_netc),
        "source_household_realized_C_H_t1": market_cha,
        "source_household_C_H_minus_aggregate_C_H": market_cha - saved_cbudget,
        "source_realized_goods_residual_t1": float(netc - market_cha),
        "source_asset_market_residual_t1": asset_agg - aprime,
        "source_labor_market_residual_t1": float(saved_resid[T - 1]),
        "source_initial_gross_capital_gap": float(np.sum(pre*k[:, None]) - q0*k0),
        "source_initial_gross_bond_gap": float(np.sum(pre*b[:, None]) - B0),
        "source_initial_total_wealth_gap": float(np.sum(pre*grid[:, None]) - (q0*k0+B0)),
        "source_pre_interpolation_min_resources": float(np.min(resources[occupied])),
        "source_pre_interpolation_violating_states": int(np.sum(failing)),
        "source_pre_interpolation_violating_mass": float(np.sum(pre[failing])),
        "source_violating_net_wealth_pre_weighted": float(np.sum((pre*grid[:, None])[failing])),
        "source_violating_net_wealth_post_weighted": float(np.sum((pre*wealth_post[:, None])[failing])),
        "source_violating_gross_foreign_liability_weighted": float(np.sum((pre*np.maximum(-b,0)[:, None])[failing])),
        "source_original_state_offgrid_mass": float(np.sum(pre[offgrid])),
        "source_post_interpolation_min_resources": float(np.min(source_resource[post_occupied])),
        "source_post_interpolation_violating_mass": float(np.sum(post[post_occupied & (source_resource <= 0)])),
        "source_first_period_max_budget_abs": float(np.max(np.abs(budget[good]))),
        "source_first_period_weighted_budget_abs": float(np.sum(post * np.abs(budget))),
        "source_first_period_mass_with_nonpositive_c": float(np.sum(post[post_occupied & (c1 <= 0)])),
        "source_first_period_mass_at_c_floor": float(np.sum(post[post_occupied & np.isclose(c1, 1e-5, rtol=1e-8, atol=1e-10)])),
        "source_first_period_mass_borrowing_violation": float(np.sum(post[post_occupied & (s1 < scalar(p, "bmin") - 1e-9)])),
        "source_pre_distribution_mass": float(pre.sum()),
        "source_post_distribution_mass": float(post.sum()),
        **euler,
    }
    # Source path home-goods residuals are computed from actual aggregated
    # consumption c_t_exp, rather than the scalar c_t_agg_exp.
    Cpath = vector(data["c_t_exp"])
    Ypath = vector(data["y_agg_exp"])
    Xpath = vector(data["tot_exp"])
    if len(Cpath) < 30 or len(Ypath) < 30:
        raise ValueError("Missing source aggregate household paths")
    periods = []
    for t in range(30):
        invest = K[t+1] - (1 - delta)*K[t]
        adj = phi/2*K[t]*((invest-delta*K[t])/K[t])**2
        goods = Ypath[t] - invest - adj - dex*Xpath[t]**(-theta_star) - Cpath[t]
        periods.append((t+1, float(goods), float(saved_resid[t]),
                        float(saved_resid[T-1+t])))
    metrics["source_first30_max_realized_goods_abs"] = float(max(abs(a) for _,a,_,_ in periods))
    metrics["source_first30_max_asset_residual_abs"] = float(max(abs(a) for _,_,a,_ in periods))
    metrics["source_first30_max_labor_residual_abs"] = float(max(abs(a) for _,_,_,a in periods))
    return metrics, periods


def write_results(metrics, periods, folder):
    folder.mkdir(parents=True, exist_ok=True)
    with (folder/"first_period_metrics.csv").open("w", newline="", encoding="utf-8") as fp:
        w = csv.writer(fp)
        w.writerow(["metric", "value"])
        w.writerows(metrics.items())
    with (folder/"first_30_periods.csv").open("w", newline="", encoding="utf-8") as fp:
        w = csv.writer(fp)
        w.writerow(["period", "goods_residual", "source_asset_residual", "source_labor_residual"])
        w.writerows(periods)
    lines = ["# Stage 6 — Original authors' saved depreciation path",
             "", "Independent first-period economic diagnostics; **no GE solve**.",
             "The MATLAB source and its actual saved policies are used, not a conditional Julia branch.",
             "",
             "| Metric | Value |", "|---|---:|"]
    for key, value in metrics.items():
        if isinstance(value, float):
            shown = f"{value:.12g}"
        else:
            shown = str(value)
        lines.append(f"| {key} | {shown} |")
    lines += ["", "The t=1 pre-interpolation feasibility test applies to the "
              "original inherited (asset,productivity) states. Actual t=1 "
              "policies/budgets apply to the MATLAB post-revaluation fine grid.",
              "A nonzero scalar goods residual is not a full GE root; "
              "a failed conditional root does not imply global nonexistence.",
              "Euler diagnostics only cover the first date. Later Euler and "
              "infinite-horizon terminal conditions remain unverified.", ""]
    report = "\n".join(lines)
    (folder/"stage6_audit.md").write_text(report, encoding="utf-8")
    print(report)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", type=Path,
                    default=Path("results/author_reference/exp_flex_lev16.mat"))
    ap.add_argument("--output", type=Path,
                    default=Path("results/stage6_author_audit"))
    args = ap.parse_args()
    if not args.source.is_file():
        ap.error("Original MATLAB saved workspace is missing: " + str(args.source))
    metrics, periods = audit(read_original(args.source))
    write_results(metrics, periods, args.output)


if __name__ == "__main__":
    main()
