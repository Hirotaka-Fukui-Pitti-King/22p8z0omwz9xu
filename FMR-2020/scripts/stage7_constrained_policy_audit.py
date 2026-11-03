#!/usr/bin/env python3
"""Stage 7: original-source constrained-policy and aggregate budget reconciliation.

Independent read-only audit of the author's saved leveraged flexible-FX .mat.
Does not solve a macroeconomic equilibrium or modify the saved household policy.
The corrected-first-period counterfactual holds savings, macro prices, and
the post-revaluation distribution FIXED; it is not a new optimal policy/GE.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

from stage6_author_depreciation_audit import (
    read_original, vector, parameter, scalar, source_portfolio,
    first_period_euler,
)


def source_constrained_coarse(coarse, incomes, gross, bmin, m, floor=1e-5):
    """Exact original MATLAB lines 68-71, prior to coarse-to-fine PCHIP."""
    uncapped = (gross * coarse[:, None] + incomes[None, :] - bmin) / m
    return uncapped, np.maximum(uncapped, floor)


def budget_consistent_binding(consumption, savings, assets, income, gross, bmin, m, mass, tolerance=1e-9):
    """Correction ONLY at positive-mass binding states, with fixed savings.

    Nonbinding decisions and distribution remain unchanged. Reject economic
    infeasibility; never clamp a negative budget to a positive consumption.
    """
    old = np.asarray(consumption, dtype=float)
    saving = np.asarray(savings, dtype=float)
    mu = np.asarray(mass, dtype=float)
    binding = (mu > 0) & (np.abs(saving-bmin) <= tolerance)
    exact = (gross * assets[:, None] + income[None, :] - bmin) / m
    if np.any(~np.isfinite(exact[binding])) or np.any(exact[binding] <= 0):
        raise ValueError("Binding household with economically infeasible exact consumption")
    new = old.copy()
    new[binding] = exact[binding]
    return new, binding, exact


def audit_stage7(data):
    p = data["params"]
    fine = vector(parameter(p, "Agrid_fine"))
    coarse = vector(parameter(p, "Agrid"))
    pre = np.asarray(data["G0_shock"], dtype=float)
    mu = np.asarray(data["In_dis_exp"], dtype=float)
    old_policy = np.asarray(data["c_tran_fine_exp"], dtype=float)
    savings_policy = np.asarray(data["dec_tran_fine_exp"], dtype=float)
    if fine.size != 769 or coarse.size != 193 or pre.shape != (769,7) or mu.shape != pre.shape:
        raise ValueError("Unexpected actual author source grids or distributions")
    if np.any(np.diff(coarse) <= 0) or np.any(np.diff(fine) <= 0):
        raise ValueError("Original asset grids must be strictly ordered")
    if abs(mu.sum()-1)>1e-8 or abs(pre.sum()-1)>1e-8:
        raise ValueError("Original distributions must have unit mass")
    if np.any(mu < 0) or np.any(pre < 0):
        raise ValueError("Original distributions cannot have negative mass")
    c_old = old_policy[:,:,0]
    s = savings_policy[:,:,0]
    t = vector(data["tot_exp"])[0]
    tpre = scalar(p,"tot_initial")
    K0 = scalar(p,"k_initial")
    q0 = scalar(p,"q_initial")
    B0 = scalar(p,"BG_shock")
    bmin = scalar(p,"bmin")
    R = 1 + vector(data["r_t_exp"])[0]
    r_all = 1 + vector(data["r_t_exp"])
    w = vector(data["w_nr_guess_exp"])[0]
    L = vector(data["kl_exp"])[199]
    profits = vector(data["Profits_exp"])[0]
    income_ex = vector(parameter(p,"ex"))
    income = ((scalar(p,"epsilon_f")-1)/scalar(p,"epsilon_f")
              *(1-scalar(p,"alpha"))*vector(data["y_agg_exp"])[0]
              *income_ex/scalar(p,"ex_mean"))+profits
    m = 1 + t ** (scalar(p,"theta")-1) * scalar(p,"omega")/(1-scalar(p,"omega"))
    if not m > 0: raise ValueError("Nonpositive home-good expenditure conversion")

    coarse_raw, coarse_floor = source_constrained_coarse(coarse,income,R,bmin,m)
    interp_floor = np.column_stack([
        PchipInterpolator(coarse,coarse_floor[:,j],extrapolate=False)(fine)
        for j in range(7)])
    interp_raw = np.column_stack([
        PchipInterpolator(coarse,coarse_raw[:,j],extrapolate=False)(fine)
        for j in range(7)])
    c_new, binding, exact = budget_consistent_binding(c_old,s,fine,income,R,bmin,m,mu)
    occupied = mu > 0
    nonbinding = occupied & ~binding
    old_budget = m*c_old + s - R*fine[:,None]-income[None,:]
    new_budget = m*c_new + s - R*fine[:,None]-income[None,:]
    source_interp_discrepancy = c_old-interp_floor
    pchip_floor_contribution = interp_floor-interp_raw
    pchip_affine_contribution = interp_raw-exact
    coarse_floor_active = coarse_raw < 1e-5
    corrected_policies = old_policy.copy()
    corrected_policies[:,:,0] = c_new
    before_euler = first_period_euler(fine,old_policy,s,mu,p,r_all)
    after_euler = first_period_euler(fine,corrected_policies,s,mu,p,r_all)
    gross_capital,foreign_bonds=source_portfolio(fine,pre,p)
    ratio=tpre/t
    wealth_revalued=gross_capital+ratio*foreign_bonds
    pre_after_wealth=float(np.sum(pre*wealth_revalued[:,None]))
    post_wealth=float(np.sum(mu*fine[:,None]))
    assumed_wealth=q0*K0+B0*ratio
    total_savings=float(np.sum(mu*s))
    B_target=vector(data["new_shock_b_agg"])[0]
    q=vector(data["q_agg_exp"])[0]
    K_next=vector(data["kl_exp"])[0]
    asset_target=B_target+q*K_next
    aggregate_labor_income = float(np.sum(mu*(income[None,:]-profits)))
    macro_labor_income=w*L
    source_budget_consumption=float(vector(data["c_t_agg_exp"])[0])
    goods_consumption=float(vector(data["net_c_exp"])[0])
    C_old=float(np.sum(mu*c_old))
    C_new=float(np.sum(mu*c_new))
    aggregate_budget_error=float(np.sum(mu*old_budget))
    R_goods_old=goods_consumption-C_old
    R_goods_new=goods_consumption-C_new
    portfolio_wedge=R*(pre_after_wealth-assumed_wealth)/m
    transport_wedge=R*(post_wealth-pre_after_wealth)/m
    labor_wedge=(aggregate_labor_income-macro_labor_income)/m
    asset_wedge=-(total_savings-asset_target)/m
    household_budget_wedge=aggregate_budget_error/m
    accounting_total=(portfolio_wedge+transport_wedge+labor_wedge
                      +asset_wedge+household_budget_wedge)
    original_total=C_old-source_budget_consumption
    # The difference must close at machine precision, not merely approximately.
    if abs(accounting_total-original_total)>1e-8:
        raise AssertionError(f"Independent aggregate budget decomposition fails: {accounting_total-original_total}")
    delta_C=C_new-C_old
    if abs(delta_C+aggregate_budget_error/m)>1e-8:
        raise AssertionError("Binding correction does not match weighted budget residual")
    if np.max(np.abs(new_budget[occupied]))>1e-8:
        raise AssertionError("Counterfactual consumption does not restore household budgets")
    if np.max(np.abs(old_budget[nonbinding]))>1e-8:
        raise AssertionError("Unconstrained source household budget identity failed")
    if np.max(np.abs(source_interp_discrepancy[binding]))>1e-7:
        raise AssertionError("Source constrained consumption is not reconstructed by PCHIP")
    if np.max(np.abs(pchip_affine_contribution[binding]))>1e-8:
        raise AssertionError("PCHIP of original affine budget is not exact")
    if np.max(np.abs((source_interp_discrepancy+pchip_floor_contribution+
                      pchip_affine_contribution)[binding] -
                     (c_old-exact)[binding]))>1e-9:
        raise AssertionError("Constrained interpolation decomposition fails")
    if np.any(c_new[occupied] <= 0) or np.any(s[occupied]<bmin-1e-9):
        raise AssertionError("Correction lost household feasibility")

    idx=np.unravel_index(np.argmax(np.where(occupied,abs(old_budget),-1)),old_budget.shape)
    worst_ai,worst_j=idx
    below=np.searchsorted(coarse,fine[worst_ai],side="right")-1
    cutoff = int(np.count_nonzero(coarse_floor_active[:,worst_j]))
    old_max=float(np.max(np.abs(old_budget[occupied])))
    rows={
        "source_x1":float(t),"source_x_pre":float(tpre),
        "source_binding_mass":float(mu[binding].sum()),
        "source_budget_error_mass_over_1e5":float(mu[occupied & (abs(old_budget)>1e-5)].sum()),
        "source_max_absolute_budget_error":old_max,
        "corrected_max_absolute_budget_error":float(np.max(abs(new_budget[occupied]))),
        "source_weighted_signed_budget_error":aggregate_budget_error,
        "corrected_weighted_signed_budget_error":float(np.sum(mu*new_budget)),
        "original_aggregate_home_consumption":C_old,
        "corrected_fixed_path_home_consumption":C_new,
        "corrected_minus_original_consumption":delta_C,
        "source_scalar_budget_consumption":source_budget_consumption,
        "goods_available_home_consumption":goods_consumption,
        "original_goods_market_residual":R_goods_old,
        "corrected_fixed_path_goods_residual":R_goods_new,
        "goods_residual_change":R_goods_new-R_goods_old,
        "original_household_minus_scalar_budget_consumption":original_total,
        "corrected_household_minus_scalar_budget_consumption":C_new-source_budget_consumption,
        "decomposition_portfolio_revaluation":portfolio_wedge,
        "decomposition_distribution_transport":transport_wedge,
        "decomposition_labor_income":labor_wedge,
        "decomposition_asset_market":asset_wedge,
        "decomposition_individual_budget_error":household_budget_wedge,
        "decomposition_total":accounting_total,
        "decomposition_check_error":accounting_total-original_total,
        "initial_gross_capital_gap":float(np.sum(pre*gross_capital[:,None])-q0*K0),
        "initial_gross_foreign_bonds_gap":float(np.sum(pre*foreign_bonds[:,None])-B0),
        "first_period_saved_asset_market_gap":total_savings-asset_target,
        "distribution_mass":float(mu.sum()),
        "coarse_floor_active_grid_count_total":int(coarse_floor_active.sum()),
        "coarse_floor_active_grid_count_worst_income_state":cutoff,
        "binding_max_abs_original_minus_PCHIP_floored":float(np.max(abs(source_interp_discrepancy[binding]))),
        "binding_max_abs_PCHIP_unfloored_minus_exact":float(np.max(abs(pchip_affine_contribution[binding]))),
        "binding_max_abs_PCHIP_floor_wedge":float(np.max(abs(pchip_floor_contribution[binding]))),
        "binding_weighted_PCHIP_floor_contribution":float(np.sum(mu[binding]*pchip_floor_contribution[binding])),
        "binding_weighted_original_minus_PCHIP_floored":float(np.sum(mu[binding]*source_interp_discrepancy[binding])),
        "worst_asset_index_1based":int(worst_ai+1),
        "worst_productivity_index_1based":int(worst_j+1),
        "worst_asset":float(fine[worst_ai]),
        "worst_source_consumption":float(c_old[idx]),
        "worst_exact_binding_consumption":float(exact[idx]),
        "worst_source_budget_error":float(old_budget[idx]),
        "worst_near_coarse_index_1based":int(below+1),
        "worst_neighbor_coarse_floor_lower":bool(coarse_floor_active[max(below,0),worst_j]),
        "worst_neighbor_coarse_floor_upper":bool(coarse_floor_active[min(below+1,len(coarse)-1),worst_j]),
        "original_euler_max":before_euler["euler_max"],
        "corrected_first_period_fixed_savings_euler_max":after_euler["euler_max"],
        "original_euler_undefined_mass":before_euler["euler_undefined_mass"],
        "corrected_first_period_euler_undefined_mass":after_euler["euler_undefined_mass"],
    }
    binding_detail=[]
    for ai,j in zip(*np.where(occupied & binding)):
        binding_detail.append((
            int(ai+1),int(j+1),float(fine[ai]),float(mu[ai,j]),
            float(c_old[ai,j]),float(exact[ai,j]),float(c_new[ai,j]),
            float(interp_floor[ai,j]),float(interp_raw[ai,j]),
            float(s[ai,j]),float(old_budget[ai,j]),float(new_budget[ai,j]),
            float(pchip_floor_contribution[ai,j]),
            float(source_interp_discrepancy[ai,j]),
        ))
    return rows,binding_detail


def save(metrics,binding_rows,out):
    out.mkdir(parents=True,exist_ok=True)
    with (out/"stage7_metrics.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f); writer.writerow(["metric","value"]);writer.writerows(metrics.items())
    with (out/"stage7_binding_household_comparison.csv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.writer(f)
        writer.writerow(["asset_index_1based","income_index_1based","asset","mass",
                         "original_consumption","exact_budget_consumption",
                         "corrected_consumption","PCHIP_coarse_floored",
                         "PCHIP_coarse_unfloored","saved_savings",
                         "original_budget_error","corrected_budget_error",
                         "PCHIP_floor_wedge","source_minus_reconstructed_PCHIP"])
        writer.writerows(binding_rows)
    lines=["# Stage 7 — Original author first-date constrained-policy audit","",
           "Source .mat read only; fixed aggregate K/L, terms, wages, returns, "
           "profits, savings and initial/post-shock distribution. "
           "**No general-equilibrium re-solve or new Euler-optimality claim.**","",
           "| Diagnostic | Numerical result |","|---|---:|"]
    for name,value in metrics.items():
        show=f"{value:.13g}" if isinstance(value,float) else str(value)
        lines.append(f"| {name} | {show} |")
    lines+=["","The corrected counterfactual enforces the fine-node budget "
             "identity ONLY on occupied binding states, leaving source savings "
             "and all aggregates fixed. PCHIP floor and portfolio wedges are "
             "measured separately. Corrected first-period Euler numbers are "
             "diagnostics, not a new optimal policy.",""]
    report="\n".join(lines)
    (out/"stage7_audit.md").write_text(report,encoding="utf-8")
    print(report)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,
                        default=Path("results/author_reference/exp_flex_lev16.mat"))
    parser.add_argument("--output",type=Path,
                        default=Path("results/stage7_constrained_policy"))
    args=parser.parse_args()
    if not args.source.is_file():
        parser.error("Missing original author .mat: "+str(args.source))
    m,rows=audit_stage7(read_original(args.source))
    save(m,rows,args.output)


if __name__=="__main__":
    main()
