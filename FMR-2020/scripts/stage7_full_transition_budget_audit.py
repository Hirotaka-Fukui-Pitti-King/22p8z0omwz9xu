#!/usr/bin/env python3
"""Read-only all-interior-dates original MATLAB household-budget audit.

Independently transports the authors' exact 769x7 fine-grid distribution,
with the source's linear asset interpolation and exogenous Markov transition.
No EGM, prices, equilibrium, or input .mat values are recomputed.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np

from stage6_author_depreciation_audit import read_original, vector, parameter, scalar
from stage7_constrained_policy_audit import budget_consistent_binding


def transport_one(mu, savings, fine, transition):
    """Conservative lottery interpolation; never clamp positive-mass savings."""
    n,ns=mu.shape
    new=np.zeros_like(mu)
    for j in range(ns):
        valid=mu[:,j]>0
        if not np.any(valid):
            continue
        assets=savings[valid,j]
        masses=mu[valid,j]
        if np.any(~np.isfinite(assets)) or np.any(assets<fine[0]-1e-9) or np.any(assets>fine[-1]+1e-9):
            raise ValueError("Original positive-mass policy outside fine asset grid")
        a=np.maximum(fine[0],np.minimum(fine[-1],assets))
        top=np.searchsorted(fine,a,side="right")
        top=np.clip(top,1,n-1)
        bottom=top-1
        right_weight=(a-fine[bottom])/(fine[top]-fine[bottom])
        if np.any(right_weight < -1e-10) or np.any(right_weight > 1+1e-10):
            raise AssertionError("Invalid positive-mass interpolation weights")
        for jp in range(ns):
            prob=transition[j,jp]
            np.add.at(new[:,jp],bottom,masses*(1-right_weight)*prob)
            np.add.at(new[:,jp],top,masses*right_weight*prob)
    if np.any(new < -1e-12) or abs(new.sum()-mu.sum())>1e-10:
        raise AssertionError("Forward distribution fails conservation")
    return new


def audit_path(data):
    p=data["params"]
    fine=vector(parameter(p,"Agrid_fine"))
    mu=np.asarray(data["In_dis_exp"],dtype=float).copy()
    old_c=np.asarray(data["c_tran_fine_exp"],dtype=float)
    savings=np.asarray(data["dec_tran_fine_exp"],dtype=float)
    T=len(vector(data["tot_exp"]))
    if T!=200 or old_c.shape!=(769,7,T) or savings.shape!=(769,7,T-1):
        raise ValueError("Unexpected original saved fine-policy dimensions or horizon")
    transition=np.asarray(parameter(p,"piex"),dtype=float)
    if transition.shape!=(7,7) or np.max(abs(transition.sum(axis=1)-1))>1e-10:
        raise ValueError("Invalid saved exogenous transition")
    if abs(mu.sum()-1)>1e-8: raise ValueError("Invalid saved shock-date distribution")
    profits=vector(data["Profits_exp"])
    rates=1+vector(data["r_t_exp"])
    x=vector(data["tot_exp"])
    output=vector(data["y_agg_exp"])
    source_aggregates=vector(data["c_t_exp"])
    bmin=scalar(p,"bmin")
    markup=(scalar(p,"epsilon_f")-1)/scalar(p,"epsilon_f")
    wage_state=vector(parameter(p,"ex"))/scalar(p,"ex_mean")
    alpha=scalar(p,"alpha")
    theta=scalar(p,"theta")
    omega=scalar(p,"omega")
    records=[]
    max_mass=0
    for t in range(T-1):
        m=1+x[t]**(theta-1)*omega/(1-omega)
        income=markup*(1-alpha)*output[t]*wage_state+profits[t]
        cons=old_c[:,:,t]
        sav=savings[:,:,t]
        occupied=mu>0
        if not np.all(np.isfinite(cons[occupied])):
            raise ValueError(f"Positive-mass source consumption undefined at date {t+1}")
        old_budget=m*cons+sav-rates[t]*fine[:,None]-income[None,:]
        corrected,binding,_=budget_consistent_binding(
            cons,sav,fine,income,rates[t],bmin,m,mu)
        new_budget=m*corrected+sav-rates[t]*fine[:,None]-income[None,:]
        C_original=float(np.sum(mu[occupied]*cons[occupied]))
        C_corrected=float(np.sum(mu[occupied]*corrected[occupied]))
        C_saved=float(source_aggregates[t])
        if abs(C_original-C_saved)>1e-6:
            raise AssertionError(f"Independent source distribution fails original saved aggregate C at {t+1}: {C_original-C_saved}")
        errmass=float(np.sum(mu[occupied & (np.abs(old_budget)>1e-5)]))
        record={
            "t":t+1,
            "original_C_H":C_original,
            "corrected_fixed_path_C_H":C_corrected,
            "change_in_C_H":C_corrected-C_original,
            "mass_binding":float(np.sum(mu[binding])),
            "mass_old_budget_error_above_1e5":errmass,
            "max_old_budget_error":float(np.max(abs(old_budget[occupied]))),
            "max_corrected_budget_error":float(np.max(abs(new_budget[occupied]))),
            "weighted_old_budget_error":float(np.sum(mu[occupied]*old_budget[occupied])),
            "mass_error":float(abs(mu.sum()-1)),
        }
        records.append(record)
        max_mass=max(max_mass,record["mass_error"])
        if np.max(abs(new_budget[occupied]))>1e-8:
            raise AssertionError(f"Date {t+1} budget correction failed")
        if np.any(sav[occupied]<bmin-1e-9):
            raise AssertionError(f"Date {t+1} source policy below borrowing constraint")
        if t<T-2:
            mu=transport_one(mu,sav,fine,transition)
    return records


def save(records,out):
    out.mkdir(parents=True,exist_ok=True)
    keys=list(records[0])
    with (out/"stage7_full_transition_periods.csv").open("w",newline="",encoding="utf-8") as fp:
        writer=csv.DictWriter(fp,fieldnames=keys)
        writer.writeheader()
        writer.writerows(records)
    cols=["max_old_budget_error","max_corrected_budget_error",
          "mass_old_budget_error_above_1e5","mass_error","change_in_C_H"]
    metrics={key+"_max_abs":max(abs(row[key]) for row in records) for key in cols}
    metrics["interior_periods_checked"]=len(records)
    metrics["periods_with_positive_budget_error_mass"]=sum(row["mass_old_budget_error_above_1e5"]>0 for row in records)
    with (out/"stage7_full_transition_summary.csv").open("w",newline="",encoding="utf-8") as fp:
        writer=csv.writer(fp);writer.writerow(["metric","value"]);writer.writerows(metrics.items())
    print("STAGE 7 AUTHOR-SAVED FULL INTERIOR TRANSITION HOUSEHOLD AUDIT")
    for name,value in metrics.items(): print(f"{name}: {value:.15g}" if isinstance(value,float) else f"{name}: {value}")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--source",type=Path,default=Path("results/author_reference/exp_flex_lev16.mat"))
    parser.add_argument("--output",type=Path,default=Path("results/stage7_constrained_policy"))
    args=parser.parse_args()
    records=audit_path(read_original(args.source))
    save(records,args.output)


if __name__=="__main__":
    main()
