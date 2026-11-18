#!/usr/bin/env python3
"""
Exact Forensic Biostatistical Computation Script for Track 1
Acinetobacter baumannii Preclinical Murine Subunit & Multi-Epitope Vaccine Meta-Analysis
(k = 22 independent experiments across 20 publications, 48 active vaccine arms, N = 574 mice)

Strict rules:
- All counts, routes, doses, strains, SYRCLE items, and exclusions strictly reflect full-text audit records.
- DerSimonian-Laird (DL) random-effects estimator on ln(RR) with +0.5 Haldane-Anscombe correction on 2x2 tables containing a zero cell.
- Hartung-Knapp-Sidik-Jonkman (HKSJ) variance adjustment with max(1, q) rule and t_{k-1}.
- 95% Prediction Intervals using t_{k-2}.
- Additive Risk Difference (RD): raw proportions for point estimates; Haldane-adjusted proportions (+0.5 / +1) for sampling variance of zero-cell tables.
- Fixed-effect Inverse-Variance RR and Mantel-Haenszel RR (uncorrected exact & continuity-corrected).
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from statistics import NormalDist
from typing import Any

import numpy as np

_NORM = NormalDist(0.0, 1.0)


def norm_cdf(x: float) -> float:
    return _NORM.cdf(x)


def norm_ppf(p: float) -> float:
    return _NORM.inv_cdf(p)


def _reg_inc_gamma_p(a: float, x: float) -> float:
    """Exact regularized lower incomplete gamma P(a, x) to 1e-15 precision."""
    if x <= 0.0:
        return 0.0
    gln = math.lgamma(a)
    if x < a + 1.0:
        ap = a
        s = 1.0 / a
        del_v = s
        for _ in range(1000):
            ap += 1.0
            del_v *= x / ap
            s += del_v
            if abs(del_v) < abs(s) * 1e-15:
                break
        return s * math.exp(-x + a * math.log(x) - gln)
    else:
        # Modified Lentz continued fraction for Q(a, x) = 1 - P(a, x)
        tiny = 1e-300
        b = x + 1.0 - a
        c = 1.0 / tiny
        d = 1.0 / b
        h = d
        for i in range(1, 1000):
            an = -i * (i - a)
            b += 2.0
            d = an * d + b
            if abs(d) < tiny:
                d = tiny
            c = b + an / c
            if abs(c) < tiny:
                c = tiny
            d = 1.0 / d
            del_v = d * c
            h *= del_v
            if abs(del_v - 1.0) < 1e-15:
                break
        q = math.exp(-x + a * math.log(x) - gln) * h
        return 1.0 - q


def chi2_cdf(x: float, df: int) -> float:
    if x <= 0.0:
        return 0.0
    return _reg_inc_gamma_p(0.5 * df, 0.5 * x)


def _betacf(a: float, b: float, x: float) -> float:
    tiny = 1e-300
    qab = a + b
    qap = a + 1.0
    qam = a - 1.0
    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 1000):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < tiny:
            d = tiny
        c = 1.0 + aa / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        del_v = d * c
        h *= del_v
        if abs(del_v - 1.0) < 1e-15:
            break
    return h


def reg_inc_beta(a: float, b: float, x: float) -> float:
    """Exact regularized incomplete beta I_x(a, b) to 1e-15 precision."""
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    bt = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log(1.0 - x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    return 1.0 - bt * _betacf(b, a, 1.0 - x) / b


def t_cdf(t_val: float, df: int) -> float:
    x = df / (df + t_val * t_val)
    ib = 0.5 * reg_inc_beta(0.5 * df, 0.5, x)
    return 1.0 - ib if t_val >= 0.0 else ib


def t_ppf(p: float, df: int) -> float:
    """Exact inverse Student-t CDF via bracketed bisection + Halley/Newton."""
    if p == 0.5:
        return 0.0
    if p < 0.5:
        return -t_ppf(1.0 - p, df)
    lo, hi = 0.0, 1.0
    while t_cdf(hi, df) < p:
        hi *= 2.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if t_cdf(mid, df) < p:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-14 * max(1.0, mid):
            break
    return 0.5 * (lo + hi)


def linregress_exact(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    n = len(x)
    mx = float(np.mean(x))
    my = float(np.mean(y))
    ssxx = float(np.sum((x - mx) ** 2))
    ssxy = float(np.sum((x - mx) * (y - my)))
    slope = ssxy / ssxx
    intercept = my - slope * mx
    resid = y - (intercept + slope * x)
    s2 = float(np.sum(resid ** 2) / (n - 2))
    se_slope = math.sqrt(s2 / ssxx)
    se_intercept = math.sqrt(s2 * (1.0 / n + (mx ** 2) / ssxx))
    return {
        "slope": slope,
        "intercept": intercept,
        "stderr": se_slope,
        "intercept_stderr": se_intercept,
    }


def fisher_exact_two_sided(table: list[list[int]]) -> tuple[None, float]:
    """Exact two-sided Fisher's exact test (sum of hypergeometric probs <= observed)."""
    a, b = table[0]
    c, d = table[1]
    r1 = a + b
    r2 = c + d
    c1 = a + c
    n = r1 + r2
    denom = math.comb(n, c1)

    def prob(k: int) -> float:
        return (math.comb(r1, k) * math.comb(r2, c1 - k)) / denom

    p_obs = prob(a)
    k_min = max(0, c1 - r2)
    k_max = min(r1, c1)
    p_sum = sum(prob(k) for k in range(k_min, k_max + 1) if prob(k) <= p_obs * (1.0 + 1e-12))
    return None, min(1.0, p_sum)


# Self-verify special functions against reference constants (1e-12 tolerance)
assert abs(t_ppf(0.975, 21) - 2.07961384472768) < 1e-12
assert abs(t_ppf(0.975, 20) - 2.08596344726586) < 1e-12
assert abs(t_ppf(0.975, 1) - math.tan(0.475 * math.pi)) < 1e-12
assert abs((1.0 - chi2_cdf(4.0, 2)) - math.exp(-2.0)) < 1e-14
assert abs((1.0 - chi2_cdf(3.841458820694124, 1)) - 0.05) < 1e-13

OUT_DIR = Path(r"d:\Graduation\systematic_review\forensic_audit_workspace\track1")
OUT_DIR.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 1. VERIFIED ARM-LEVEL & EXPERIMENT-LEVEL DATASET (22 EXPERIMENTS, 48 ARMS)
# =============================================================================

EXPERIMENTS_META: list[dict[str, Any]] = [
    {
        "exp_id": "Badmasti_2015",
        "publication": "Badmasti 2015",
        "citation_short": "Badmasti et al., 2015",
        "doi": "10.1016/j.molimm.2015.07.031",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "OmpA; Other single; Multivalent",
        "antigens_exact": "AbOmpA(8-346aa), Bap(1-487aa), or AbOmpA(8-346aa) + Bap(1-487aa) [OMV(PagL)+Bap arm excluded as OMV]",
        "vacc_route": "NS in audit excerpt",
        "vacc_dose_schedule": "NS in audit excerpt",
        "adjuvant_exact": "Alum",
        "adjuvant_subgroup": "Alum",
        "mouse_strain": "C57BL/6",
        "mouse_sex": "Female",
        "mouse_age": "6-8 wk",
        "challenge_strain": "ATCC 19606",
        "challenge_route": "IP",
        "challenge_inoculum": "1.14 x 10^5 CFU (+ porcine mucin)",
        "followup_duration": "7 days",
        "ctrl_events": 0,
        "ctrl_n": 7,
        "count_reporting_basis": "Derived from stated death counts + group n=7 (AbOmpA: 2 deaths; Bap: 4 deaths; AbOmpA+Bap: 0 deaths vs ATCC 19606; controls all died within 1 d)",
        "inferred_control_flag": False,
        "caveats_notes": "Primary ATCC 19606 sepsis cohort used (all 3 subunit arms tested concurrently); OMV(PagL)+Bap arm excluded as OMV platform; separate clinical isolate AB-44 cohort tested subset arms only (1 death in AbOmpA+Bap).",
        "arms": [
            {"arm_name": "AbOmpA(8-346aa)", "platform": "OmpA", "ev": 5, "nv": 7, "basis": "2 deaths by day 2 of n=7"},
            {"arm_name": "Bap(1-487aa)", "platform": "Other single", "ev": 3, "nv": 7, "basis": "4 deaths by day 4 of n=7"},
            {"arm_name": "AbOmpA + Bap", "platform": "Multivalent/chimeric", "ev": 7, "nv": 7, "basis": "0 deaths vs ATCC 19606 of n=7"},
        ],
    },
    {
        "exp_id": "Singh_2017",
        "publication": "Singh 2017",
        "citation_short": "Singh et al., 2017",
        "doi": "10.1038/s41598-017-12789-3",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "BamA / Oma87",
        "antigens_exact": "rBamA",
        "vacc_route": "IP",
        "vacc_dose_schedule": "20 ug (+ 2% Al(OH)3)",
        "adjuvant_exact": "2% Al(OH)3 (Alum)",
        "adjuvant_subgroup": "Alum",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Female",
        "mouse_age": "6-8 wk",
        "challenge_strain": "MDR clinical strain P-562",
        "challenge_route": "IN",
        "challenge_inoculum": "10^9 CFU",
        "followup_duration": "7 days",
        "ctrl_events": 0,
        "ctrl_n": 10,
        "count_reporting_basis": "Derived from stated 80% survival + n=10/group; all controls died within 24 h",
        "inferred_control_flag": False,
        "caveats_notes": "IP active immunization followed by IN pulmonary challenge; separate passive serum transfer and organ CFU cohorts not synthesized.",
        "arms": [
            {"arm_name": "rBamA", "platform": "BamA / Oma87", "ev": 8, "nv": 10, "basis": "80% of n=10"},
        ],
    },
    {
        "exp_id": "Rasooli_2020",
        "publication": "Rasooli 2020",
        "citation_short": "Rasooli et al., 2020",
        "doi": "10.1007/s10989-020-10056-0",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "BamA / Oma87",
        "antigens_exact": "rOma87 (BamA ortholog)",
        "vacc_route": "SC",
        "vacc_dose_schedule": "SC + CFA/IFA",
        "adjuvant_exact": "Freund's (CFA/IFA)",
        "adjuvant_subgroup": "Freund's",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Female",
        "mouse_age": "6-8 wk (22-25 g)",
        "challenge_strain": "ATCC 19606",
        "challenge_route": "IP",
        "challenge_inoculum": "~2xLD (+ mucin)",
        "followup_duration": "4 days (96 h)",
        "ctrl_events": 0,
        "ctrl_n": 6,
        "count_reporting_basis": "Derived/inferred from overall stated n=8 test, n=6 control + prose 'all vaccinated survived' and 'all controls died within 24 h' at ~2xLD",
        "inferred_control_flag": True,
        "caveats_notes": "Overall n=8 test and n=6 control stated without explicit per-challenge-dose breakdown (~7xLD also tested with 25% vaccinated survival); 96-h (4-day) follow-up; Draft Table 1 wrote '2 x LD100' whereas full text states ~2xLD / 2xLD50. Tested by omission in sensitivity analysis.",
        "arms": [
            {"arm_name": "rOma87 (~2xLD)", "platform": "BamA / Oma87", "ev": 8, "nv": 8, "basis": "All vaccinated survived, overall n=8"},
        ],
    },
    {
        "exp_id": "Araujo_Conde_2021",
        "publication": "Araujo/Conde 2021",
        "citation_short": "Araujo et al., 2021 (cited as Conde et al., 2021 in draft)",
        "doi": "10.1016/j.micinf.2021.104801",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "BamA / Oma87",
        "antigens_exact": "Insoluble rAbBamA inclusion bodies",
        "vacc_route": "IM",
        "vacc_dose_schedule": "25 ug + Alhydrogel",
        "adjuvant_exact": "Alhydrogel (Alum)",
        "adjuvant_subgroup": "Alum",
        "mouse_strain": "C57BL/6",
        "mouse_sex": "Female",
        "mouse_age": "8 wk",
        "challenge_strain": "Clinical isolate AB162",
        "challenge_route": "IP",
        "challenge_inoculum": "1 x 10^4 CFU in 2.5% mucin",
        "followup_duration": "7 days",
        "ctrl_events": 3,
        "ctrl_n": 10,
        "count_reporting_basis": "Derived from stated 70% vaccine vs 30% control survival ('40% increase in survival') + n=10/group",
        "inferred_control_flag": False,
        "caveats_notes": "First author is Araujo AEV (Conde L is second author); non-zero control survival (3/10 = 30.0%); Draft Section 4.1.1 erroneously states '50% survival' instead of 70%.",
        "arms": [
            {"arm_name": "Insoluble rAbBamA", "platform": "BamA / Oma87", "ev": 7, "nv": 10, "basis": "70% of n=10"},
        ],
    },
    {
        "exp_id": "Mansouri_2023",
        "publication": "Mansouri 2023",
        "citation_short": "Mansouri et al., 2023",
        "doi": "10.1016/j.intimp.2023.110650",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "BamA; Other single; Multivalent",
        "antigens_exact": "Oma87, Bap, or Oma87 + Bap",
        "vacc_route": "SC",
        "vacc_dose_schedule": "SC + CFA/IFA",
        "adjuvant_exact": "Freund's (CFA/IFA)",
        "adjuvant_subgroup": "Freund's",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Female",
        "mouse_age": "NS (20-25 g)",
        "challenge_strain": "ATCC 19606",
        "challenge_route": "IP",
        "challenge_inoculum": "Sub-LD100 sepsis challenge (50% control survival)",
        "followup_duration": "4 days (96 h)",
        "ctrl_events": 7,
        "ctrl_n": 14,
        "count_reporting_basis": "Derived from stated 100% survival in each active arm vs 50% shared control + 4 groups of n=14",
        "inferred_control_flag": False,
        "caveats_notes": "50% control survival (7/14) creates ceiling effect on RR (~2.0 max) and carries the single largest inverse-variance weight (46.8%) on lnRR; 96-h follow-up.",
        "arms": [
            {"arm_name": "Oma87", "platform": "BamA / Oma87", "ev": 14, "nv": 14, "basis": "100% of n=14"},
            {"arm_name": "Bap", "platform": "Other single", "ev": 14, "nv": 14, "basis": "100% of n=14"},
            {"arm_name": "Oma87 + Bap", "platform": "Multivalent/chimeric", "ev": 14, "nv": 14, "basis": "100% of n=14"},
        ],
    },
    {
        "exp_id": "Huang_2015",
        "publication": "Huang 2015",
        "citation_short": "Huang et al., 2015",
        "doi": "10.1016/j.vaccine.2015.07.031",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "OmpW",
        "antigens_exact": "Trx-OmpW",
        "vacc_route": "SC",
        "vacc_dose_schedule": "50 ug + Alum",
        "adjuvant_exact": "Alum",
        "adjuvant_subgroup": "Alum",
        "mouse_strain": "ICR",
        "mouse_sex": "Female",
        "mouse_age": "6-8 wk",
        "challenge_strain": "Clinical isolate Ab1",
        "challenge_route": "IP",
        "challenge_inoculum": "1 x 10^6 CFU (+ mucin)",
        "followup_duration": "7 days",
        "ctrl_events": 1,
        "ctrl_n": 6,
        "count_reporting_basis": "Derived from stated 100% vs 16.7% survival + n=6/group",
        "inferred_control_flag": False,
        "caveats_notes": "Non-zero control survival (1/6 = 16.7%); zero non-event cell in vaccine group (6/6) triggers Haldane +0.5.",
        "arms": [
            {"arm_name": "Trx-OmpW", "platform": "OmpW", "ev": 6, "nv": 6, "basis": "100% of n=6"},
        ],
    },
    {
        "exp_id": "Abdollahi_2023",
        "publication": "Abdollahi 2023",
        "citation_short": "Abdollahi & Raoufi, 2023",
        "doi": "10.1016/j.micpath.2023.106317",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "OmpW",
        "antigens_exact": "Whole OmpW2 or OmpW2 epitopes 1+3",
        "vacc_route": "SC",
        "vacc_dose_schedule": "10 ug + CFA/IFA",
        "adjuvant_exact": "Freund's (CFA/IFA)",
        "adjuvant_subgroup": "Freund's",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Female",
        "mouse_age": "6-8 wk",
        "challenge_strain": "A. baumannii (IN LD100)",
        "challenge_route": "IN",
        "challenge_inoculum": "LD100",
        "followup_duration": "NS",
        "ctrl_events": 0,
        "ctrl_n": 8,
        "count_reporting_basis": "Derived/inferred from challenge section n=8/group + 'completely protective' (100%) + control 0/8 inferred from LD100 definition",
        "inferred_control_flag": True,
        "caveats_notes": "Methods internal discrepancy: immunization §2.16 states n=5 whereas challenge §2.17.2 states n=8/group; control 0/8 inferred from LD100; follow-up duration NS. Tested by omission in sensitivity analysis.",
        "arms": [
            {"arm_name": "OmpW2 whole protein", "platform": "OmpW", "ev": 8, "nv": 8, "basis": "Completely protective, challenge n=8"},
            {"arm_name": "OmpW2 epitopes 1+3", "platform": "OmpW", "ev": 8, "nv": 8, "basis": "Completely protective, challenge n=8"},
        ],
    },
    {
        "exp_id": "Huang_2016",
        "publication": "Huang 2016",
        "citation_short": "Huang et al., 2016",
        "doi": "10.1038/srep20724",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "Other single",
        "antigens_exact": "Trx-Omp22 (5, 10, 20, and 50 ug dose arms)",
        "vacc_route": "SC",
        "vacc_dose_schedule": "5, 10, 20, or 50 ug + Alum",
        "adjuvant_exact": "Alum",
        "adjuvant_subgroup": "Alum",
        "mouse_strain": "ICR",
        "mouse_sex": "Female",
        "mouse_age": "6-8 wk",
        "challenge_strain": "Clinical isolate Ab1",
        "challenge_route": "IP",
        "challenge_inoculum": "1 x 10^6 CFU",
        "followup_duration": "7 days",
        "ctrl_events": 0,
        "ctrl_n": 6,
        "count_reporting_basis": "Derived from stated 100% (50 ug), 33% (20 ug), 33% (10 ug), 0% (5 ug) vs 0% control + n=6/group",
        "inferred_control_flag": False,
        "caveats_notes": "All 4 dose-ranging arms sharing 1 control group pooled per protocol (6/6 + 2/6 + 2/6 + 0/6 = 10/24); optimal 50 ug dose alone achieved 6/6 (100%). Distinct from Huang et al., 2016 Sci Rep 6:37242 (E. coli OMV-Omp22, excluded).",
        "arms": [
            {"arm_name": "Trx-Omp22 (50 ug)", "platform": "Other single", "ev": 6, "nv": 6, "basis": "100% of n=6"},
            {"arm_name": "Trx-Omp22 (20 ug)", "platform": "Other single", "ev": 2, "nv": 6, "basis": "33% of n=6"},
            {"arm_name": "Trx-Omp22 (10 ug)", "platform": "Other single", "ev": 2, "nv": 6, "basis": "33% of n=6"},
            {"arm_name": "Trx-Omp22 (5 ug)", "platform": "Other single", "ev": 0, "nv": 6, "basis": "0% of n=6"},
        ],
    },
    {
        "exp_id": "Akbari_2022",
        "publication": "Akbari 2022",
        "citation_short": "Akbari et al., 2022",
        "doi": "10.1016/j.intimp.2022.108731",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "Other single; Multivalent",
        "antigens_exact": "LCL-scaffolded BauA loop 7 (BauAL7P3), Omp34 loop 3 (Omp34L3P1), or hybrid BauAL7P3Omp34L3P1",
        "vacc_route": "SC",
        "vacc_dose_schedule": "20 ug SC on days 1, 14, 28 + CFA/IFA",
        "adjuvant_exact": "Freund's (CFA/IFA)",
        "adjuvant_subgroup": "Freund's",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Male",
        "mouse_age": "4-6 wk (20-25 g)",
        "challenge_strain": "MDR clinical isolate ABI022",
        "challenge_route": "IP",
        "challenge_inoculum": "4.76 x 10^4 CFU",
        "followup_duration": "7 days (168 h)",
        "ctrl_events": 0,
        "ctrl_n": 7,
        "count_reporting_basis": "Derived from exact percentages 42.86%, 42.86%, 71.43% vs 0% + n=7/group",
        "inferred_control_flag": False,
        "caveats_notes": "Also appraised in Track 2 (score 7/14); passive sera survival (66.67% combined, 33.33% single) reported separately and excluded from active synthesis.",
        "arms": [
            {"arm_name": "LCL BauA loop 7", "platform": "Other single", "ev": 3, "nv": 7, "basis": "42.86% of n=7"},
            {"arm_name": "LCL Omp34 loop 3", "platform": "Other single", "ev": 3, "nv": 7, "basis": "42.86% of n=7"},
            {"arm_name": "LCL BauA loop 7 + Omp34 loop 3", "platform": "Multivalent/chimeric", "ev": 5, "nv": 7, "basis": "71.43% of n=7"},
        ],
    },
    {
        "exp_id": "Chaudhuri_2022",
        "publication": "Chaudhuri 2022",
        "citation_short": "Chaudhuri et al., 2022",
        "doi": "10.3389/fimmu.2022.933445",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "Other single; Multivalent",
        "antigens_exact": "rBauA and 6 LCL-scaffolded BauA surface-loop constructs (loop 5, loop 7, loop 8, loops 7+5, loops 8+5, loops 8+7+5)",
        "vacc_route": "SC",
        "vacc_dose_schedule": "SC + CFA/IFA",
        "adjuvant_exact": "Freund's (CFA/IFA)",
        "adjuvant_subgroup": "Freund's",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Male",
        "mouse_age": "6-8 wk (20-25 g)",
        "challenge_strain": "ATCC 19606",
        "challenge_route": "IP",
        "challenge_inoculum": "1.8 x 10^6 CFU",
        "followup_duration": "7 days (168 h)",
        "ctrl_events": 0,
        "ctrl_n": 6,
        "count_reporting_basis": "Derived from stated group outcomes/percentages + n=6/group",
        "inferred_control_flag": False,
        "caveats_notes": "Multi-loop constructs (7+5, 8+5, 8+7+5) combine multiple surface loops from a single antigen (BauA) on the LCL beta-barrel scaffold; none exceeded whole rBauA (6/6) or single LCL loop 7 (6/6).",
        "arms": [
            {"arm_name": "rBauA", "platform": "Other single", "ev": 6, "nv": 6, "basis": "100% of n=6"},
            {"arm_name": "LCL BauA loop 5", "platform": "Other single", "ev": 4, "nv": 6, "basis": "66.7% of n=6"},
            {"arm_name": "LCL BauA loop 7", "platform": "Other single", "ev": 6, "nv": 6, "basis": "100% of n=6"},
            {"arm_name": "LCL BauA loop 8", "platform": "Other single", "ev": 3, "nv": 6, "basis": "50.0% of n=6"},
            {"arm_name": "LCL BauA loops 7+5", "platform": "Multivalent/chimeric", "ev": 5, "nv": 6, "basis": "83.3% of n=6"},
            {"arm_name": "LCL BauA loops 8+5", "platform": "Multivalent/chimeric", "ev": 4, "nv": 6, "basis": "66.7% of n=6"},
            {"arm_name": "LCL BauA loops 8+7+5", "platform": "Multivalent/chimeric", "ev": 5, "nv": 6, "basis": "83.3% of n=6"},
        ],
    },
    {
        "exp_id": "Qamsari_2020",
        "publication": "Qamsari 2020",
        "citation_short": "Qamsari et al., 2020",
        "doi": "10.3389/fimmu.2020.00158",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "Other single; Multivalent",
        "antigens_exact": "rZnuD and 5 LCL-scaffolded ZnuD surface-loop constructs (loop 2, loop 5, loop 7, loop 11, 4-loop hybrid 2+5+7+11)",
        "vacc_route": "SC",
        "vacc_dose_schedule": "SC + CFA/IFA",
        "adjuvant_exact": "Freund's (CFA/IFA)",
        "adjuvant_subgroup": "Freund's",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Male",
        "mouse_age": "6-8 wk (20-25 g)",
        "challenge_strain": "ATCC 19606",
        "challenge_route": "IP",
        "challenge_inoculum": "1.8 x 10^6 CFU",
        "followup_duration": "7 days (168 h)",
        "ctrl_events": 0,
        "ctrl_n": 8,
        "count_reporting_basis": "Explicit x/n literally printed in paper (5/8, 3/8, 2/8, 4/8, 2/6, 8/8 vs 0/8)",
        "inferred_control_flag": False,
        "caveats_notes": "Loop 11 arm had n=6 instead of n=8 because 2 mice were 'inadvertently removed prior to challenge'. Draft Section 3.3.2 wrote '2/8-5/8' without noting loop 11 was 2/6 (33.3%).",
        "arms": [
            {"arm_name": "rZnuD", "platform": "Other single", "ev": 5, "nv": 8, "basis": "Printed 5/8"},
            {"arm_name": "LCL ZnuD loop 2", "platform": "Other single", "ev": 3, "nv": 8, "basis": "Printed 3/8"},
            {"arm_name": "LCL ZnuD loop 5", "platform": "Other single", "ev": 2, "nv": 8, "basis": "Printed 2/8"},
            {"arm_name": "LCL ZnuD loop 7", "platform": "Other single", "ev": 4, "nv": 8, "basis": "Printed 4/8"},
            {"arm_name": "LCL ZnuD loop 11", "platform": "Other single", "ev": 2, "nv": 6, "basis": "Printed 2/6 (2 mice removed pre-challenge)"},
            {"arm_name": "LCL ZnuD 4-loop hybrid", "platform": "Multivalent/chimeric", "ev": 8, "nv": 8, "basis": "Printed 8/8"},
        ],
    },
    {
        "exp_id": "Ren_2019",
        "publication": "Ren 2019",
        "citation_short": "Ren et al., 2019",
        "doi": "10.4414/smw.2019.20052",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "Multivalent",
        "antigens_exact": "rMEP: Ata(aa 420-621) backbone + 4 B-cell and 2 T-cell epitopes from FilF and NucAb (GGGGSGGGG linker; NO OmpA)",
        "vacc_route": "SC",
        "vacc_dose_schedule": "30 ug SC on days 0, 14, 21 + CFA/IFA",
        "adjuvant_exact": "Freund's (CFA prime, IFA boosts; NO Alum)",
        "adjuvant_subgroup": "Freund's",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Female",
        "mouse_age": "6-8 wk",
        "challenge_strain": "ATCC 19606",
        "challenge_route": "IP",
        "challenge_inoculum": "2 x 10^8 CFU",
        "followup_duration": "10 days",
        "ctrl_events": 0,
        "ctrl_n": 9,
        "count_reporting_basis": "Derived from stated 88.9% survival + n=9/group; all controls died within 72 h",
        "inferred_control_flag": False,
        "caveats_notes": "CRITICAL DRAFT ERROR: Draft Table 1 and Table 3 falsely list OmpA ('FilF, OmpA, NucAb') and Alum; full text proves rMEP is Ata(420-621) + FilF + NucAb with Freund's only. Also in Track 2; 72 mice procured vs 63 allocated unexplained.",
        "arms": [
            {"arm_name": "rMEP (Ata + FilF + NucAb)", "platform": "Multivalent/chimeric", "ev": 8, "nv": 9, "basis": "88.9% of n=9"},
        ],
    },
    {
        "exp_id": "Hagag_2022",
        "publication": "Hagag 2022",
        "citation_short": "Hagag et al., 2022",
        "doi": "10.1007/s00253-022-12231-3",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "Other single; Multivalent",
        "antigens_exact": "Wza, YiaD, or Wza + YiaD protein arms [unconjugated CPS and CPS-containing pentavalent arms excluded]",
        "vacc_route": "SC",
        "vacc_dose_schedule": "SC + Alum",
        "adjuvant_exact": "Alum",
        "adjuvant_subgroup": "Alum",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Female",
        "mouse_age": "NS (20-25 g)",
        "challenge_strain": "ATCC 19606",
        "challenge_route": "IP",
        "challenge_inoculum": "5 x 10^8 CFU (+ mucin)",
        "followup_duration": "7 days",
        "ctrl_events": 2,
        "ctrl_n": 10,
        "count_reporting_basis": "Derived from stated 70% (Wza), 60% (YiaD), 100% (Wza+YiaD) vs 20% control + n=10/group",
        "inferred_control_flag": False,
        "caveats_notes": "Unconjugated CPS and pentavalent (containing unconjugated CPS) arms excluded per protocol; non-zero control survival (2/10 = 20%); no zero cells at pooled experiment level.",
        "arms": [
            {"arm_name": "Wza", "platform": "Other single", "ev": 7, "nv": 10, "basis": "70% of n=10"},
            {"arm_name": "YiaD", "platform": "Other single", "ev": 6, "nv": 10, "basis": "60% of n=10"},
            {"arm_name": "Wza + YiaD", "platform": "Multivalent/chimeric", "ev": 10, "nv": 10, "basis": "100% of n=10"},
        ],
    },
    {
        "exp_id": "Prather_2026_exp1",
        "publication": "Prather 2026",
        "citation_short": "Prather et al., 2026, exp. 1",
        "doi": "10.3390/vaccines14030275",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "Multivalent",
        "antigens_exact": "rAMEV1 (TrxA + pOmpA, pPlc1, pBamA, pBauA, pBlp2) or rAMEV2 (TrxA + pNlpE, pNucAB, pOmp38, pTonB, pZnuD)",
        "vacc_route": "SC",
        "vacc_dose_schedule": "3 doses SC + TiterMax Gold",
        "adjuvant_exact": "TiterMax Gold",
        "adjuvant_subgroup": "Other/None",
        "mouse_strain": "C57BL/6",
        "mouse_sex": "Female",
        "mouse_age": "7-8 wk",
        "challenge_strain": "Hypervirulent AB5075",
        "challenge_route": "IP",
        "challenge_inoculum": "2.07 x 10^6 CFU (~4 x LD50)",
        "followup_duration": "20 days",
        "ctrl_events": 0,
        "ctrl_n": 10,
        "count_reporting_basis": "Derived from stated 60% (AMEV1) and 80% (AMEV2) survival + n=10/group; all mock died within 2 days",
        "inferred_control_flag": False,
        "caveats_notes": "Independent experiment 1 in C57BL/6 mice with TiterMax Gold; AMEV1 contains pBamA and pOmpA (no LptD), whereas AMEV2 contains pOmp38 (no BamA, no LptD).",
        "arms": [
            {"arm_name": "rAMEV1", "platform": "Multivalent/chimeric", "ev": 6, "nv": 10, "basis": "60% of n=10"},
            {"arm_name": "rAMEV2 (exp. 1)", "platform": "Multivalent/chimeric", "ev": 8, "nv": 10, "basis": "80% of n=10"},
        ],
    },
    {
        "exp_id": "Prather_2026_exp2",
        "publication": "Prather 2026",
        "citation_short": "Prather et al., 2026, exp. 2",
        "doi": "10.3390/vaccines14030275",
        "cohort_stage": "Primary_14_pubs",
        "platform_summary": "Multivalent",
        "antigens_exact": "rAMEV2 (TrxA + pNlpE, pNucAB, pOmp38, pTonB, pZnuD) or rAMEV5 (TrxA + pNlpE, pNucAB, pTonB, pPlc1, pBauA, pBlp2, pOmp38; NO BamA, NO LptD)",
        "vacc_route": "SC",
        "vacc_dose_schedule": "10 ug SC on days 0, 14, 28 + AddaS03",
        "adjuvant_exact": "AddaS03",
        "adjuvant_subgroup": "Other/None",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Female",
        "mouse_age": "7-8 wk",
        "challenge_strain": "Hypervirulent AB5075",
        "challenge_route": "IP",
        "challenge_inoculum": "4.25 x 10^6 CFU",
        "followup_duration": "7 days",
        "ctrl_events": 0,
        "ctrl_n": 10,
        "count_reporting_basis": "Derived from stated 70% (AMEV2) and 90% (AMEV5) survival + n=10/group; all mock died within 2 days",
        "inferred_control_flag": False,
        "caveats_notes": "CRITICAL DRAFT ERROR: Draft Section 4.1.4 and Table 3 claim AMEV5 includes BamA and LptD; full text proves AMEV5 is TrxA + pNlpE, pNucAB, pTonB, pPlc1, pBauA, pBlp2, pOmp38 (NO BamA, NO LptD).",
        "arms": [
            {"arm_name": "rAMEV2 (exp. 2)", "platform": "Multivalent/chimeric", "ev": 7, "nv": 10, "basis": "70% of n=10"},
            {"arm_name": "rAMEV5", "platform": "Multivalent/chimeric", "ev": 9, "nv": 10, "basis": "90% of n=10"},
        ],
    },
    {
        "exp_id": "Singh_2016",
        "publication": "Singh 2016",
        "citation_short": "Singh et al., 2016",
        "doi": "10.3389/fmicb.2016.00158",
        "cohort_stage": "Supplementary_6_pubs",
        "platform_summary": "Other single",
        "antigens_exact": "FilF",
        "vacc_route": "SC",
        "vacc_dose_schedule": "20 ug SC on days 1, 14, 21 + CFA/IFA",
        "adjuvant_exact": "Freund's (CFA prime, IFA boosts)",
        "adjuvant_subgroup": "Freund's",
        "mouse_strain": "NR in supplementary audit",
        "mouse_sex": "NR in supplementary audit",
        "mouse_age": "NR in supplementary audit",
        "challenge_strain": "ATCC 19606",
        "challenge_route": "IT",
        "challenge_inoculum": "10^8 CFU (day 29)",
        "followup_duration": "7 days",
        "ctrl_events": 0,
        "ctrl_n": 10,
        "count_reporting_basis": "Derived from stated 50% survival through day 7 + n=10/group; all adjuvant controls died within 24-48 h",
        "inferred_control_flag": False,
        "caveats_notes": "Supplementary full-text addition; distinct from Singh et al., 2017 (rBamA, IN).",
        "arms": [
            {"arm_name": "FilF", "platform": "Other single", "ev": 5, "nv": 10, "basis": "50% of n=10"},
        ],
    },
    {
        "exp_id": "Garg_2016",
        "publication": "Garg 2016",
        "citation_short": "Garg et al., 2016",
        "doi": "10.1016/j.ijmm.2015.10.005",
        "cohort_stage": "Supplementary_6_pubs",
        "platform_summary": "Other single",
        "antigens_exact": "NucAb",
        "vacc_route": "IP",
        "vacc_dose_schedule": "25 ug IP on days 0, 14, 21 + CFA/IFA",
        "adjuvant_exact": "Freund's (CFA prime, IFA boosts)",
        "adjuvant_subgroup": "Freund's",
        "mouse_strain": "NR in supplementary audit",
        "mouse_sex": "NR in supplementary audit",
        "mouse_age": "NR in supplementary audit",
        "challenge_strain": "ATCC 19606",
        "challenge_route": "IT",
        "challenge_inoculum": "10^8 CFU (day 28)",
        "followup_duration": "7 days",
        "ctrl_events": 0,
        "ctrl_n": 10,
        "count_reporting_basis": "Derived from stated 20% survival at 48 h (surviving mice remained healthy through day 7; abstract reports 20% active survival) + n=10/group; 100% control mortality by 24 h",
        "inferred_control_flag": False,
        "caveats_notes": "Methods describe PBS+adjuvant control whereas Results/Fig. 4 call control 'unimmunized'; passive transfer (40% vs 0%, n=10/group) excluded. Tested by omission in sensitivity analysis.",
        "arms": [
            {"arm_name": "NucAb", "platform": "Other single", "ev": 2, "nv": 10, "basis": "20% at 48 h healthy through day 7, n=10"},
        ],
    },
    {
        "exp_id": "Skerniskyte_2019",
        "publication": "Skerniskyte 2019",
        "citation_short": "Skerniskyte et al., 2019",
        "doi": "10.1186/s12866-019-1615-3",
        "cohort_stage": "Supplementary_6_pubs",
        "platform_summary": "Other single",
        "antigens_exact": "Blp1 C-terminal fragment",
        "vacc_route": "IM",
        "vacc_dose_schedule": "2 ug IM on days 0, 14, 28 + CFA/IFA",
        "adjuvant_exact": "Freund's (CFA prime, IFA boosts)",
        "adjuvant_subgroup": "Freund's",
        "mouse_strain": "NR in supplementary audit",
        "mouse_sex": "NR in supplementary audit",
        "mouse_age": "NR in supplementary audit",
        "challenge_strain": "Clinical isolate AbIC I",
        "challenge_route": "IP",
        "challenge_inoculum": "1 x 10^8 CFU (day 42)",
        "followup_duration": "7 days",
        "ctrl_events": 0,
        "ctrl_n": 5,
        "count_reporting_basis": "Derived from stated 60% survival vs 'no survival in the control group' at day 7 + n=5/group",
        "inferred_control_flag": False,
        "caveats_notes": "Smallest group size (n=5/arm); separate passive serum transfer (n=3/group, 100% vs 0% at day 45 after 5 x 10^7 CFU) excluded.",
        "arms": [
            {"arm_name": "Blp1 C-terminal fragment", "platform": "Other single", "ev": 3, "nv": 5, "basis": "60% of n=5"},
        ],
    },
    {
        "exp_id": "Behrouz_2024",
        "publication": "Behrouz 2024",
        "citation_short": "Behrouz et al., 2024",
        "doi": "10.1038/s41598-024-79013-x",
        "cohort_stage": "Supplementary_6_pubs",
        "platform_summary": "Other single",
        "antigens_exact": "FliC-Omp22 fusion (Omp22 inserted into heterologous flagellin FliC) or Omp22 alone [active FliC-only carrier arm survival count graph-only, omitted]",
        "vacc_route": "SC",
        "vacc_dose_schedule": "10 ug SC on days 0, 14, 28 (no external adjuvant)",
        "adjuvant_exact": "None external (intrinsic FliC TLR5 agonist in fusion arm)",
        "adjuvant_subgroup": "Other/None",
        "mouse_strain": "NR in supplementary audit",
        "mouse_sex": "NR in supplementary audit",
        "mouse_age": "NR in supplementary audit",
        "challenge_strain": "MDR clinical isolate 58ST",
        "challenge_route": "IP",
        "challenge_inoculum": "1.9 x 10^7 CFU in 10% porcine mucin (day 56)",
        "followup_duration": "7 days",
        "ctrl_events": 0,
        "ctrl_n": 8,
        "count_reporting_basis": "Explicit x/n literally printed for active vaccine arms (6/8 and 3/8); PBS control 0/8 derived from 'all succumbed within two days' + n=8",
        "inferred_control_flag": False,
        "caveats_notes": "Active FliC-only arm (heterologous carrier control, n=8) had no printed x/n or % in text and was omitted; FliC-Omp22 presents a single A. baumannii antigen (Omp22) in flagellin carrier (classified under Other single, with sensitivity check under Multivalent/chimeric).",
        "arms": [
            {"arm_name": "FliC-Omp22 fusion", "platform": "Other single", "ev": 6, "nv": 8, "basis": "Printed 6/8 (75%)"},
            {"arm_name": "Omp22", "platform": "Other single", "ev": 3, "nv": 8, "basis": "Printed 3/8"},
        ],
    },
    {
        "exp_id": "Freitas_2024_exp_a",
        "publication": "Freitas 2024",
        "citation_short": "Freitas et al., 2024, exp. a",
        "doi": "10.1016/j.vaccine.2024.04.094",
        "cohort_stage": "Supplementary_6_pubs",
        "platform_summary": "Other single",
        "antigens_exact": "rCAM87009.1 (hypothetical adhesin)",
        "vacc_route": "IM",
        "vacc_dose_schedule": "25 ug rCAM87009.1 + 25 ug bio-AgNP IM on days 0, 14, 21",
        "adjuvant_exact": "Biogenic silver nanoparticles (bio-AgNP, 25 ug)",
        "adjuvant_subgroup": "Other/None",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Female",
        "mouse_age": "NR in supplementary audit",
        "challenge_strain": "MDR clinical isolate Ab02",
        "challenge_route": "IP",
        "challenge_inoculum": "DL100 (100 uL at ~1 x 10^5 CFU/mL; day 28)",
        "followup_duration": "7 days",
        "ctrl_events": 0,
        "ctrl_n": 7,
        "count_reporting_basis": "Derived from stated 100% protection + n=7/group; all 7 PBS+bio-AgNP controls died by day 3",
        "inferred_control_flag": False,
        "caveats_notes": "Independent formulation-matched control group (PBS + bio-AgNP, n=7) distinct from Alhydrogel control in exp. b (28 female BALB/c total across 4 groups of 7).",
        "arms": [
            {"arm_name": "rCAM87009.1 + bio-AgNP", "platform": "Other single", "ev": 7, "nv": 7, "basis": "100% of n=7"},
        ],
    },
    {
        "exp_id": "Freitas_2024_exp_b",
        "publication": "Freitas 2024",
        "citation_short": "Freitas et al., 2024, exp. b",
        "doi": "10.1016/j.vaccine.2024.04.094",
        "cohort_stage": "Supplementary_6_pubs",
        "platform_summary": "Other single",
        "antigens_exact": "rCAM87009.1 (hypothetical adhesin)",
        "vacc_route": "IM",
        "vacc_dose_schedule": "25 ug rCAM87009.1 + 15% Alhydrogel IM on days 0, 14, 21",
        "adjuvant_exact": "15% Alhydrogel (Alum)",
        "adjuvant_subgroup": "Alum",
        "mouse_strain": "BALB/c",
        "mouse_sex": "Female",
        "mouse_age": "NR in supplementary audit",
        "challenge_strain": "MDR clinical isolate Ab02",
        "challenge_route": "IP",
        "challenge_inoculum": "DL100 (100 uL at ~1 x 10^5 CFU/mL; day 28)",
        "followup_duration": "7 days",
        "ctrl_events": 0,
        "ctrl_n": 7,
        "count_reporting_basis": "Derived from stated 100% protection + n=7/group; PBS+Alhydrogel controls had 4 deaths on day 3 and remaining 3 on day 4",
        "inferred_control_flag": False,
        "caveats_notes": "Independent formulation-matched control group (PBS + Alhydrogel, n=7).",
        "arms": [
            {"arm_name": "rCAM87009.1 + Alhydrogel", "platform": "Other single", "ev": 7, "nv": 7, "basis": "100% of n=7"},
        ],
    },
    {
        "exp_id": "Jeffreys_2024",
        "publication": "Jeffreys 2024",
        "citation_short": "Jeffreys et al., 2024",
        "doi": "10.3390/vaccines12040358",
        "cohort_stage": "Supplementary_6_pubs",
        "platform_summary": "Multivalent",
        "antigens_exact": "rAMEV2 (TrxA leader + 5 peptides from NlpE, NucAB, Omp38, TonB, ZnuD; NO AMEV1)",
        "vacc_route": "SC",
        "vacc_dose_schedule": "10 ug SC on days 0, 14, 28 + AddaS03",
        "adjuvant_exact": "AddaS03",
        "adjuvant_subgroup": "Other/None",
        "mouse_strain": "C57BL/6 (PMC11054912 §3.4/3.6)",
        "mouse_sex": "NR in supplementary audit",
        "mouse_age": "7-8 wk (PMC11054912 §3.4)",
        "challenge_strain": "Clinical isolate Ci79",
        "challenge_route": "IN",
        "challenge_inoculum": "10^8 CFU (Fig. 7A caption: 1.01 x 10^8 CFU in 50 uL)",
        "followup_duration": "30 days",
        "ctrl_events": 0,
        "ctrl_n": 10,
        "count_reporting_basis": "Derived from stated 60% survival over 30 days + n=10/group; all mock succumbed by day 3",
        "inferred_control_flag": False,
        "caveats_notes": "Also in Track 2 (score 9/14; Draft Table 3 erroneously lists 'AMEV1 & AMEV2' whereas Jeffreys 2024 constructed AMEV2 only). Separate IP AB5075 pilot (80% survival, n=10 vaccinated) excluded because PBS mock group n was not stated.",
        "arms": [
            {"arm_name": "rAMEV2 (IN main)", "platform": "Multivalent/chimeric", "ev": 6, "nv": 10, "basis": "60% of n=10 over 30 d"},
        ],
    },
]


# =============================================================================
# 2. STATISTICAL FUNCTIONS
# =============================================================================

def rr_effect(ev: float, nv: float, ec: float, nc: float) -> dict[str, Any]:
    """Compute ln(RR) and sampling variance with +0.5 Haldane-Anscombe if any cell == 0."""
    bv = nv - ev
    bc = nc - ec
    zero_cell = (ev == 0) or (bv == 0) or (ec == 0) or (bc == 0)
    if zero_cell:
        a, b, c, d = ev + 0.5, bv + 0.5, ec + 0.5, bc + 0.5
    else:
        a, b, c, d = float(ev), float(bv), float(ec), float(bc)
    n1 = a + b
    n0 = c + d
    rr = (a / n1) / (c / n0)
    ln_rr = math.log(rr)
    var = (1.0 / a) - (1.0 / n1) + (1.0 / c) - (1.0 / n0)
    se = math.sqrt(var)
    z95 = norm_ppf(0.975)
    ci_low = math.exp(ln_rr - z95 * se)
    ci_high = math.exp(ln_rr + z95 * se)
    return {
        "zero_cell": zero_cell,
        "rr": rr,
        "ln_rr": ln_rr,
        "var": var,
        "se": se,
        "ci_low": ci_low,
        "ci_high": ci_high,
    }


def rd_effect(ev: float, nv: float, ec: float, nc: float, haldane_point_est: bool = False) -> dict[str, Any]:
    """
    Additive Risk Difference:
    - If haldane_point_est=False: raw proportions (ev/nv - ec/nc) for point estimate,
      Haldane-adjusted proportions (+0.5 / +1) for sampling variance when any 2x2 cell == 0.
    - If haldane_point_est=True (standard RevMan / meta::metabin incr=0.5 convention):
      +0.5 added to all 4 cells of zero-cell 2x2 tables for BOTH point estimate and variance.
    """
    bv = nv - ev
    bc = nc - ec
    zero_cell = (ev == 0) or (bv == 0) or (ec == 0) or (bc == 0)
    p1_raw = ev / nv
    p0_raw = ec / nc
    if zero_cell:
        n1_adj = nv + 1.0
        n0_adj = nc + 1.0
        p1 = (ev + 0.5) / n1_adj
        p0 = (ec + 0.5) / n0_adj
        var = (p1 * (1.0 - p1) / n1_adj) + (p0 * (1.0 - p0) / n0_adj)
        rd = (p1 - p0) if haldane_point_est else (p1_raw - p0_raw)
    else:
        var = (p1_raw * (1.0 - p1_raw) / nv) + (p0_raw * (1.0 - p0_raw) / nc)
        rd = p1_raw - p0_raw
    se = math.sqrt(var)
    z95 = norm_ppf(0.975)
    return {
        "zero_cell": zero_cell,
        "rd": rd,
        "var": var,
        "se": se,
        "ci_low": rd - z95 * se,
        "ci_high": rd + z95 * se,
    }


def pool_dl(y: np.ndarray, v: np.ndarray, exp_scale: bool = True) -> dict[str, Any]:
    """DerSimonian-Laird + HKSJ (with max(1, q)) + Prediction Interval + Fixed-Effect IV."""
    k = len(y)
    w = 1.0 / v
    sum_w = float(np.sum(w))
    mu_fe = float(np.sum(w * y) / sum_w)
    se_fe = math.sqrt(1.0 / sum_w)
    z95 = float(norm_ppf(0.975))

    if k == 1:
        est = math.exp(mu_fe) if exp_scale else mu_fe
        lo = math.exp(mu_fe - z95 * se_fe) if exp_scale else (mu_fe - z95 * se_fe)
        hi = math.exp(mu_fe + z95 * se_fe) if exp_scale else (mu_fe + z95 * se_fe)
        return {
            "k": 1,
            "est_dl": est,
            "wald_low": lo,
            "wald_high": hi,
            "hksj_low": None,
            "hksj_high": None,
            "Q": None,
            "df": 0,
            "p_Q": None,
            "I2": None,
            "tau2_dl": None,
            "weights_dl_pct": [100.0],
        }

    df = k - 1
    Q = float(np.sum(w * (y - mu_fe) ** 2))
    p_Q = float(1.0 - chi2_cdf(Q, df))
    C = float(sum_w - np.sum(w ** 2) / sum_w)
    tau2 = max(0.0, (Q - df) / C) if C > 0 else 0.0
    I2 = max(0.0, (Q - df) / Q * 100.0) if Q > 0 else 0.0

    w_star = 1.0 / (v + tau2)
    sum_w_star = float(np.sum(w_star))
    mu_dl = float(np.sum(w_star * y) / sum_w_star)
    se_dl = math.sqrt(1.0 / sum_w_star)
    z_dl = mu_dl / se_dl
    p_wald = float(2.0 * (1.0 - norm_cdf(abs(z_dl))))

    # HKSJ with max(1, q)
    q_raw = float(np.sum(w_star * (y - mu_dl) ** 2) / df)
    q_trunc = max(1.0, q_raw)
    se_hksj = se_dl * math.sqrt(q_trunc)
    se_hksj_raw = se_dl * math.sqrt(max(1e-15, q_raw))
    t_crit_hksj = float(t_ppf(0.975, df))
    t_stat = mu_dl / se_hksj
    p_hksj = float(2.0 * (1.0 - t_cdf(abs(t_stat), df)))

    # Prediction interval (standard df = k - 2 when k >= 3; also record df = k - 1 for audit check)
    if k >= 3:
        t_crit_pi_k2 = float(t_ppf(0.975, k - 2))
        se_pi = math.sqrt(tau2 + se_dl ** 2)
        pi_k2_low_raw = mu_dl - t_crit_pi_k2 * se_pi
        pi_k2_high_raw = mu_dl + t_crit_pi_k2 * se_pi

        se_pi_hksj = math.sqrt(tau2 + se_hksj ** 2)
        pi_k2_hksj_low_raw = mu_dl - t_crit_pi_k2 * se_pi_hksj
        pi_k2_hksj_high_raw = mu_dl + t_crit_pi_k2 * se_pi_hksj

        pi_k1_low_raw = mu_dl - t_crit_hksj * se_pi
        pi_k1_high_raw = mu_dl + t_crit_hksj * se_pi
    else:
        pi_k2_low_raw = pi_k2_high_raw = None
        pi_k2_hksj_low_raw = pi_k2_hksj_high_raw = None
        pi_k1_low_raw = pi_k1_high_raw = None

    tx = math.exp if exp_scale else (lambda x: x)
    weights_pct = [float(ws / sum_w_star * 100.0) for ws in w_star]

    return {
        "k": k,
        "mu_dl_raw": mu_dl,
        "se_dl_raw": se_dl,
        "est_dl": tx(mu_dl),
        "wald_low": tx(mu_dl - z95 * se_dl),
        "wald_high": tx(mu_dl + z95 * se_dl),
        "z_wald": z_dl,
        "p_wald": p_wald,
        "q_hksj_raw": q_raw,
        "q_hksj_used": q_trunc,
        "se_hksj_raw": se_hksj,
        "hksj_low": tx(mu_dl - t_crit_hksj * se_hksj) if k >= 3 else None,
        "hksj_high": tx(mu_dl + t_crit_hksj * se_hksj) if k >= 3 else None,
        "hksj_low_if_k2": tx(mu_dl - t_crit_hksj * se_hksj),
        "hksj_high_if_k2": tx(mu_dl + t_crit_hksj * se_hksj),
        "hksj_untruncated_low": tx(mu_dl - t_crit_hksj * se_hksj_raw),
        "hksj_untruncated_high": tx(mu_dl + t_crit_hksj * se_hksj_raw),
        "t_hksj": t_stat,
        "df_hksj": df,
        "p_hksj": p_hksj,
        "Q": Q,
        "df": df,
        "p_Q": p_Q,
        "I2": I2,
        "tau2_dl": tau2,
        "C": C,
        "pi_df_k2_low": tx(pi_k2_low_raw) if pi_k2_low_raw is not None else None,
        "pi_df_k2_high": tx(pi_k2_high_raw) if pi_k2_high_raw is not None else None,
        "pi_df_k2_hksjse_low": tx(pi_k2_hksj_low_raw) if pi_k2_hksj_low_raw is not None else None,
        "pi_df_k2_hksjse_high": tx(pi_k2_hksj_high_raw) if pi_k2_hksj_high_raw is not None else None,
        "pi_df_k1_low": tx(pi_k1_low_raw) if pi_k1_low_raw is not None else None,
        "pi_df_k1_high": tx(pi_k1_high_raw) if pi_k1_high_raw is not None else None,
        "est_fe_iv": tx(mu_fe),
        "fe_iv_low": tx(mu_fe - z95 * se_fe),
        "fe_iv_high": tx(mu_fe + z95 * se_fe),
        "z_fe_iv": mu_fe / se_fe,
        "p_fe_iv": float(2.0 * (1.0 - norm_cdf(abs(mu_fe / se_fe)))),
        "weights_dl_pct": weights_pct,
    }


def mantel_haenszel_rr(rows: list[dict[str, Any]], continuity_correction: bool = False) -> dict[str, Any]:
    """
    Mantel-Haenszel pooled Risk Ratio with Greenland-Robins (1985) variance.
    Note: Exact uncorrected MH RR requires NO continuity correction as long as sum(c_i * n1_i / N_i) > 0!
    """
    num = 0.0
    den = 0.0
    var_num = 0.0
    for r in rows:
        ev, nv, ec, nc = float(r["ev"]), float(r["nv"]), float(r["ec"]), float(r["nc"])
        bv, bc = nv - ev, nc - ec
        zero_cell = (ev == 0) or (bv == 0) or (ec == 0) or (bc == 0)
        if continuity_correction and zero_cell:
            a, b, c, d = ev + 0.5, bv + 0.5, ec + 0.5, bc + 0.5
        else:
            a, b, c, d = ev, bv, ec, bc
        n1 = a + b
        n0 = c + d
        N = n1 + n0
        m1 = a + c
        num += a * n0 / N
        den += c * n1 / N
        var_num += (m1 * n1 * n0 - a * c * N) / (N ** 2)
    rr_mh = num / den
    var_ln_rr_mh = var_num / (num * den)
    se_ln_rr_mh = math.sqrt(var_ln_rr_mh)
    z95 = float(norm_ppf(0.975))
    ln_rr = math.log(rr_mh)
    z_val = ln_rr / se_ln_rr_mh
    return {
        "continuity_correction": continuity_correction,
        "rr_mh": rr_mh,
        "ci_low": math.exp(ln_rr - z95 * se_ln_rr_mh),
        "ci_high": math.exp(ln_rr + z95 * se_ln_rr_mh),
        "se_ln_rr": se_ln_rr_mh,
        "z": z_val,
        "p": float(2.0 * (1.0 - norm_cdf(abs(z_val)))),
    }


def egger_test(y: np.ndarray, se: np.ndarray) -> dict[str, Any]:
    """Egger's linear regression of SND = y/se on precision = 1/se."""
    k = len(y)
    precision = 1.0 / se
    snd = y / se
    res = linregress_exact(precision, snd)
    # In linregress(x=precision, y=snd): snd = intercept + slope * precision
    intercept = float(res["intercept"])
    se_intercept = float(res["intercept_stderr"])
    t_stat = intercept / se_intercept
    df = k - 2
    p_val = float(2.0 * (1.0 - t_cdf(abs(t_stat), df)))
    return {
        "k": k,
        "df": df,
        "intercept": intercept,
        "se_intercept": se_intercept,
        "t": t_stat,
        "p": p_val,
        "slope_fe_limit": float(res["slope"]),
        "se_slope": float(res["stderr"]),
    }


def q_between_subgroups(subgroup_results: list[dict[str, Any]]) -> dict[str, Any]:
    """Mixed-effects Q_between test across pooled subgroup DL estimates (with k >= 2 or included groups)."""
    mus = np.array([g["mu_dl_raw"] for g in subgroup_results], dtype=float)
    ses = np.array([g["se_dl_raw"] for g in subgroup_results], dtype=float)
    ws = 1.0 / (ses ** 2)
    mu_bar = float(np.sum(ws * mus) / np.sum(ws))
    Q_bet = float(np.sum(ws * (mus - mu_bar) ** 2))
    df = len(subgroup_results) - 1
    p_bet = float(1.0 - chi2_cdf(Q_bet, df)) if df > 0 else None
    return {"Q_between": Q_bet, "df": df, "p_between": p_bet}


# =============================================================================
# 3. BUILD EXPERIMENT-LEVEL & ARM-LEVEL TABLES AND RUN PRIMARY META-ANALYSIS
# =============================================================================

exp_rows: list[dict[str, Any]] = []
arm_rows: list[dict[str, Any]] = []

arm_idx = 0
for idx, exp in enumerate(EXPERIMENTS_META, start=1):
    ev_sum = sum(a["ev"] for a in exp["arms"])
    nv_sum = sum(a["nv"] for a in exp["arms"])
    ec = exp["ctrl_events"]
    nc = exp["ctrl_n"]
    rr_res = rr_effect(ev_sum, nv_sum, ec, nc)
    rd_raw_res = rd_effect(ev_sum, nv_sum, ec, nc, haldane_point_est=False)
    rd_cc_res = rd_effect(ev_sum, nv_sum, ec, nc, haldane_point_est=True)
    platforms_present = sorted(set(a["platform"] for a in exp["arms"]))

    exp_row = {
        "exp_num": idx,
        "exp_id": exp["exp_id"],
        "publication": exp["publication"],
        "citation_short": exp["citation_short"],
        "doi": exp["doi"],
        "cohort_stage": exp["cohort_stage"],
        "platform_summary": exp["platform_summary"],
        "antigens_exact": exp["antigens_exact"],
        "vacc_route": exp["vacc_route"],
        "vacc_dose_schedule": exp["vacc_dose_schedule"],
        "adjuvant_exact": exp["adjuvant_exact"],
        "adjuvant_subgroup": exp["adjuvant_subgroup"],
        "mouse_strain": exp["mouse_strain"],
        "mouse_sex": exp["mouse_sex"],
        "mouse_age": exp["mouse_age"],
        "challenge_strain": exp["challenge_strain"],
        "challenge_route": exp["challenge_route"],
        "challenge_inoculum": exp["challenge_inoculum"],
        "followup_duration": exp["followup_duration"],
        "arms_count": len(exp["arms"]),
        "ev": ev_sum,
        "nv": nv_sum,
        "vacc_pct": round(ev_sum / nv_sum * 100.0, 2),
        "ec": ec,
        "nc": nc,
        "ctrl_pct": round(ec / nc * 100.0, 2),
        "count_reporting_basis": exp["count_reporting_basis"],
        "has_zero_cell": rr_res["zero_cell"],
        "ctrl_zero_events": (ec == 0),
        "inferred_control_flag": exp["inferred_control_flag"],
        "RR": rr_res["rr"],
        "RR_ci_low": rr_res["ci_low"],
        "RR_ci_high": rr_res["ci_high"],
        "lnRR": rr_res["ln_rr"],
        "var_lnRR": rr_res["var"],
        "se_lnRR": rr_res["se"],
        "RD_raw_pp": rd_raw_res["rd"] * 100.0,
        "RD_raw_ci_low_pp": rd_raw_res["ci_low"] * 100.0,
        "RD_raw_ci_high_pp": rd_raw_res["ci_high"] * 100.0,
        "RD_haldane_pp": rd_cc_res["rd"] * 100.0,
        "RD_haldane_ci_low_pp": rd_cc_res["ci_low"] * 100.0,
        "RD_haldane_ci_high_pp": rd_cc_res["ci_high"] * 100.0,
        "var_RD": rd_raw_res["var"],
        "se_RD": rd_raw_res["se"],
        "caveats_notes": exp["caveats_notes"],
    }
    exp_rows.append(exp_row)

    for a in exp["arms"]:
        arm_idx += 1
        a_rr = rr_effect(a["ev"], a["nv"], ec, nc)
        arm_rows.append({
            "arm_index": arm_idx,
            "exp_id": exp["exp_id"],
            "publication": exp["publication"],
            "doi": exp["doi"],
            "arm_name": a["arm_name"],
            "platform_subgroup": a["platform"],
            "challenge_route": exp["challenge_route"],
            "adjuvant_subgroup": exp["adjuvant_subgroup"],
            "vacc_events": a["ev"],
            "vacc_n": a["nv"],
            "vacc_pct": round(a["ev"] / a["nv"] * 100.0, 2),
            "ctrl_events_shared": ec,
            "ctrl_n_shared": nc,
            "arms_in_experiment": len(exp["arms"]),
            "platforms_in_experiment": len(platforms_present),
            "count_reporting_basis": a["basis"],
            "arm_RR_vs_full_ctrl": round(a_rr["rr"], 4),
            "arm_RR_ci_low": round(a_rr["ci_low"], 4),
            "arm_RR_ci_high": round(a_rr["ci_high"], 4),
        })

# Primary pooling (k = 22)
y_rr_all = np.array([r["lnRR"] for r in exp_rows])
v_rr_all = np.array([r["var_lnRR"] for r in exp_rows])
se_rr_all = np.array([r["se_lnRR"] for r in exp_rows])
overall_rr = pool_dl(y_rr_all, v_rr_all, exp_scale=True)

y_rd_raw_all = np.array([r["RD_raw_pp"] / 100.0 for r in exp_rows])
y_rd_cc_all = np.array([r["RD_haldane_pp"] / 100.0 for r in exp_rows])
v_rd_all = np.array([r["var_RD"] for r in exp_rows])
overall_rd_raw = pool_dl(y_rd_raw_all, v_rd_all, exp_scale=False)
overall_rd_cc = pool_dl(y_rd_cc_all, v_rd_all, exp_scale=False)

for i, r in enumerate(exp_rows):
    r["weight_DL_RR_pct"] = overall_rr["weights_dl_pct"][i]
    r["weight_DL_RD_raw_pct"] = overall_rd_raw["weights_dl_pct"][i]
    r["weight_DL_RD_haldane_pct"] = overall_rd_cc["weights_dl_pct"][i]

mh_uncorrected = mantel_haenszel_rr(exp_rows, continuity_correction=False)
mh_corrected = mantel_haenszel_rr(exp_rows, continuity_correction=True)
egger_all = egger_test(y_rr_all, se_rr_all)

var_ec0 = [r["var_lnRR"] for r in exp_rows if r["ctrl_zero_events"]]
var_ec_pos = [r["var_lnRR"] for r in exp_rows if not r["ctrl_zero_events"]]


# =============================================================================
# 4. SUBGROUP ANALYSES (PLATFORM 30 UNITS, CHALLENGE ROUTE, ADJUVANT)
# =============================================================================

PLATFORM_ORDER = [
    "OmpA",
    "BamA / Oma87",
    "OmpW",
    "Other single",
    "Multivalent/chimeric",
]


def build_platform_units(split_mode: str = "equal_by_platform", reclassify_behrouz_fusion: bool = False) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    units: list[dict[str, Any]] = []
    unit_id = 0
    for exp in EXPERIMENTS_META:
        arms_local = []
        for a in exp["arms"]:
            plat = a["platform"]
            if reclassify_behrouz_fusion and exp["exp_id"] == "Behrouz_2024" and a["arm_name"] == "FliC-Omp22 fusion":
                plat = "Multivalent/chimeric"
            arms_local.append({**a, "platform_eff": plat})

        by_plat: dict[str, list[dict[str, Any]]] = {}
        for a in arms_local:
            by_plat.setdefault(a["platform_eff"], []).append(a)

        m_plats = len(by_plat)
        total_arms_exp = len(arms_local)
        for plat in PLATFORM_ORDER:
            if plat not in by_plat:
                continue
            unit_id += 1
            p_arms = by_plat[plat]
            ev_u = sum(x["ev"] for x in p_arms)
            nv_u = sum(x["nv"] for x in p_arms)
            if split_mode == "equal_by_platform":
                ec_u = exp["ctrl_events"] / m_plats
                nc_u = exp["ctrl_n"] / m_plats
            elif split_mode == "proportional_to_arms":
                frac = len(p_arms) / total_arms_exp
                ec_u = exp["ctrl_events"] * frac
                nc_u = exp["ctrl_n"] * frac
            else:
                raise ValueError(split_mode)
            eff = rr_effect(ev_u, nv_u, ec_u, nc_u)
            units.append({
                "unit_id": unit_id,
                "exp_id": exp["exp_id"],
                "publication": exp["publication"],
                "platform_subgroup": plat,
                "arms_pooled_in_unit": "; ".join(x["arm_name"] for x in p_arms),
                "arms_count_in_unit": len(p_arms),
                "platforms_in_exp": m_plats,
                "vacc_events": ev_u,
                "vacc_n": nv_u,
                "vacc_pct": round(ev_u / nv_u * 100.0, 2),
                "ctrl_events_orig": exp["ctrl_events"],
                "ctrl_n_orig": exp["ctrl_n"],
                "ctrl_events_split": round(ec_u, 4),
                "ctrl_n_split": round(nc_u, 4),
                "has_zero_cell": eff["zero_cell"],
                "RR": eff["rr"],
                "RR_ci_low": eff["ci_low"],
                "RR_ci_high": eff["ci_high"],
                "lnRR": eff["ln_rr"],
                "var_lnRR": eff["var"],
                "se_lnRR": eff["se"],
            })

    summary_by_plat: dict[str, Any] = {}
    pooled_for_qbet = []
    pooled_for_qbet_incl_ompa = []
    for plat in PLATFORM_ORDER:
        sub = [u for u in units if u["platform_subgroup"] == plat]
        y_sub = np.array([u["lnRR"] for u in sub])
        v_sub = np.array([u["var_lnRR"] for u in sub])
        res = pool_dl(y_sub, v_sub, exp_scale=True)
        for u, wt in zip(sub, res["weights_dl_pct"]):
            u["subgroup_weight_pct"] = wt
        ev_tot = sum(u["vacc_events"] for u in sub)
        nv_tot = sum(u["vacc_n"] for u in sub)
        pubs = sorted(set(u["publication"] for u in sub))
        arms_tot = sum(u["arms_count_in_unit"] for u in sub)
        summary_by_plat[plat] = {
            "k_units": len(sub),
            "publications_count": len(pubs),
            "arms_count": arms_tot,
            "vacc_events": ev_tot,
            "vacc_n": nv_tot,
            "vacc_pct": ev_tot / nv_tot * 100.0,
            "pooled": res,
        }
        if len(sub) >= 2:
            pooled_for_qbet.append(res)
        pooled_for_qbet_incl_ompa.append(
            res if len(sub) >= 2 else {"mu_dl_raw": sub[0]["lnRR"], "se_dl_raw": sub[0]["se_lnRR"]}
        )

    q_bet_4 = q_between_subgroups(pooled_for_qbet)
    q_bet_5 = q_between_subgroups(pooled_for_qbet_incl_ompa)
    return units, {
        "total_units": len(units),
        "total_arms": sum(u["arms_count_in_unit"] for u in units),
        "subgroups": summary_by_plat,
        "Q_between_excluding_k1_OmpA": q_bet_4,
        "Q_between_all_5_groups": q_bet_5,
    }


platform_units_30, platform_summary_equal = build_platform_units("equal_by_platform", False)
_, platform_summary_armprop = build_platform_units("proportional_to_arms", False)
_, platform_summary_behrouz_split = build_platform_units("equal_by_platform", True)

# Challenge route subgroups
route_summary: dict[str, Any] = {}
route_list_for_q = []
for route in ["IP", "IN", "IT"]:
    sub = [r for r in exp_rows if r["challenge_route"] == route]
    res_rr = pool_dl(np.array([r["lnRR"] for r in sub]), np.array([r["var_lnRR"] for r in sub]), True)
    res_rd_cc = pool_dl(np.array([r["RD_haldane_pp"] / 100.0 for r in sub]), np.array([r["var_RD"] for r in sub]), False)
    res_rd_raw = pool_dl(np.array([r["RD_raw_pp"] / 100.0 for r in sub]), np.array([r["var_RD"] for r in sub]), False)
    ev_s, nv_s = sum(r["ev"] for r in sub), sum(r["nv"] for r in sub)
    ec_s, nc_s = sum(r["ec"] for r in sub), sum(r["nc"] for r in sub)
    pubs = sorted(set(r["publication"] for r in sub))
    route_summary[route] = {
        "k": len(sub),
        "publications_count": len(pubs),
        "vacc_events": ev_s,
        "vacc_n": nv_s,
        "vacc_pct": ev_s / nv_s * 100.0,
        "ctrl_events": ec_s,
        "ctrl_n": nc_s,
        "ctrl_pct": ec_s / nc_s * 100.0,
        "RR": res_rr,
        "RD_haldane_both": res_rd_cc,
        "RD_raw_pt_haldane_var": res_rd_raw,
    }
    route_list_for_q.append(res_rr)
route_q_between = q_between_subgroups(route_list_for_q)

# Adjuvant subgroups
adj_summary: dict[str, Any] = {}
adj_list_for_q = []
for adj in ["Alum", "Freund's", "Other/None"]:
    sub = [r for r in exp_rows if r["adjuvant_subgroup"] == adj]
    res_rr = pool_dl(np.array([r["lnRR"] for r in sub]), np.array([r["var_lnRR"] for r in sub]), True)
    res_rd_cc = pool_dl(np.array([r["RD_haldane_pp"] / 100.0 for r in sub]), np.array([r["var_RD"] for r in sub]), False)
    res_rd_raw = pool_dl(np.array([r["RD_raw_pp"] / 100.0 for r in sub]), np.array([r["var_RD"] for r in sub]), False)
    ev_s, nv_s = sum(r["ev"] for r in sub), sum(r["nv"] for r in sub)
    ec_s, nc_s = sum(r["ec"] for r in sub), sum(r["nc"] for r in sub)
    pubs = sorted(set(r["publication"] for r in sub))
    adj_summary[adj] = {
        "k": len(sub),
        "publications_count": len(pubs),
        "vacc_events": ev_s,
        "vacc_n": nv_s,
        "vacc_pct": ev_s / nv_s * 100.0,
        "ctrl_events": ec_s,
        "ctrl_n": nc_s,
        "ctrl_pct": ec_s / nc_s * 100.0,
        "RR": res_rr,
        "RD_haldane_both": res_rd_cc,
        "RD_raw_pt_haldane_var": res_rd_raw,
    }
    adj_list_for_q.append(res_rr)
adj_q_between = q_between_subgroups(adj_list_for_q)


# =============================================================================
# 5. SENSITIVITY ANALYSES & LEAVE-ONE-OUT (20 PUBLICATIONS & 22 EXPERIMENTS)
# =============================================================================

# Leave-one-publication-out (20 unique publications in order of appearance)
pub_order: list[str] = []
for r in exp_rows:
    if r["publication"] not in pub_order:
        pub_order.append(r["publication"])

loo_pub_rows: list[dict[str, Any]] = []
for pub in pub_order:
    rem = [r for r in exp_rows if r["publication"] != pub]
    omitted = [r for r in exp_rows if r["publication"] == pub]
    res_rr = pool_dl(np.array([r["lnRR"] for r in rem]), np.array([r["var_lnRR"] for r in rem]), True)
    res_rd_cc = pool_dl(np.array([r["RD_haldane_pp"] / 100.0 for r in rem]), np.array([r["var_RD"] for r in rem]), False)
    res_rd_raw = pool_dl(np.array([r["RD_raw_pp"] / 100.0 for r in rem]), np.array([r["var_RD"] for r in rem]), False)
    ev_r = sum(r["ev"] for r in rem)
    nv_r = sum(r["nv"] for r in rem)
    ec_r = sum(r["ec"] for r in rem)
    nc_r = sum(r["nc"] for r in rem)
    loo_pub_rows.append({
        "omitted_publication": pub,
        "doi": omitted[0]["doi"],
        "experiments_omitted": len(omitted),
        "k_remaining": len(rem),
        "mice_remaining": nv_r + nc_r,
        "vacc_surv_remaining": f"{ev_r}/{nv_r}",
        "ctrl_surv_remaining": f"{ec_r}/{nc_r}",
        "RR_DL": res_rr["est_dl"],
        "RR_wald_low": res_rr["wald_low"],
        "RR_wald_high": res_rr["wald_high"],
        "RR_hksj_low": res_rr["hksj_low"],
        "RR_hksj_high": res_rr["hksj_high"],
        "q_hksj": res_rr["q_hksj_raw"],
        "Q": res_rr["Q"],
        "p_Q": res_rr["p_Q"],
        "I2_pct": res_rr["I2"],
        "tau2_DL": res_rr["tau2_dl"],
        "RD_haldane_DL_pp": res_rd_cc["est_dl"] * 100.0,
        "RD_haldane_hksj_low_pp": res_rd_cc["hksj_low"] * 100.0,
        "RD_haldane_hksj_high_pp": res_rd_cc["hksj_high"] * 100.0,
        "I2_RD_haldane_pct": res_rd_cc["I2"],
        "RD_raw_DL_pp": res_rd_raw["est_dl"] * 100.0,
        "RD_raw_hksj_low_pp": res_rd_raw["hksj_low"] * 100.0,
        "RD_raw_hksj_high_pp": res_rd_raw["hksj_high"] * 100.0,
        "I2_RD_raw_pct": res_rd_raw["I2"],
    })

# Leave-one-experiment-out (22 experiments)
loo_exp_rows: list[dict[str, Any]] = []
for r_omit in exp_rows:
    rem = [r for r in exp_rows if r["exp_id"] != r_omit["exp_id"]]
    res_rr = pool_dl(np.array([r["lnRR"] for r in rem]), np.array([r["var_lnRR"] for r in rem]), True)
    loo_exp_rows.append({
        "omitted_exp_id": r_omit["exp_id"],
        "RR_DL": res_rr["est_dl"],
        "RR_hksj_low": res_rr["hksj_low"],
        "RR_hksj_high": res_rr["hksj_high"],
    })


def run_subset_sensitivity(label: str, subset: list[dict[str, Any]]) -> dict[str, Any]:
    res_rr = pool_dl(np.array([r["lnRR"] for r in subset]), np.array([r["var_lnRR"] for r in subset]), True)
    res_rd_cc = pool_dl(np.array([r["RD_haldane_pp"] / 100.0 for r in subset]), np.array([r["var_RD"] for r in subset]), False)
    res_rd_raw = pool_dl(np.array([r["RD_raw_pp"] / 100.0 for r in subset]), np.array([r["var_RD"] for r in subset]), False)
    ev_s, nv_s = sum(r["ev"] for r in subset), sum(r["nv"] for r in subset)
    ec_s, nc_s = sum(r["ec"] for r in subset), sum(r["nc"] for r in subset)
    mh_raw = mantel_haenszel_rr(subset, False) if sum(r["ec"] for r in subset) > 0 else None
    mh_cc = mantel_haenszel_rr(subset, True)
    eg = egger_test(np.array([r["lnRR"] for r in subset]), np.array([r["se_lnRR"] for r in subset])) if len(subset) >= 3 else None
    return {
        "label": label,
        "k": len(subset),
        "publications_count": len(set(r["publication"] for r in subset)),
        "total_mice": nv_s + nc_s,
        "vacc_events": ev_s,
        "vacc_n": nv_s,
        "vacc_pct": ev_s / nv_s * 100.0,
        "ctrl_events": ec_s,
        "ctrl_n": nc_s,
        "ctrl_pct": ec_s / nc_s * 100.0,
        "RR": res_rr,
        "RD_haldane_both": res_rd_cc,
        "RD_raw_pt_haldane_var": res_rd_raw,
        "MH_uncorrected": mh_raw,
        "MH_continuity_corrected": mh_cc,
        "Egger": eg,
    }


sens_omit_inferred_ctrl = run_subset_sensitivity(
    "Omit inferred-control experiments (Rasooli 2020, Abdollahi 2023)",
    [r for r in exp_rows if not r["inferred_control_flag"]],
)
sens_ip_only = run_subset_sensitivity(
    "Intraperitoneal (IP) challenge only",
    [r for r in exp_rows if r["challenge_route"] == "IP"],
)
sens_omit_garg = run_subset_sensitivity(
    "Omit Garg 2016 (48-h survival stated healthy to day 7)",
    [r for r in exp_rows if r["exp_id"] != "Garg_2016"],
)
sens_pre_supp_15 = run_subset_sensitivity(
    "Pre-supplementary 14-publication / 15-experiment dataset",
    [r for r in exp_rows if r["cohort_stage"] == "Primary_14_pubs"],
)
sens_ec_pos_4 = run_subset_sensitivity(
    "Experiments with >=1 surviving control (ec > 0: Araujo/Conde 2021, Mansouri 2023, Huang 2015, Hagag 2022)",
    [r for r in exp_rows if not r["ctrl_zero_events"]],
)
sens_ec_zero_18 = run_subset_sensitivity(
    "Experiments with zero surviving controls (ec = 0, k = 18)",
    [r for r in exp_rows if r["ctrl_zero_events"]],
)
sens_7d_plus_known = run_subset_sensitivity(
    "Experiments with explicitly stated follow-up >= 7 days (omitting 96-h Rasooli 2020 & Mansouri 2023 and NS Abdollahi 2023, k = 19)",
    [r for r in exp_rows if r["exp_id"] not in {"Rasooli_2020", "Mansouri_2023", "Abdollahi_2023"}],
)


# =============================================================================
# 6. WITHIN-STUDY HEAD-TO-HEAD MULTIVALENT / FUSION VS SINGLE-ANTIGEN COMPARISONS
# =============================================================================

HEAD_TO_HEAD_STUDIES = [
    ("Badmasti_2015", "Multivalent vs single subunits"),
    ("Akbari_2022", "Dual-loop LCL chimera vs single-loop LCL constructs"),
    ("Qamsari_2020", "4-loop ZnuD LCL hybrid vs whole rZnuD & single-loop LCL constructs"),
    ("Hagag_2022", "Wza+YiaD bivalent protein pool vs Wza or YiaD alone"),
    ("Mansouri_2023", "Oma87+Bap bivalent vs Oma87 or Bap alone"),
    ("Chaudhuri_2022", "Multi-loop BauA LCL hybrids vs whole rBauA & single-loop LCL constructs"),
    ("Behrouz_2024", "FliC-Omp22 flagellin-fusion chimera vs Omp22 alone"),
]

within_study_comparisons: list[dict[str, Any]] = []
for exp_id, comp_type in HEAD_TO_HEAD_STUDIES:
    exp = next(e for e in EXPERIMENTS_META if e["exp_id"] == exp_id)
    if exp_id == "Behrouz_2024":
        multi_arms = [a for a in exp["arms"] if "FliC" in a["arm_name"]]
        single_arms = [a for a in exp["arms"] if "FliC" not in a["arm_name"]]
    else:
        multi_arms = [a for a in exp["arms"] if a["platform"] == "Multivalent/chimeric"]
        single_arms = [a for a in exp["arms"] if a["platform"] != "Multivalent/chimeric"]

    ev_m, nv_m = sum(a["ev"] for a in multi_arms), sum(a["nv"] for a in multi_arms)
    ev_s, nv_s = sum(a["ev"] for a in single_arms), sum(a["nv"] for a in single_arms)
    best_single = max(single_arms, key=lambda a: a["ev"] / a["nv"])
    worst_single = min(single_arms, key=lambda a: a["ev"] / a["nv"])

    # Exact 2-sided Fisher's exact test: pooled multivalent/fusion vs pooled single arms
    _, p_fisher_pooled = fisher_exact_two_sided([[ev_m, nv_m - ev_m], [ev_s, nv_s - ev_s]])
    # Exact Fisher test: best multivalent arm vs best single arm
    best_multi = max(multi_arms, key=lambda a: a["ev"] / a["nv"])
    _, p_fisher_best = fisher_exact_two_sided([
        [best_multi["ev"], best_multi["nv"] - best_multi["ev"]],
        [best_single["ev"], best_single["nv"] - best_single["ev"]],
    ])
    rr_m_vs_s = rr_effect(ev_m, nv_m, ev_s, nv_s)
    rd_m_vs_s = rd_effect(ev_m, nv_m, ev_s, nv_s)

    if (best_multi["ev"] / best_multi["nv"]) > (best_single["ev"] / best_single["nv"]):
        direction = "Multivalent/fusion numerically higher than every single-antigen arm"
    elif math.isclose(best_multi["ev"] / best_multi["nv"], best_single["ev"] / best_single["nv"]):
        direction = "Tie / ceiling (100% in all arms)"
    else:
        direction = "Multivalent/multi-loop did NOT exceed best single-antigen arm"

    within_study_comparisons.append({
        "exp_id": exp_id,
        "publication": exp["publication"],
        "comparison_type": comp_type,
        "multi_arms_detail": [f"{a['arm_name']}: {a['ev']}/{a['nv']} ({a['ev']/a['nv']*100:.1f}%)" for a in multi_arms],
        "single_arms_detail": [f"{a['arm_name']}: {a['ev']}/{a['nv']} ({a['ev']/a['nv']*100:.1f}%)" for a in single_arms],
        "multi_pooled": f"{ev_m}/{nv_m} ({ev_m/nv_m*100:.1f}%)",
        "single_pooled": f"{ev_s}/{nv_s} ({ev_s/nv_s*100:.1f}%)",
        "best_multi": f"{best_multi['arm_name']} {best_multi['ev']}/{best_multi['nv']} ({best_multi['ev']/best_multi['nv']*100:.1f}%)",
        "best_single": f"{best_single['arm_name']} {best_single['ev']}/{best_single['nv']} ({best_single['ev']/best_single['nv']*100:.1f}%)",
        "worst_single": f"{worst_single['arm_name']} {worst_single['ev']}/{worst_single['nv']} ({worst_single['ev']/worst_single['nv']*100:.1f}%)",
        "RR_multi_vs_single_pooled": rr_m_vs_s["rr"],
        "RR_ci_low": rr_m_vs_s["ci_low"],
        "RR_ci_high": rr_m_vs_s["ci_high"],
        "RD_multi_vs_single_pp": rd_m_vs_s["rd"] * 100.0,
        "fisher_p_pooled": float(p_fisher_pooled),
        "fisher_p_best_vs_best": float(p_fisher_best),
        "direction": direction,
    })


# =============================================================================
# 7. SYRCLE 10-DOMAIN + REPORTING QUALITY MATRIX (ALL 20 PUBLICATIONS)
# =============================================================================
# Coding strictly from INLINE_USER_AUDIT_TABLES.md & A._baumannii_mouse_vaccine_challenge_report.md:
# 'Unclear (NS)' = not stated in full text; 'Reported (method NS)', 'Reported (histopathology only)', etc.

SYRCLE_20: list[dict[str, Any]] = [
    {
        "publication": "Badmasti 2015",
        "doi": "10.1016/j.molimm.2015.07.031",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age stated; weight/baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS; survival is objective)",
        "D8_incomplete_outcome_data": "Unclear (exclusions NS)",
        "D9_selective_reporting": "Low/Unclear (all ATCC 19606 arms reported; AB-44 tested subset arms)",
        "D10_other_bias": "Unclear (sample-size calculation NS)",
        "strain_sex_reported": "Yes (C57BL/6, female)",
        "age_reported": "Yes (6-8 wk)",
        "followup_reported": "Yes (7 d)",
        "humane_endpoint_reported": "No (NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Unclear/NS in excerpt",
        "ethics_approval_number": "No (NS)",
    },
    {
        "publication": "Singh 2017",
        "doi": "10.1038/s41598-017-12789-3",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age stated; weight/baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Unclear (exclusions NS)",
        "D9_selective_reporting": "Low/Unclear (no protocol; survival & CFU cohorts reported)",
        "D10_other_bias": "Unclear (sample-size calculation NS)",
        "strain_sex_reported": "Yes (BALB/c, female)",
        "age_reported": "Yes (6-8 wk)",
        "followup_reported": "Yes (7 d)",
        "humane_endpoint_reported": "No (NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes (Animal Ethics Committee)",
        "ethics_approval_number": "No (NS)",
    },
    {
        "publication": "Rasooli 2020",
        "doi": "10.1007/s10989-020-10056-0",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age/weight range stated; group baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Unclear (overall n=8/n=6 stated; subgroup n by challenge dose NS)",
        "D9_selective_reporting": "Unclear (control survival inferred from prose, not tabulated)",
        "D10_other_bias": "Unclear (sample-size calculation NS)",
        "strain_sex_reported": "Yes (BALB/c, female)",
        "age_reported": "Yes (6-8 wk; 22-25 g)",
        "followup_reported": "Yes (4 d / 96 h)",
        "humane_endpoint_reported": "No (NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes (Shahed University)",
        "ethics_approval_number": "No (NS)",
    },
    {
        "publication": "Araujo/Conde 2021",
        "doi": "10.1016/j.micinf.2021.104801",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age stated; weight/baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Unclear (exclusions NS)",
        "D9_selective_reporting": "Low/Unclear",
        "D10_other_bias": "Unclear (sample-size calculation NS)",
        "strain_sex_reported": "Yes (C57BL/6, female)",
        "age_reported": "Yes (8 wk)",
        "followup_reported": "Yes (7 d)",
        "humane_endpoint_reported": "Yes (severe symptoms/suffering/inability to access food or water)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes",
        "ethics_approval_number": "Yes (LW-21/18)",
    },
    {
        "publication": "Mansouri 2023",
        "doi": "10.1016/j.intimp.2023.110650",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/weight stated; age & baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Unclear (exclusions NS)",
        "D9_selective_reporting": "Low/Unclear",
        "D10_other_bias": "Unclear (sample-size calculation NS; shared control across 3 arms)",
        "strain_sex_reported": "Yes (BALB/c, female)",
        "age_reported": "No (NS; weight 20-25 g stated)",
        "followup_reported": "Yes (4 d / 96 h)",
        "humane_endpoint_reported": "No (NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes (Shahed University)",
        "ethics_approval_number": "Yes (IR.SHAHED.REC.1400.071)",
    },
    {
        "publication": "Huang 2015",
        "doi": "10.1016/j.vaccine.2015.07.031",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age stated; weight/baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (SPF housing stated; random housing NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Unclear (exclusions NS)",
        "D9_selective_reporting": "Low/Unclear",
        "D10_other_bias": "Unclear (sample-size calculation NS)",
        "strain_sex_reported": "Yes (ICR, female)",
        "age_reported": "Yes (6-8 wk)",
        "followup_reported": "Yes (7 d)",
        "humane_endpoint_reported": "No (NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes (committee named)",
        "ethics_approval_number": "No (NS)",
    },
    {
        "publication": "Abdollahi 2023",
        "doi": "10.1016/j.micpath.2023.106317",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age stated; baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Concern/Discrepancy (Methods §2.16 states n=5 vs §2.17.2 states n=8/group)",
        "D9_selective_reporting": "Unclear (control survival not explicitly enumerated; inferred from LD100)",
        "D10_other_bias": "Unclear (sample-size calculation NS)",
        "strain_sex_reported": "Yes (BALB/c, female)",
        "age_reported": "Yes (6-8 wk)",
        "followup_reported": "No (NS)",
        "humane_endpoint_reported": "No (NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Unclear/NS",
        "ethics_approval_number": "No (NS)",
    },
    {
        "publication": "Huang 2016",
        "doi": "10.1038/srep20724",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age stated; weight/baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Unclear (exclusions NS)",
        "D9_selective_reporting": "Low/Unclear (all 4 dose arms reported)",
        "D10_other_bias": "Unclear (sample-size calculation NS; shared control across 4 doses)",
        "strain_sex_reported": "Yes (ICR, female)",
        "age_reported": "Yes (6-8 wk)",
        "followup_reported": "Yes (7 d)",
        "humane_endpoint_reported": "No (NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes",
        "ethics_approval_number": "Yes (SYXK (dian) 2010-0007)",
    },
    {
        "publication": "Akbari 2022",
        "doi": "10.1016/j.intimp.2022.108731",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age/weight stated; baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Unclear (exclusions NS)",
        "D9_selective_reporting": "Low/Unclear",
        "D10_other_bias": "Unclear (sample-size calculation NS)",
        "strain_sex_reported": "Yes (BALB/c, male)",
        "age_reported": "Yes (4-6 wk; 20-25 g)",
        "followup_reported": "Yes (7 d / 168 h)",
        "humane_endpoint_reported": "No (NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes (Shahed University)",
        "ethics_approval_number": "No (NS)",
    },
    {
        "publication": "Chaudhuri 2022",
        "doi": "10.3389/fimmu.2022.933445",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age/weight stated; baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Unclear (exclusions NS)",
        "D9_selective_reporting": "Low/Unclear (all 7 arms reported)",
        "D10_other_bias": "Unclear (sample-size calculation NS; shared control across 7 arms)",
        "strain_sex_reported": "Yes (BALB/c, male)",
        "age_reported": "Yes (6-8 wk; 20-25 g)",
        "followup_reported": "Yes (7 d / 168 h)",
        "humane_endpoint_reported": "Yes (euthanized when 'severely ill'; operational criteria NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes (Shahed University)",
        "ethics_approval_number": "Yes (1399.073)",
    },
    {
        "publication": "Qamsari 2020",
        "doi": "10.3389/fimmu.2020.00158",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age/weight stated; baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Noted attrition (2 mice in loop 11 arm 'inadvertently removed prior to challenge', n=6 vs n=8)",
        "D9_selective_reporting": "Low (explicit x/n printed for all arms)",
        "D10_other_bias": "Unclear (sample-size calculation NS)",
        "strain_sex_reported": "Yes (BALB/c, male)",
        "age_reported": "Yes (6-8 wk; 20-25 g)",
        "followup_reported": "Yes (7 d / 168 h)",
        "humane_endpoint_reported": "Yes (euthanized when 'severely sick'; operational criteria NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes (Shahed University)",
        "ethics_approval_number": "No (NS)",
    },
    {
        "publication": "Ren 2019",
        "doi": "10.4414/smw.2019.20052",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Reported ('randomly divided'; sequence generation method NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age stated; baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Reported (random selection for immunoassay/survival cohorts)",
        "D7_outcome_assessor_blinding": "Partial (blinded for histopathology only; survival objective)",
        "D8_incomplete_outcome_data": "Concern/Discrepancy (72 mice procured vs 63 allocated across cohorts unexplained)",
        "D9_selective_reporting": "Low/Unclear",
        "D10_other_bias": "Unclear (sample-size calculation NS)",
        "strain_sex_reported": "Yes (BALB/c, female)",
        "age_reported": "Yes (6-8 wk)",
        "followup_reported": "Yes (10 d)",
        "humane_endpoint_reported": "No (NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes",
        "ethics_approval_number": "Yes (SCXKCHUAN2018-18)",
    },
    {
        "publication": "Hagag 2022",
        "doi": "10.1007/s00253-022-12231-3",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/weight stated; age & baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Partial (veterinary pathologist blinded for liver histopathology only)",
        "D8_incomplete_outcome_data": "Unclear (exclusions NS)",
        "D9_selective_reporting": "Low/Unclear",
        "D10_other_bias": "Unclear (sample-size calculation NS)",
        "strain_sex_reported": "Yes (BALB/c, female)",
        "age_reported": "No (NS; weight 20-25 g stated)",
        "followup_reported": "Yes (7 d)",
        "humane_endpoint_reported": "No (NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes",
        "ethics_approval_number": "Yes (2021-284)",
    },
    {
        "publication": "Prather 2026",
        "doi": "10.3390/vaccines14030275",
        "appraisal_scope": "Full_10_domain",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Unclear (strain/sex/age stated; baseline comparability NS)",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Unclear (NS)",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Unclear (NS)",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Unclear (exclusions NS)",
        "D9_selective_reporting": "Low/Unclear",
        "D10_other_bias": "Unclear (sample-size calculation NS)",
        "strain_sex_reported": "Yes (C57BL/6 exp 1 & BALB/c exp 2, female)",
        "age_reported": "Yes (7-8 wk)",
        "followup_reported": "Yes (20 d exp 1; 7 d exp 2)",
        "humane_endpoint_reported": "No (NS)",
        "sample_size_calc_reported": "No (NS)",
        "ethics_committee_stated": "Yes (UTSA IACUC)",
        "ethics_approval_number": "Yes (MU070)",
    },
    # 6 Supplementary publications (core items extracted in A._baumannii_mouse_vaccine_challenge_report.md)
    {
        "publication": "Singh 2016",
        "doi": "10.3389/fmicb.2016.00158",
        "appraisal_scope": "Supplementary_core_4_items",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Not appraised in supplementary stage",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Not appraised in supplementary stage",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Not appraised in supplementary stage",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Not appraised in supplementary stage",
        "D9_selective_reporting": "Not appraised in supplementary stage",
        "D10_other_bias": "Not appraised in supplementary stage",
        "strain_sex_reported": "Not extracted in supplementary audit",
        "age_reported": "Not extracted in supplementary audit",
        "followup_reported": "Yes (7 d)",
        "humane_endpoint_reported": "Not extracted in supplementary audit",
        "sample_size_calc_reported": "Not extracted in supplementary audit",
        "ethics_committee_stated": "Yes (Animal Ethics Committee of Panjab University, Chandigarh)",
        "ethics_approval_number": "No (NS)",
    },
    {
        "publication": "Garg 2016",
        "doi": "10.1016/j.ijmm.2015.10.005",
        "appraisal_scope": "Supplementary_core_4_items",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Not appraised in supplementary stage",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Not appraised in supplementary stage",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Not appraised in supplementary stage",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Not appraised in supplementary stage",
        "D9_selective_reporting": "Not appraised in supplementary stage",
        "D10_other_bias": "Not appraised in supplementary stage",
        "strain_sex_reported": "Not extracted in supplementary audit",
        "age_reported": "Not extracted in supplementary audit",
        "followup_reported": "Yes (7 d)",
        "humane_endpoint_reported": "Not extracted in supplementary audit",
        "sample_size_calc_reported": "Not extracted in supplementary audit",
        "ethics_committee_stated": "Yes (Animal Ethics Committee of Panjab University)",
        "ethics_approval_number": "No (NS)",
    },
    {
        "publication": "Skerniskyte 2019",
        "doi": "10.1186/s12866-019-1615-3",
        "appraisal_scope": "Supplementary_core_4_items",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Not appraised in supplementary stage",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Not appraised in supplementary stage",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Not appraised in supplementary stage",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Not appraised in supplementary stage",
        "D9_selective_reporting": "Not appraised in supplementary stage",
        "D10_other_bias": "Not appraised in supplementary stage",
        "strain_sex_reported": "Not extracted in supplementary audit",
        "age_reported": "Not extracted in supplementary audit",
        "followup_reported": "Yes (7 d)",
        "humane_endpoint_reported": "Not extracted in supplementary audit",
        "sample_size_calc_reported": "Not extracted in supplementary audit",
        "ethics_committee_stated": "Yes (Lithuanian State Food and Veterinary Service)",
        "ethics_approval_number": "Yes (G2-72)",
    },
    {
        "publication": "Behrouz 2024",
        "doi": "10.1038/s41598-024-79013-x",
        "appraisal_scope": "Supplementary_core_4_items",
        "D1_sequence_generation": "Unclear (NS; 'allocated to experimental and control groups')",
        "D2_baseline_characteristics": "Not appraised in supplementary stage",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Not appraised in supplementary stage",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Not appraised in supplementary stage",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Not appraised in supplementary stage (active FliC-only arm graph-only)",
        "D9_selective_reporting": "Not appraised in supplementary stage",
        "D10_other_bias": "Not appraised in supplementary stage",
        "strain_sex_reported": "Not extracted in supplementary audit",
        "age_reported": "Not extracted in supplementary audit",
        "followup_reported": "Yes (7 d)",
        "humane_endpoint_reported": "Not extracted in supplementary audit",
        "sample_size_calc_reported": "Not extracted in supplementary audit",
        "ethics_committee_stated": "Yes (Shahed University Research Ethics Committees)",
        "ethics_approval_number": "Yes (R.SHAHED.REC.1400.197)",
    },
    {
        "publication": "Freitas 2024",
        "doi": "10.1016/j.vaccine.2024.04.094",
        "appraisal_scope": "Supplementary_core_4_items",
        "D1_sequence_generation": "Unclear (NS)",
        "D2_baseline_characteristics": "Not appraised in supplementary stage",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Not appraised in supplementary stage",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Not appraised in supplementary stage",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Not appraised in supplementary stage",
        "D9_selective_reporting": "Not appraised in supplementary stage",
        "D10_other_bias": "Not appraised in supplementary stage",
        "strain_sex_reported": "Yes (BALB/c, female)",
        "age_reported": "Not extracted in supplementary audit",
        "followup_reported": "Yes (7 d)",
        "humane_endpoint_reported": "Not extracted in supplementary audit",
        "sample_size_calc_reported": "Not extracted in supplementary audit",
        "ethics_committee_stated": "Yes (UFPel CEUA)",
        "ethics_approval_number": "Yes (007034/2021-16)",
    },
    {
        "publication": "Jeffreys 2024",
        "doi": "10.3390/vaccines12040358",
        "appraisal_scope": "Supplementary_core_4_items",
        "D1_sequence_generation": "Reported in main experiment ('randomly divided'; sequence method NS)",
        "D2_baseline_characteristics": "Not appraised in supplementary stage",
        "D3_allocation_concealment": "Unclear (NS)",
        "D4_random_housing": "Not appraised in supplementary stage",
        "D5_caregiver_blinding": "Unclear (NS)",
        "D6_random_outcome_assessment": "Not appraised in supplementary stage",
        "D7_outcome_assessor_blinding": "Unclear (NS)",
        "D8_incomplete_outcome_data": "Not appraised in supplementary stage",
        "D9_selective_reporting": "Not appraised in supplementary stage (pilot IP mock n NS)",
        "D10_other_bias": "Not appraised in supplementary stage",
        "strain_sex_reported": "Strain Yes in PMC (C57BL/6); sex NS in supp audit",
        "age_reported": "Yes in PMC (7-8 wk)",
        "followup_reported": "Yes (30 d)",
        "humane_endpoint_reported": "Not extracted in supplementary audit",
        "sample_size_calc_reported": "Not extracted in supplementary audit",
        "ethics_committee_stated": "Yes (UTSA IACUC)",
        "ethics_approval_number": "Yes (MU070)",
    },
]


# =============================================================================
# 8. COMPLETE TABLE S2: EXCLUDED STUDIES (8 ORIGINAL + 28 SUPPLEMENTARY SURVIVAL + 41 OTHER SCREENED)
# =============================================================================

EXCLUDED_STUDIES: list[dict[str, str]] = [
    # --- Group A: 8 Excluded Original Track 1 Publications (Re-extraction Audit) ---
    {
        "record_id": "EX_ORIG_01",
        "study": "Luo et al., 2012",
        "doi": "10.1371/journal.pone.0029446",
        "antigen_or_platform": "rOmpA subunit (diabetic mouse model)",
        "exclusion_stage_group": "Original_Track1_Reextraction_8",
        "exclusion_category": "Graph-only survival (no printed x/n or %)",
        "reason_detail": "Survival x/n counts and percentages not stated in text/tables; survival curves only.",
    },
    {
        "record_id": "EX_ORIG_02",
        "study": "Lin et al., 2013",
        "doi": "10.1016/j.vaccine.2012.11.008",
        "antigen_or_platform": "rOmpA dose/immune polarization",
        "exclusion_stage_group": "Original_Track1_Reextraction_8",
        "exclusion_category": "No in vivo bacterial challenge survival experiment",
        "reason_detail": "Evaluated rOmpA dose-dependent immune polarization and epitope mapping only; no bacterial challenge arm.",
    },
    {
        "record_id": "EX_ORIG_03",
        "study": "Zhang et al., 2016",
        "doi": "10.1016/j.micpath.2016.04.019",
        "antigen_or_platform": "Purified OmpA (intranasal mucosal vaccination)",
        "exclusion_stage_group": "Original_Track1_Reextraction_8",
        "exclusion_category": "Graph-only / qualitative survival (no numerical x/n or %)",
        "reason_detail": "No numerical survival counts or percentages stated in text; reports only 'higher survival, p < 0.05' with survival curve.",
    },
    {
        "record_id": "EX_ORIG_04",
        "study": "Fereshteh et al., 2022",
        "doi": "10.1016/j.micpath.2022.105882",
        "antigen_or_platform": "Recombinant DcaP-like protein and AbOmpA",
        "exclusion_stage_group": "Original_Track1_Reextraction_8",
        "exclusion_category": "Internally inconsistent survival counts/percentages",
        "reason_detail": "States n=8/group yet reports 80%/100% survival alongside 'one death per group' (7/8 = 87.5%, incompatible with 80%).",
    },
    {
        "record_id": "EX_ORIG_05",
        "study": "Hu et al., 2023",
        "doi": "10.1016/j.micpath.2023.106195",
        "antigen_or_platform": "Unverified citation (putative OmpA/Omp22/OmpK)",
        "exclusion_stage_group": "Original_Track1_Reextraction_8",
        "exclusion_category": "Unverifiable / mismatched bibliographic record",
        "reason_detail": "Citation could not be matched to an eligible mouse challenge trial (listed DOI resolves to an unrelated review/article).",
    },
    {
        "record_id": "EX_ORIG_06",
        "study": "Zhu et al., 2025",
        "doi": "10.1002/eji.70019",
        "antigen_or_platform": "Immunogenicity-guided subunit candidate",
        "exclusion_stage_group": "Original_Track1_Reextraction_8",
        "exclusion_category": "Preprint only / peer-reviewed full text not verified",
        "reason_detail": "Only a preprint record could be verified at audit cutoff, not the peer-reviewed published article.",
    },
    {
        "record_id": "EX_ORIG_07",
        "study": "Huang et al., 2014",
        "doi": "10.1371/journal.pone.0100727",
        "antigen_or_platform": "A. baumannii outer membrane vesicles (OMVs)",
        "exclusion_stage_group": "Original_Track1_Reextraction_8",
        "exclusion_category": "Non-subunit vaccine platform (OMV)",
        "reason_detail": "Evaluated native A. baumannii OMV vaccine rather than a defined recombinant protein or peptide subunit.",
    },
    {
        "record_id": "EX_ORIG_08",
        "study": "Du et al., 2021",
        "doi": "10.2147/IJN.S296527",
        "antigen_or_platform": "Chitosan-PLGA-rOmp22 multi-epitope peptide",
        "exclusion_stage_group": "Original_Track1_Reextraction_8",
        "exclusion_category": "Percentages without stated survival group size (n)",
        "reason_detail": "Survival percentages reported without stating survival-cohort group sizes (denominators).",
    },

    # --- Group B: Exact 28 Unique Supplementary Protein/Peptide Survival-Reporting Records Excluded ---
    {
        "record_id": "EX_SUPP_SURV_01",
        "study": "Guo et al., 2018",
        "doi": "10.3967/bes2018.019",
        "antigen_or_platform": "OmpK, Omp22, and fused OmpK/Omp22",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Percentages without stated survival group size (n)",
        "reason_detail": "Reports 66.7% (fusion), 25% (OmpK), 37.5% (Omp22) vs 0% control at day 8, but survival-cohort n is not stated.",
    },
    {
        "record_id": "EX_SUPP_SURV_02",
        "study": "Yang et al., 2019",
        "doi": "10.1111/sji.12769",
        "antigen_or_platform": "MF59-adjuvanted OmpK/Omp22 fusion (intratracheal)",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Missing control survival count/% (abstract-only accessible)",
        "reason_detail": "Accessible abstract states n=12/group and 83.3% (MF59) vs 41.7% (fusion alone), but control survival is not stated.",
    },
    {
        "record_id": "EX_SUPP_SURV_03",
        "study": "Ramezanalizadeh et al., 2020 (Vaccine)",
        "doi": "10.1016/j.vaccine.2020.06.052",
        "antigen_or_platform": "Recombinant pilus proteins CsuA/B and FimA",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Missing control survival count/% (abstract-only accessible)",
        "reason_detail": "Accessible abstract gives active survival 5/8, 3/8, 4/8, but control group survival and n are not reported.",
    },
    {
        "record_id": "EX_SUPP_SURV_04",
        "study": "Ramezanalizadeh et al., 2022",
        "doi": "10.1038/s41598-022-24163-z",
        "antigen_or_platform": "Quadrivalent CsuA/B + FimA + BauA + HemTR",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Graph-only control survival",
        "reason_detail": "Active survival stated as 87.5% (n=8/group) for ABI038 and ABI022, but control survival is plotted in Fig. 7 without numeric text/legend value.",
    },
    {
        "record_id": "EX_SUPP_SURV_05",
        "study": "Alipouri et al., 2022",
        "doi": "10.1186/s43141-022-00325-4",
        "antigen_or_platform": "Recombinant VgrG",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Percentages without dose-specific group size (n)",
        "reason_detail": "12 vaccinated and 12 control mice total split across 3 challenge doses (1.9, 2.4, 3.2 x 10^8 CFU); dose-specific n not stated (72-h follow-up).",
    },
    {
        "record_id": "EX_SUPP_SURV_06",
        "study": "Masoumkhani et al., 2020",
        "doi": "10.15761/TiM.1000236",
        "antigen_or_platform": "Recombinant VacJ lipoprotein",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Qualitative survival only / control survival NR",
        "reason_detail": "n=6/group; IN challenge survival reported qualitatively without count/%, IN control survival NR, IP challenge had 0% survival.",
    },
    {
        "record_id": "EX_SUPP_SURV_07",
        "study": "Ramezanalizadeh et al., 2021",
        "doi": "10.2217/fmb-2020-0133",
        "antigen_or_platform": "HemTR and BauA siderophore receptors",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Incomplete abstract-only record (no numerical survival)",
        "reason_detail": "Accessible abstract reports delayed mortality and better combination protection without numerical group counts or n.",
    },
    {
        "record_id": "EX_SUPP_SURV_08",
        "study": "Li et al., 2017",
        "doi": "10.3892/mmr.2017.6688",
        "antigen_or_platform": "His-SmpA and His-PLD",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Non-lethal active challenge (0 deaths in any group at 72 h) / separate survival n NR",
        "reason_detail": "No mice in any active or control group died within 72 h; 15-day lethal challenge survival was conducted only for passive serum transfer.",
    },
    {
        "record_id": "EX_SUPP_SURV_09",
        "study": "Sefidi et al., 2016",
        "doi": "10.5662/wjm.v6.i3.190",
        "antigen_or_platform": "Bap with or without Pseudomonas flagellin",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Percentages without challenge-arm n and control survival NR",
        "reason_detail": "Reports 100% (Bap) and 80% (Bap-flagellin) at 48 h, but A. baumannii challenge-arm n and control survival are not stated.",
    },
    {
        "record_id": "EX_SUPP_SURV_10",
        "study": "Esmaeilkhani et al., 2019",
        "doi": "NR (Avicenna J Med Biotechnol 2019;11:180-186)",
        "antigen_or_platform": "BauA cork, loop, and cork-loop cocktail",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Qualitative survival only / 48-h follow-up",
        "reason_detail": "n=5/group across multiple routes and challenge doses; text states immunized mice 'survived in all groups' without formulation-specific counts (48-h follow-up).",
    },
    {
        "record_id": "EX_SUPP_SURV_11",
        "study": "Song et al., 2020",
        "doi": "10.1155/2020/1975875",
        "antigen_or_platform": "MacB-derived B/CD8/CD4 epitopes (RAE)",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Graph-only per-group survival",
        "reason_detail": "n=10/challenge group; per-group endpoint survival is available only in Fig. 10 (text states only all PBS died within 6 d and RAE died within 14 d).",
    },
    {
        "record_id": "EX_SUPP_SURV_12",
        "study": "Sun et al., 2022",
        "doi": "10.3389/fimmu.2022.884555",
        "antigen_or_platform": "39-aa Ata-alpha peptide fused to CTB",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Percentages without stated survival group size (n)",
        "reason_detail": "Reports 90% survival for CTB-Ata-alpha + alum vs 0% controls within 2 days, but group sizes (n) are not reported.",
    },
    {
        "record_id": "EX_SUPP_SURV_13",
        "study": "Mirali et al., 2023",
        "doi": "10.1016/j.micpath.2023.106262",
        "antigen_or_platform": "Recombinant Omp34 and BauA cocktail",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Internally inconsistent control survival percentage vs group size",
        "reason_detail": "Methods state n=10/group, yet Results report active control survival as 42.8% (non-integer for n=10; corresponds to 3/7, unresolved).",
    },
    {
        "record_id": "EX_SUPP_SURV_14",
        "study": "Raoufi & Abdollahi, 2024",
        "doi": "10.1016/j.intimp.2024.112972",
        "antigen_or_platform": "OprB porin and epitope constructs",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Unreconciled group size (n=5 vs n=8) and control survival NR numerically",
        "reason_detail": "Immunization states n=5/group while challenge states n=8/group; active % given at 2xLD100, but control survival is not numerically stated.",
    },
    {
        "record_id": "EX_SUPP_SURV_15",
        "study": "Raoufi et al., 2026",
        "doi": "10.1016/j.ijbiomac.2026.155110",
        "antigen_or_platform": "TonB-dependent cobalamin (vitamin B12) receptor and loop constructs",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Group-level ranges only / internally inconsistent with n=5",
        "reason_detail": "Reports combined ranges only (60-100% at LD100; 50-75% text vs 50-70% legend at 2xLD100) without formulation-specific values; 50%/75% incompatible with n=5.",
    },
    {
        "record_id": "EX_SUPP_SURV_16",
        "study": "Chiang et al., 2015",
        "doi": "10.1080/21645515.2015.1010910",
        "antigen_or_platform": "Three reverse-vaccinology recombinant outer membrane proteins",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Percentages without stated group size / incomplete abstract-only record",
        "reason_detail": "Accessible record reports 60% partial protection in pneumonia model without n, control survival, or per-antigen counts.",
    },
    {
        "record_id": "EX_SUPP_SURV_17",
        "study": "Sangroodi et al., 2015",
        "doi": "10.3906/SAG-1407-45",
        "antigen_or_platform": "Recombinant BauA cork and beta-barrel domains",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Incomplete abstract-only record (no n or numerical counts)",
        "reason_detail": "Accessible record states immunized mice survived and controls died after 24 h, without group sizes, challenge dose, or per-domain counts.",
    },
    {
        "record_id": "EX_SUPP_SURV_18",
        "study": "Cai et al., 2016",
        "doi": "NR (Immunological Journal 2016:104-108)",
        "antigen_or_platform": "Recombinant A1S1969 outer membrane protein",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Percentages without stated survival group size (n)",
        "reason_detail": "Reports 56.6% vaccinated vs 20% control survival, but group sizes, adjuvant, vaccination route, and follow-up are not stated.",
    },
    {
        "record_id": "EX_SUPP_SURV_19",
        "study": "Aghajani et al., 2019",
        "doi": "10.1111/apm.12992",
        "antigen_or_platform": "Siderophore receptors BfnH and BauA",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Incomplete abstract-only record (active survival qualitative only)",
        "reason_detail": "Active groups described qualitatively ('partial protection and increased survival') without % or n; 50% values were passive transfer.",
    },
    {
        "record_id": "EX_SUPP_SURV_20",
        "study": "Mahmoudi et al., 2020",
        "doi": "10.1111/apm.13061",
        "antigen_or_platform": "Conserved recombinant fimbrial protein",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Incomplete abstract-only record (no n or control survival)",
        "reason_detail": "Abstract states all actively immunized mice survived ATCC 19606 challenge, without n, control survival, route, or follow-up.",
    },
    {
        "record_id": "EX_SUPP_SURV_21",
        "study": "Pazoki et al., 2021",
        "doi": "10.1016/j.molimm.2021.04.026",
        "antigen_or_platform": "VgrG and LCL-scaffolded VgrG segment",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Incomplete abstract-only record (no numerical outcomes or n)",
        "reason_detail": "Abstract states survival increased and LCL-VgrG protected better, without numerical counts, percentages, or group sizes.",
    },
    {
        "record_id": "EX_SUPP_SURV_22",
        "study": "Raoufi et al., 2021",
        "doi": "10.1016/j.micpath.2021.105346",
        "antigen_or_platform": "DcaP porin and epitope-based subunit",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Incomplete abstract-only record (no numerical outcomes or n)",
        "reason_detail": "Abstract states mortality decreased without exact survival values, group sizes, or control counts.",
    },
    {
        "record_id": "EX_SUPP_SURV_23",
        "study": "Zeng et al., 2023",
        "doi": "10.1016/j.ygeno.2023.110590",
        "antigen_or_platform": "Recombinant peptidoglycan-associated lipoprotein (Pal)",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Percentages without stated group size / control survival NR",
        "reason_detail": "Accessible record reports 40% protection in pneumonia model without n, control survival, or regimen details.",
    },
    {
        "record_id": "EX_SUPP_SURV_24",
        "study": "Erami et al., 2021",
        "doi": "10.1016/j.micpath.2021.105291",
        "antigen_or_platform": "Recombinant Omp34",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Incomplete abstract-only record (no group survival values)",
        "reason_detail": "Survival at 5xLD50 listed as an endpoint in abstract, but group survival counts/percentages and n are not reported.",
    },
    {
        "record_id": "EX_SUPP_SURV_25",
        "study": "Bolourchi et al., 2019",
        "doi": "10.1016/j.micpath.2019.03.031",
        "antigen_or_platform": "Recombinant PKF and AbOmpA",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Incomplete abstract-only record (no group survival values)",
        "reason_detail": "Survival rates studied in sepsis model, but accessible abstract gives no numerical values or group sizes.",
    },
    {
        "record_id": "EX_SUPP_SURV_26",
        "study": "Hosseingholi et al., 2015",
        "doi": "NR (BioTechnology: An Indian Journal Vol. 11)",
        "antigen_or_platform": "Recombinant phospholipase D (PLD)",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Percentages without stated n / unverified peer-review status",
        "reason_detail": "Reports 95% active and 80% passive survival without group size n; journal indexing/peer-review status unverified.",
    },
    {
        "record_id": "EX_SUPP_SURV_27",
        "study": "Fashami et al., 2025",
        "doi": "10.1016/j.medmic.2025.100164",
        "antigen_or_platform": "Intranasal rOmp34 and rOmp34L3x5",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Incomplete abstract-only record (no survival values in abstract)",
        "reason_detail": "Abstract states survival analysis was performed after Ab/TU/ColR pulmonary challenge, without outcomes or group sizes.",
    },
    {
        "record_id": "EX_SUPP_SURV_28",
        "study": "Mehdinejadiani et al., 2021",
        "doi": "10.18502/ijaai.v20i5.7404",
        "antigen_or_platform": "OmpA epitope-based peptide vaccines",
        "exclusion_stage_group": "Supplementary_Protein_Peptide_Survival_28",
        "exclusion_category": "Graph-only per-group survival",
        "reason_detail": "n=5 per active group, but per-group survival proportions are plotted only in Fig. 8 without printed x/n or % in prose.",
    },

    # --- Group C: Additional Screened Supplementary Records Excluded by Platform or Non-Survival Endpoint ---
    {"record_id": "EX_PLAT_01", "study": "Huang et al., 2016 (OMV-Omp22)", "doi": "10.1038/srep37242", "antigen_or_platform": "E. coli-derived OMVs displaying Omp22", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (OMV/OMC)", "reason_detail": "OMV delivery platform rather than purified protein/peptide subunit alone."},
    {"record_id": "EX_PLAT_02", "study": "Pulido et al., 2019", "doi": "10.1016/j.vaccine.2019.11.043", "antigen_or_platform": "LPS-free OMV vaccine", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (OMV/OMC)", "reason_detail": "OMV vaccine platform."},
    {"record_id": "EX_PLAT_03", "study": "Li et al., 2020", "doi": "10.3389/fimmu.2020.01069", "antigen_or_platform": "A. baumannii OMV vaccine", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (OMV/OMC)", "reason_detail": "OMV vaccine platform."},
    {"record_id": "EX_PLAT_04", "study": "Higham et al., 2023", "doi": "10.1016/j.jinf.2023.02.035", "antigen_or_platform": "Intranasal OMV vaccine", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (OMV/OMC)", "reason_detail": "OMV vaccine platform."},
    {"record_id": "EX_PLAT_05", "study": "Bjanes et al., 2022", "doi": "10.1002/anbr.202200130", "antigen_or_platform": "OMV-coated nanoparticle vaccine", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (OMV/OMC)", "reason_detail": "OMV-coated nanoparticle platform."},
    {"record_id": "EX_PLAT_06", "study": "Pulido et al., 2018", "doi": "10.1016/j.vaccine.2018.05.113", "antigen_or_platform": "LPS-free outer membrane complex (OMC)", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (OMV/OMC)", "reason_detail": "Undefined outer membrane complex rather than recombinant subunit."},
    {"record_id": "EX_PLAT_07", "study": "McConnell et al., 2010", "doi": "10.1128/IAI.00741-10", "antigen_or_platform": "Outer membrane complex (OMC)", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (OMV/OMC)", "reason_detail": "Outer membrane complex platform."},
    {"record_id": "EX_PLAT_08", "study": "Lei et al., 2019", "doi": "10.1007/s11033-019-04994-2", "antigen_or_platform": "OmpA + Pal DNA vaccine", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (DNA/mRNA)", "reason_detail": "Nucleic-acid DNA vaccine platform."},
    {"record_id": "EX_PLAT_09", "study": "Ansari et al., 2019", "doi": "10.22038/ijbms.2019.30799.7427", "antigen_or_platform": "ompA DNA vaccine", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (DNA/mRNA)", "reason_detail": "DNA vaccine platform (immunogenicity only)."},
    {"record_id": "EX_PLAT_10", "study": "Hosseinnezhad-Lazarjani et al., 2023", "doi": "10.1016/j.vaccine.2023.02.046", "antigen_or_platform": "Chitosan-csuC DNA nanovaccine", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (DNA/mRNA)", "reason_detail": "DNA nanovaccine platform."},
    {"record_id": "EX_PLAT_11", "study": "Liu et al., 2025", "doi": "10.3390/vaccines13070764", "antigen_or_platform": "SmpA-PLD and OmpK-Omp22 mRNA-LNP", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (DNA/mRNA)", "reason_detail": "mRNA-LNP vaccine platform."},
    {"record_id": "EX_PLAT_12", "study": "Higham et al., 2025", "doi": "10.1038/s41541-025-01202-0", "antigen_or_platform": "Intranasal mRNA vaccine", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (DNA/mRNA)", "reason_detail": "mRNA platform and bacterial burden endpoint."},
    {"record_id": "EX_PLAT_13", "study": "Ma et al., 2025", "doi": "10.1038/s41598-024-84823-0", "antigen_or_platform": "Multi-epitope mRNA design (in silico)", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "In silico mRNA design (Track 2 only, no mouse challenge)", "reason_detail": "Purely computational mRNA vaccine design (included in Track 2)."},
    {"record_id": "EX_PLAT_14", "study": "Zou et al., 2026", "doi": "10.3390/vaccines14110921", "antigen_or_platform": "KL49 capsular polysaccharide glycoconjugate", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (Polysaccharide conjugate)", "reason_detail": "Capsular polysaccharide glycoconjugate."},
    {"record_id": "EX_PLAT_15", "study": "Rudenko et al., 2022", "doi": "10.1128/spectrum.01674-22", "antigen_or_platform": "K9 capsular oligosaccharide conjugate", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (Polysaccharide conjugate)", "reason_detail": "Capsular polysaccharide conjugate."},
    {"record_id": "EX_PLAT_16", "study": "Wei et al., 2021", "doi": "10.1021/acscentsci.1c00656", "antigen_or_platform": "Pseudaminic acid-CRM197 conjugate", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (Polysaccharide conjugate)", "reason_detail": "Carbohydrate-CRM197 glycoconjugate."},
    {"record_id": "EX_PLAT_17", "study": "Lee et al., 2018", "doi": "10.1021/jacs.8b04078", "antigen_or_platform": "Exopolysaccharide glycoconjugate", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (Polysaccharide conjugate)", "reason_detail": "Exopolysaccharide glycoconjugate."},
    {"record_id": "EX_PLAT_18", "study": "Sianturi et al., 2022", "doi": "10.1002/anie.202209556", "antigen_or_platform": "Semisynthetic glycoconjugate", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (Polysaccharide conjugate)", "reason_detail": "Semisynthetic glycoconjugate without mouse survival challenge."},
    {"record_id": "EX_PLAT_19", "study": "Bentancor et al., 2011", "doi": "10.1128/IAI.05653-11", "antigen_or_platform": "PNAG-tetanus toxoid conjugate / passive antibody", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (Polysaccharide conjugate / passive)", "reason_detail": "Polysaccharide conjugate and passive antibody CFU endpoints."},
    {"record_id": "EX_PLAT_20", "study": "Dollery et al., 2021", "doi": "10.3390/vaccines9020096", "antigen_or_platform": "Radiation-inactivated whole-cell vaccine", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (Whole-cell / live-attenuated)", "reason_detail": "Whole-cell inactivated bacterin."},
    {"record_id": "EX_PLAT_21", "study": "Dollery et al., 2022", "doi": "10.3389/fimmu.2022.941010", "antigen_or_platform": "Gamma- or UV-inactivated whole-cell vaccine", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (Whole-cell / live-attenuated)", "reason_detail": "Whole-cell inactivated bacterin."},
    {"record_id": "EX_PLAT_22", "study": "Khan et al., 2022", "doi": "10.3390/pharmaceutics14071357", "antigen_or_platform": "Liposome-encapsulated whole-cell antigens", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (Whole-cell / live-attenuated)", "reason_detail": "Whole-cell antigen extract."},
    {"record_id": "EX_PLAT_23", "study": "Shu et al., 2016", "doi": "10.1038/srep22332", "antigen_or_platform": "Antibiotic-exposed inactivated whole-cell vaccine", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (Whole-cell / live-attenuated)", "reason_detail": "Inactivated whole-cell vaccine; in vitro bactericidal assays only."},
    {"record_id": "EX_PLAT_24", "study": "Sharmin et al., 2021", "doi": "10.1101/2021.08.30.458307", "antigen_or_platform": "Heat-inactivated whole-cell vaccine (bioRxiv)", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform & preprint", "reason_detail": "Whole-cell bacterin and unrefereed preprint."},
    {"record_id": "EX_PLAT_25", "study": "Bakht et al., 2025", "doi": "10.1016/j.vaccine.2025.127298", "antigen_or_platform": "Live-attenuated delta-hfq mutant", "exclusion_stage_group": "Supplementary_Platform_Exclusion_25", "exclusion_category": "Non-subunit vaccine platform (Whole-cell / live-attenuated)", "reason_detail": "Live-attenuated mutant vaccine."},

    # --- Group D: 16 Screened Subunit/Candidate Studies Without Verified Active Survival Endpoint ---
    {"record_id": "EX_NOSURV_01", "study": "Ni et al., 2017", "doi": "10.3390/ijms18020458", "antigen_or_platform": "Reverse-vaccinology resistance determinants", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "No in vivo mouse survival experiment", "reason_detail": "In silico reverse-vaccinology prediction only."},
    {"record_id": "EX_NOSURV_02", "study": "Bentancor et al., 2012", "doi": "10.1128/IAI.06096-11", "antigen_or_platform": "Trimeric autotransporter Ata (passive antisera)", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "Passive immunization / bacterial burden only", "reason_detail": "Tested anti-Ata antibodies and lung CFU burden rather than active vaccination survival."},
    {"record_id": "EX_NOSURV_03", "study": "Fattahian et al., 2011", "doi": "10.1016/j.micpath.2011.09.004", "antigen_or_platform": "Biofilm-associated protein (Bap)", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "Bacterial burden endpoint only ( check/no verified survival counts )", "reason_detail": "Reported bacterial burden outcomes rather than usable survival counts in audit record."},
    {"record_id": "EX_NOSURV_04", "study": "Gholizadeh et al., 2022", "doi": "10.4103/1735-5362.350237", "antigen_or_platform": "OMP-PLGA nanoparticles", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "Bacterial burden endpoint only (no survival)", "reason_detail": "Mice euthanized for organ CFU enumeration; no survival endpoint."},
    {"record_id": "EX_NOSURV_05", "study": "Timm et al., 2025", "doi": "10.1038/s41467-025-62402-9", "antigen_or_platform": "Adhesin Abp2D (CAUTI model)", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "Bacterial burden endpoint only (no survival)", "reason_detail": "Non-lethal CAUTI model measuring bladder/catheter bacterial titers at 24 h."},
    {"record_id": "EX_NOSURV_06", "study": "Chiang et al., 2020", "doi": "10.3389/fmicb.2020.00536", "antigen_or_platform": "rAdeK + tigecycline combination", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "Bacterial burden endpoint only (no survival)", "reason_detail": "Post-challenge endpoints were tissue bacterial counts with tigecycline, not survival."},
    {"record_id": "EX_NOSURV_07", "study": "Palethorpe et al., 2026", "doi": "10.1371/journal.ppat.1013958", "antigen_or_platform": "Microarray-identified recombinant protein antigens (Ag7, etc.)", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "Bacterial burden endpoint only (no survival)", "reason_detail": "Active vaccination and passive mAb experiments measured spleen/lung CFU burden rather than survival."},
    {"record_id": "EX_NOSURV_08", "study": "Tamehri et al., 2022", "doi": "10.1016/j.micpath.2022.105874", "antigen_or_platform": "BauA + OmpA combination", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "Bacterial burden / histopathology only (no verified survival)", "reason_detail": "Reported organ bacterial burdens and histopathology rather than verified survival counts."},
    {"record_id": "EX_NOSURV_09", "study": "Abdallah et al., 2023", "doi": "10.21608/aprh.2023.177396.1202", "antigen_or_platform": "ABC transporter protein antigen", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "No verified survival endpoint in accessible record", "reason_detail": "Abstract mentions protection lasting 24 h without establishing lethal survival counts."},
    {"record_id": "EX_NOSURV_10", "study": "Sabzi et al., 2024", "doi": "10.1016/j.ijpharm.2024.123961", "antigen_or_platform": "rOmp22-polydopamine nanoadjuvant", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "No verified survival endpoint in accessible record", "reason_detail": "Accessible abstract does not report survival counts."},
    {"record_id": "EX_NOSURV_11", "study": "Yang et al., 2025", "doi": "10.1021/acsinfecdis.4c01049", "antigen_or_platform": "Intranasal chitosan-PLGA-rOmp22", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "No verified survival endpoint in accessible record", "reason_detail": "Accessible abstract describes pneumonia protection without survival counts/n."},
    {"record_id": "EX_NOSURV_12", "study": "Davoudi et al., 2021", "doi": "10.18502/ijm.v13i5.7429", "antigen_or_platform": "Recombinant multi-epitope protein", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "No in vivo mouse survival experiment", "reason_detail": "Expression and antibody recognition only; no mouse challenge survival."},
    {"record_id": "EX_NOSURV_13", "study": "Mehdinejadiani et al., 2019", "doi": "10.18502/ijaai.v18i6.2178", "antigen_or_platform": "OmpA antigenic peptides", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "No in vivo mouse survival experiment", "reason_detail": "Peptide immunogenicity in mice; no lethal bacterial challenge."},
    {"record_id": "EX_NOSURV_14", "study": "Banisaeed et al., 2026", "doi": "10.1016/j.vaccine.2026.129205", "antigen_or_platform": "LCL-scaffolded OmpA loop 3", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "No verified survival endpoint in accessible record", "reason_detail": "Abstract mentions protection without numerical survival outcomes; full text unavailable."},
    {"record_id": "EX_NOSURV_15", "study": "Lu et al., 2026", "doi": "10.1039/d6bm00143b", "antigen_or_platform": "Biodegradable polymeric nanoparticle vaccine", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "Unresolved vaccine platform & missing n", "reason_detail": "Abstract reports 100% survival without identifying whether antigen is a defined protein/peptide subunit or stating n."},
    {"record_id": "EX_NOSURV_16", "study": "Hamzavi et al., 2026 / Bashiri et al., 2025", "doi": "10.1038/s41541-026-01609-3; 10.1007/s00284-025-04275-1", "antigen_or_platform": "Omp34 loop-based immunogens", "exclusion_stage_group": "Supplementary_No_Verified_Survival_16", "exclusion_category": "No accessible abstract/full text at audit cutoff", "reason_detail": "Indexed titles only without accessible abstract or full text to verify survival reporting."},
]


# =============================================================================
# 9. WRITE ALL CSVs AND JSON
# =============================================================================

def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


write_csv(OUT_DIR / "track1_experiments_22.csv", exp_rows)
write_csv(OUT_DIR / "track1_arms_all.csv", arm_rows)
write_csv(OUT_DIR / "track1_platform_units_30.csv", platform_units_30)
write_csv(OUT_DIR / "track1_leave_one_out_20.csv", loo_pub_rows)
write_csv(OUT_DIR / "track1_syrcle_20.csv", SYRCLE_20)
write_csv(OUT_DIR / "table_s2_excluded_studies.csv", EXCLUDED_STUDIES)

# SYRCLE exact tallies
syrcle_14 = [s for s in SYRCLE_20 if s["appraisal_scope"] == "Full_10_domain"]
syrcle_counts = {
    "total_publications": len(SYRCLE_20),
    "full_10_domain_publications": len(syrcle_14),
    "supplementary_core_publications": len(SYRCLE_20) - len(syrcle_14),
    "random_allocation_stated_of_20": sum(1 for s in SYRCLE_20 if s["D1_sequence_generation"].startswith("Reported")),
    "sequence_generation_method_stated_of_20": 0,
    "allocation_concealment_stated_of_20": 0,
    "caregiver_blinding_stated_of_20": 0,
    "random_housing_stated_of_14": 0,
    "baseline_comparability_stated_of_14": 0,
    "sample_size_calculation_stated_of_14": 0,
    "random_outcome_selection_stated_of_14": sum(1 for s in syrcle_14 if s["D6_random_outcome_assessment"].startswith("Reported")),
    "blinded_outcome_assessment_histopath_only_of_20": sum(1 for s in SYRCLE_20 if s["D7_outcome_assessor_blinding"].startswith("Partial")),
    "unexplained_exclusions_or_n_discrepancies_of_14": sum(1 for s in syrcle_14 if ("Discrepancy" in s["D8_incomplete_outcome_data"] or "attrition" in s["D8_incomplete_outcome_data"])),
    "strain_and_sex_stated_of_14": sum(1 for s in syrcle_14 if s["strain_sex_reported"].startswith("Yes")),
    "age_stated_of_14": sum(1 for s in syrcle_14 if s["age_reported"].startswith("Yes")),
    "followup_stated_of_14": sum(1 for s in syrcle_14 if s["followup_reported"].startswith("Yes")),
    "followup_stated_of_20": sum(1 for s in SYRCLE_20 if s["followup_reported"].startswith("Yes")),
    "humane_endpoint_stated_of_14": sum(1 for s in syrcle_14 if s["humane_endpoint_reported"].startswith("Yes")),
    "ethics_approval_number_stated_of_20": sum(1 for s in SYRCLE_20 if s["ethics_approval_number"].startswith("Yes")),
    "ethics_committee_or_approval_mentioned_of_20": sum(1 for s in SYRCLE_20 if s["ethics_committee_stated"].startswith("Yes")),
    "ethics_neither_committee_nor_number_in_audit_of_20": [
        s["publication"] for s in SYRCLE_20 if not s["ethics_committee_stated"].startswith("Yes")
    ],
}

stats_payload = {
    "dataset_totals": {
        "publications": len(pub_order),
        "independent_experiments": len(exp_rows),
        "individual_vaccine_arms": len(arm_rows),
        "platform_units": len(platform_units_30),
        "total_mice": sum(r["nv"] + r["nc"] for r in exp_rows),
        "vaccinated_mice": sum(r["nv"] for r in exp_rows),
        "vaccinated_survivors": sum(r["ev"] for r in exp_rows),
        "vaccinated_survival_pct": sum(r["ev"] for r in exp_rows) / sum(r["nv"] for r in exp_rows) * 100.0,
        "control_mice": sum(r["nc"] for r in exp_rows),
        "control_survivors": sum(r["ec"] for r in exp_rows),
        "control_survival_pct": sum(r["ec"] for r in exp_rows) / sum(r["nc"] for r in exp_rows) * 100.0,
        "experiments_with_zero_cell_2x2": sum(1 for r in exp_rows if r["has_zero_cell"]),
        "experiments_with_zero_control_survivors": sum(1 for r in exp_rows if r["ctrl_zero_events"]),
        "experiments_with_surviving_controls": [r["publication"] for r in exp_rows if not r["ctrl_zero_events"]],
        "var_lnRR_range_in_ec0_experiments": [min(var_ec0), max(var_ec0)],
        "var_lnRR_range_in_ec_pos_experiments": [min(var_ec_pos), max(var_ec_pos)],
        "study_RD_raw_pp_range": [min(r["RD_raw_pp"] for r in exp_rows), max(r["RD_raw_pp"] for r in exp_rows)],
        "study_RD_haldane_pp_range": [min(r["RD_haldane_pp"] for r in exp_rows), max(r["RD_haldane_pp"] for r in exp_rows)],
        "study_total_n_range": [min(r["nv"] + r["nc"] for r in exp_rows), max(r["nv"] + r["nc"] for r in exp_rows)],
        "study_vacc_n_range": [min(r["nv"] for r in exp_rows), max(r["nv"] for r in exp_rows)],
    },
    "overall_RR_k22": overall_rr,
    "overall_RD_haldane_both_k22": {
        **overall_rd_cc,
        "est_dl_pp": overall_rd_cc["est_dl"] * 100.0,
        "wald_low_pp": overall_rd_cc["wald_low"] * 100.0,
        "wald_high_pp": overall_rd_cc["wald_high"] * 100.0,
        "hksj_low_pp": overall_rd_cc["hksj_low"] * 100.0,
        "hksj_high_pp": overall_rd_cc["hksj_high"] * 100.0,
        "pi_df20_low_pp": overall_rd_cc["pi_df_k2_low"] * 100.0,
        "pi_df20_high_pp": overall_rd_cc["pi_df_k2_high"] * 100.0,
        "pi_df21_low_pp": overall_rd_cc["pi_df_k1_low"] * 100.0,
        "pi_df21_high_pp": overall_rd_cc["pi_df_k1_high"] * 100.0,
    },
    "overall_RD_raw_pt_haldane_var_k22": {
        **overall_rd_raw,
        "est_dl_pp": overall_rd_raw["est_dl"] * 100.0,
        "wald_low_pp": overall_rd_raw["wald_low"] * 100.0,
        "wald_high_pp": overall_rd_raw["wald_high"] * 100.0,
        "hksj_low_pp": overall_rd_raw["hksj_low"] * 100.0,
        "hksj_high_pp": overall_rd_raw["hksj_high"] * 100.0,
        "pi_df20_low_pp": overall_rd_raw["pi_df_k2_low"] * 100.0,
        "pi_df20_high_pp": overall_rd_raw["pi_df_k2_high"] * 100.0,
        "pi_df21_low_pp": overall_rd_raw["pi_df_k1_low"] * 100.0,
        "pi_df21_high_pp": overall_rd_raw["pi_df_k1_high"] * 100.0,
    },
    "fixed_effect_and_mantel_haenszel_k22": {
        "FE_inverse_variance_RR": {
            "rr": overall_rr["est_fe_iv"],
            "ci_low": overall_rr["fe_iv_low"],
            "ci_high": overall_rr["fe_iv_high"],
            "z": overall_rr["z_fe_iv"],
            "p": overall_rr["p_fe_iv"],
        },
        "MH_exact_uncorrected_RR": mh_uncorrected,
        "MH_continuity_corrected_RR": mh_corrected,
    },
    "egger_regression_k22": egger_all,
    "subgroups": {
        "platform_30_units_equal_ctrl_split": platform_summary_equal,
        "platform_30_units_arm_proportional_ctrl_split_check": platform_summary_armprop,
        "platform_31_units_behrouz_flic_omp22_as_chimeric_sensitivity": platform_summary_behrouz_split,
        "challenge_route": {"groups": route_summary, "Q_between": route_q_between},
        "adjuvant": {"groups": adj_summary, "Q_between": adj_q_between},
    },
    "sensitivity_analyses": {
        "leave_one_publication_out_20": {
            "min_RR_publication": min(loo_pub_rows, key=lambda x: x["RR_DL"]),
            "max_RR_publication": max(loo_pub_rows, key=lambda x: x["RR_DL"]),
            "lowest_hksj_lower_bound": min(loo_pub_rows, key=lambda x: x["RR_hksj_low"]),
            "all_20_rows": loo_pub_rows,
        },
        "leave_one_experiment_out_22": {
            "min_RR_experiment": min(loo_exp_rows, key=lambda x: x["RR_DL"]),
            "max_RR_experiment": max(loo_exp_rows, key=lambda x: x["RR_DL"]),
            "all_22_rows": loo_exp_rows,
        },
        "omit_inferred_controls_k20": sens_omit_inferred_ctrl,
        "ip_only_k17": sens_ip_only,
        "omit_garg_2016_k21": sens_omit_garg,
        "pre_supplementary_k15": sens_pre_supp_15,
        "non_zero_controls_only_k4": sens_ec_pos_4,
        "zero_controls_only_k18": sens_ec_zero_18,
        "stated_followup_ge_7d_k19": sens_7d_plus_known,
    },
    "within_study_head_to_head": within_study_comparisons,
    "syrcle_summary_counts": syrcle_counts,
    "excluded_studies_counts": {
        "original_track1_excluded_8": sum(1 for x in EXCLUDED_STUDIES if x["exclusion_stage_group"] == "Original_Track1_Reextraction_8"),
        "supplementary_protein_peptide_survival_28": sum(1 for x in EXCLUDED_STUDIES if x["exclusion_stage_group"] == "Supplementary_Protein_Peptide_Survival_28"),
        "supplementary_platform_exclusions_25": sum(1 for x in EXCLUDED_STUDIES if x["exclusion_stage_group"] == "Supplementary_Platform_Exclusion_25"),
        "supplementary_no_verified_survival_16": sum(1 for x in EXCLUDED_STUDIES if x["exclusion_stage_group"] == "Supplementary_No_Verified_Survival_16"),
        "total_table_s2_rows": len(EXCLUDED_STUDIES),
    },
}

with (OUT_DIR / "track1_stats.json").open("w", encoding="utf-8") as f:
    json.dump(stats_payload, f, indent=2)

print("=== TRACK 1 EXACT VERIFICATION SUMMARY ===")
print(f"Experiments: k = {len(exp_rows)} across {len(pub_order)} publications | Arms = {len(arm_rows)} | Platform units = {len(platform_units_30)}")
print(f"Total mice: N = {stats_payload['dataset_totals']['total_mice']} "
      f"(Vacc: {stats_payload['dataset_totals']['vaccinated_survivors']}/{stats_payload['dataset_totals']['vaccinated_mice']} = "
      f"{stats_payload['dataset_totals']['vaccinated_survival_pct']:.2f}%; "
      f"Ctrl: {stats_payload['dataset_totals']['control_survivors']}/{stats_payload['dataset_totals']['control_mice']} = "
      f"{stats_payload['dataset_totals']['control_survival_pct']:.2f}%)")
print(f"var(lnRR) range in ec=0 (k=18): {min(var_ec0):.4f} to {max(var_ec0):.4f} | in ec>0 (k=4): {min(var_ec_pos):.4f} to {max(var_ec_pos):.4f}")
print(f"Overall RR_DL = {overall_rr['est_dl']:.4f} "
      f"[Wald 95% CI: {overall_rr['wald_low']:.4f}-{overall_rr['wald_high']:.4f}; "
      f"HKSJ 95% CI: {overall_rr['hksj_low']:.4f}-{overall_rr['hksj_high']:.4f}; "
      f"PI(df=20): {overall_rr['pi_df_k2_low']:.4f}-{overall_rr['pi_df_k2_high']:.4f}; "
      f"PI(df=21): {overall_rr['pi_df_k1_low']:.4f}-{overall_rr['pi_df_k1_high']:.4f}]")
print(f"Heterogeneity RR: Q = {overall_rr['Q']:.4f} (df={overall_rr['df']}, p={overall_rr['p_Q']:.4f}), "
      f"I2 = {overall_rr['I2']:.2f}%, tau2 = {overall_rr['tau2_dl']:.6f}, q_HKSJ = {overall_rr['q_hksj_raw']:.4f}")
print(f"FE-IV RR = {overall_rr['est_fe_iv']:.4f} [{overall_rr['fe_iv_low']:.4f}-{overall_rr['fe_iv_high']:.4f}]; "
      f"MH uncorrected RR = {mh_uncorrected['rr_mh']:.4f} [{mh_uncorrected['ci_low']:.4f}-{mh_uncorrected['ci_high']:.4f}]; "
      f"MH +0.5 RR = {mh_corrected['rr_mh']:.4f} [{mh_corrected['ci_low']:.4f}-{mh_corrected['ci_high']:.4f}]")
print(f"Overall RD (Haldane pt+var, RevMan/Draft numbers) = +{overall_rd_cc['est_dl']*100:.2f} pp "
      f"[Wald: {overall_rd_cc['wald_low']*100:.2f}-{overall_rd_cc['wald_high']*100:.2f}; "
      f"HKSJ: {overall_rd_cc['hksj_low']*100:.2f}-{overall_rd_cc['hksj_high']*100:.2f}; "
      f"PI(df=20): {overall_rd_cc['pi_df_k2_low']*100:.2f}-{overall_rd_cc['pi_df_k2_high']*100:.2f}] | "
      f"Q_RD = {overall_rd_cc['Q']:.4f} (p={overall_rd_cc['p_Q']:.6f}), I2_RD = {overall_rd_cc['I2']:.2f}%, tau2_RD = {overall_rd_cc['tau2_dl']:.6f}")
print(f"Overall RD (Raw pt + Haldane var, Draft prose formula) = +{overall_rd_raw['est_dl']*100:.2f} pp "
      f"[Wald: {overall_rd_raw['wald_low']*100:.2f}-{overall_rd_raw['wald_high']*100:.2f}; "
      f"HKSJ: {overall_rd_raw['hksj_low']*100:.2f}-{overall_rd_raw['hksj_high']*100:.2f}; "
      f"PI(df=20): {overall_rd_raw['pi_df_k2_low']*100:.2f}-{overall_rd_raw['pi_df_k2_high']*100:.2f}] | "
      f"Q_RD = {overall_rd_raw['Q']:.4f} (p={overall_rd_raw['p_Q']:.6f}), I2_RD = {overall_rd_raw['I2']:.2f}%, tau2_RD = {overall_rd_raw['tau2_dl']:.6f}")
print(f"Egger test: intercept = {egger_all['intercept']:.4f}, SE = {egger_all['se_intercept']:.4f}, "
      f"t = {egger_all['t']:.4f} (df={egger_all['df']}), p = {egger_all['p']:.4e}")

print("\n--- Challenge Route Subgroups ---")
for rt, rdata in route_summary.items():
    pr = rdata["RR"]
    h_str = f"{pr['hksj_low']:.2f}-{pr['hksj_high']:.2f}" if pr["hksj_low"] else "NE"
    print(f"  {rt}: k={rdata['k']}, pubs={rdata['publications_count']}, vacc={rdata['vacc_events']}/{rdata['vacc_n']} ({rdata['vacc_pct']:.1f}%), "
          f"ctrl={rdata['ctrl_events']}/{rdata['ctrl_n']} ({rdata['ctrl_pct']:.1f}%), RR={pr['est_dl']:.2f} [{pr['wald_low']:.2f}-{pr['wald_high']:.2f}], HKSJ={h_str}")
print("  Route Q_between:", route_q_between)

print("\n--- Adjuvant Subgroups ---")
for ad, adata in adj_summary.items():
    pr = adata["RR"]
    h_str = f"{pr['hksj_low']:.2f}-{pr['hksj_high']:.2f}" if pr["hksj_low"] else "NE"
    print(f"  {ad}: k={adata['k']}, pubs={adata['publications_count']}, vacc={adata['vacc_events']}/{adata['vacc_n']} ({adata['vacc_pct']:.1f}%), "
          f"ctrl={adata['ctrl_events']}/{adata['ctrl_n']} ({adata['ctrl_pct']:.1f}%), RR={pr['est_dl']:.2f} [{pr['wald_low']:.2f}-{pr['wald_high']:.2f}], HKSJ={h_str}")
print("  Adjuvant Q_between:", adj_q_between)

print("\n--- Sensitivity Analyses ---")
for s_key in ["omit_inferred_controls_k20", "ip_only_k17", "omit_garg_2016_k21", "pre_supplementary_k15", "non_zero_controls_only_k4", "zero_controls_only_k18", "stated_followup_ge_7d_k19"]:
    sd = stats_payload["sensitivity_analyses"][s_key]
    pr = sd["RR"]
    rd_cc = sd["RD_haldane_both"]
    rd_rw = sd["RD_raw_pt_haldane_var"]
    print(f"  {s_key} (k={sd['k']}, N={sd['total_mice']}): RR={pr['est_dl']:.2f} [Wald {pr['wald_low']:.2f}-{pr['wald_high']:.2f}; HKSJ {pr['hksj_low']:.2f}-{pr['hksj_high']:.2f}], "
          f"I2={pr['I2']:.1f}% | RD_cc=+{rd_cc['est_dl']*100:.1f} [HKSJ {rd_cc['hksj_low']*100:.1f}-{rd_cc['hksj_high']*100:.1f}] | "
          f"RD_raw=+{rd_rw['est_dl']*100:.1f} [HKSJ {rd_rw['hksj_low']*100:.1f}-{rd_rw['hksj_high']*100:.1f}]")

print("\n--- Leave-One-Publication-Out Extremes (20 pubs) ---")
loo_min = stats_payload["sensitivity_analyses"]["leave_one_publication_out_20"]["min_RR_publication"]
loo_max = stats_payload["sensitivity_analyses"]["leave_one_publication_out_20"]["max_RR_publication"]
loo_low = stats_payload["sensitivity_analyses"]["leave_one_publication_out_20"]["lowest_hksj_lower_bound"]
print(f"  Min RR: omit {loo_min['omitted_publication']} -> RR={loo_min['RR_DL']:.2f} [HKSJ {loo_min['RR_hksj_low']:.2f}-{loo_min['RR_hksj_high']:.2f}]")
print(f"  Max RR: omit {loo_max['omitted_publication']} -> RR={loo_max['RR_DL']:.2f} [HKSJ {loo_max['RR_hksj_low']:.2f}-{loo_max['RR_hksj_high']:.2f}]")
print(f"  Lowest HKSJ lower bound: omit {loo_low['omitted_publication']} -> {loo_low['RR_hksj_low']:.2f}")

print("\n--- Within-Study Head-to-Head Comparisons ---")
for w in within_study_comparisons:
    print(f"  {w['publication']}: multi={w['multi_pooled']} vs single={w['single_pooled']} | best_multi={w['best_multi']} vs best_single={w['best_single']} | Fisher_p_pooled={w['fisher_p_pooled']:.4f} | {w['direction']}")
print("\n--- SYRCLE & Ethics Counts ---")
print(json.dumps(syrcle_counts, indent=2))
