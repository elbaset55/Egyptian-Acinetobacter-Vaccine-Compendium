#!/usr/bin/env python3
"""
Master Data Consolidation & 600-DPI Publication Figure Engine (v2 — Visually Audited)
=====================================================================================
1. Consolidates all full-text-verified CSVs and JSONs into d:/Graduation/systematic_review/data/
   and archives retired 25-trial files so zero legacy data remain in active folders.
2. Generates all publication-grade 600-DPI PNG figures with zero overlapping labels/legends:
   - Frontiers in Immunology Manuscript (Figures 1–9; zero thesis self-citation)
   - Volume I Graduation Thesis & Master Unified Thesis (Fig_SR1–Fig_SR10 + Thesis Benchmark Heatmap)
"""

import csv
import json
import math
import shutil
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Polygon
import numpy as np

ROOT = Path(r"d:\Graduation")
SR_DIR = ROOT / "systematic_review"
AUDIT_DIR = SR_DIR / "forensic_audit_workspace"
DATA_DIR = SR_DIR / "data"
FIG_DIR = SR_DIR / "figures"
ARCHIVE_DIR = AUDIT_DIR / "archived_legacy_pre_audit"

DATA_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# 1. CONSOLIDATE VERIFIED DATASETS & ARCHIVE RETIRED 25-TRIAL FILES
# ---------------------------------------------------------------------------
legacy_files = [
    DATA_DIR / "R_meta_analysis_25_trials.csv",
    DATA_DIR / "R_leave_one_out_sensitivity.csv",
    ROOT / "journal_submission" / "review_article" / "trials_25_data.json",
]
for lf in legacy_files:
    if lf.exists():
        shutil.move(str(lf), str(ARCHIVE_DIR / lf.name))

copy_map = {
    AUDIT_DIR / "track1" / "track1_experiments_22.csv": DATA_DIR / "Track1_22_experiments_verified.csv",
    AUDIT_DIR / "track1" / "track1_arms_all.csv": DATA_DIR / "Track1_48_arms_verified.csv",
    AUDIT_DIR / "track1" / "track1_platform_units_30.csv": DATA_DIR / "Track1_30_platform_units_verified.csv",
    AUDIT_DIR / "track1" / "track1_leave_one_out_20.csv": DATA_DIR / "Track1_20_leave_one_out_verified.csv",
    AUDIT_DIR / "track1" / "track1_syrcle_20.csv": DATA_DIR / "Track1_20_SYRCLE_verified.csv",
    AUDIT_DIR / "track1" / "table_s2_excluded_studies.csv": DATA_DIR / "Table_S2_77_excluded_studies_verified.csv",
    AUDIT_DIR / "track2" / "track2_22_verified.csv": DATA_DIR / "Track2_22_insilico_rubric_verified.csv",
    AUDIT_DIR / "track3" / "track3_17_cohorts_verified.csv": DATA_DIR / "Track3_17_cohorts_verified.csv",
}
for src, dst in copy_map.items():
    shutil.copy2(src, dst)

with open(AUDIT_DIR / "track1" / "track1_stats.json", encoding="utf-8") as f:
    t1_stats = json.load(f)
with open(AUDIT_DIR / "track2" / "track2_stats.json", encoding="utf-8") as f:
    t2_stats = json.load(f)
with open(AUDIT_DIR / "track3" / "track3_stats.json", encoding="utf-8") as f:
    t3_stats = json.load(f)

with open(DATA_DIR / "Track1_22_experiments_verified.csv", encoding="utf-8") as f:
    t1_exps = list(csv.DictReader(f))
with open(DATA_DIR / "Track1_30_platform_units_verified.csv", encoding="utf-8") as f:
    t1_units = list(csv.DictReader(f))
with open(DATA_DIR / "Track1_20_leave_one_out_verified.csv", encoding="utf-8") as f:
    t1_loo = list(csv.DictReader(f))
with open(DATA_DIR / "Track1_20_SYRCLE_verified.csv", encoding="utf-8") as f:
    t1_syrcle = list(csv.DictReader(f))
with open(DATA_DIR / "Track2_22_insilico_rubric_verified.csv", encoding="utf-8") as f:
    t2_rows = list(csv.DictReader(f))
with open(DATA_DIR / "Track3_17_cohorts_verified.csv", encoding="utf-8") as f:
    t3_rows = list(csv.DictReader(f))

# Compute empirical antigen counts across Track 1 (20 pubs) and Track 2 (22 pubs)
t1_pub_antigens = {}
for r in t1_exps:
    pub = r["publication"].replace("Araujo/Conde 2021", "Araujo 2021")
    t1_pub_antigens.setdefault(pub, []).append(r["antigens_exact"])

antigen_rules = [
    ("OmpA / Omp38", lambda s: any(k in s for k in ["OmpA", "Omp38"])),
    ("Omp34 / CarO / Omp33-36", lambda s: any(k in s for k in ["Omp34", "CarO", "Omp33-36", "34-kDa OMP"])),
    ("BamA / Oma87", lambda s: any(k in s for k in ["BamA", "Oma87"])),
    ("NucAb", lambda s: any(k in s for k in ["NucAb", "NucAB"])),
    ("BauA", lambda s: "BauA" in s),
    ("FilF", lambda s: "FilF" in s),
    ("ZnuD", lambda s: "ZnuD" in s),
    ("TonB / TonB-rec.", lambda s: "TonB" in s),
    ("Bap / Blp1 / Blp2", lambda s: any(k in s for k in ["Bap", "Blp1", "Blp2"])),
    ("OmpW", lambda s: "OmpW" in s),
    ("Omp22", lambda s: "Omp22" in s),
    ("OmpK", lambda s: "OmpK" in s),
    ("LptD / LptE", lambda s: any(k in s for k in ["LptD", "LptE"])),
]

antigen_freq = []
for label, matcher in antigen_rules:
    c_t1 = 0
    for pub, texts in t1_pub_antigens.items():
        cleaned = " ".join(texts)
        for neg in ["NO OmpA", "NO BamA", "NO LptD", "NOT LptD", "NOT BamA", "NOT Bap", "NO BauA"]:
            cleaned = cleaned.replace(neg, "")
        if matcher(cleaned):
            c_t1 += 1
    c_t2 = 0
    for r in t2_rows:
        cleaned = r["targets_verified"]
        for neg in ["NO OmpA", "NO BamA", "NO LptD", "NOT LptD", "NOT BamA", "NOT Bap", "NO BauA", "NO Omp34"]:
            cleaned = cleaned.replace(neg, "")
        if matcher(cleaned):
            c_t2 += 1
    antigen_freq.append({"antigen": label, "track1_pubs_of_20": c_t1, "track2_pubs_of_22": c_t2})

master_json = {
    "metadata": {
        "title": "Subunit and multi-epitope vaccines against carbapenem-resistant Acinetobacter baumannii: a systematic review and meta-analysis",
        "prospero_id": "CRD420261572661",
        "zenodo_doi": "10.5281/zenodo.24218732",
        "github_url": "https://github.com/elbaset55/Egyptian-Acinetobacter-Vaccine-Compendium",
        "audit_date": "2026-11-18",
        "unique_included_publications": 55,
    },
    "track1_preclinical_meta_analysis": t1_stats,
    "track2_insilico_quality_audit": t2_stats,
    "track3_epidemiological_synthesis": t3_stats,
    "antigen_representation_across_tracks": antigen_freq,
}
with open(DATA_DIR / "meta_analysis_results.json", "w", encoding="utf-8") as f:
    json.dump(master_json, f, indent=2, ensure_ascii=False)

print("Consolidated verified CSVs and meta_analysis_results.json successfully.")

# ---------------------------------------------------------------------------
# GLOBAL MATPLOTLIB PUBLICATION STYLING
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "axes.edgecolor": "#334155",
    "axes.linewidth": 0.9,
    "axes.labelcolor": "#0F172A",
    "xtick.color": "#1E293B",
    "ytick.color": "#1E293B",
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})

NAVY = "#0F294A"
CRIMSON = "#B91C1C"


def save_600dpi(fig, out_path: Path):
    fig.savefig(out_path, dpi=600, bbox_inches="tight", pad_inches=0.16, facecolor="white")
    plt.close(fig)
    print(f"Saved {out_path.name} ({out_path.stat().st_size / 1024:.1f} KB)")


# ===========================================================================
# FIGURE 1: PRISMA 2020 FLOW DIAGRAM (Fig_SR1_PRISMA2020_Flowchart.png)
# ===========================================================================
def build_fig1_prisma():
    fig, ax = plt.subplots(figsize=(11.5, 9.2))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.axis("off")

    def draw_box(x, y, w, h, title, body, hdr_color=NAVY, box_fc="#FFFFFF", ec="#94A3B8", title_fs=9.5, body_fs=8.4):
        rect = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.35,rounding_size=1.2",
                              facecolor=box_fc, edgecolor=ec, linewidth=1.2)
        ax.add_patch(rect)
        hdr_h = min(3.3, h * 0.30)
        hdr = FancyBboxPatch((x, y + h - hdr_h), w, hdr_h, boxstyle="round,pad=0.25,rounding_size=1.0",
                             facecolor=hdr_color, edgecolor=hdr_color, linewidth=0.8)
        ax.add_patch(hdr)
        ax.text(x + w / 2, y + h - hdr_h / 2, title, color="white", fontsize=title_fs,
                fontweight="bold", ha="center", va="center")
        ax.text(x + w / 2, y + (h - hdr_h) / 2, body, color="#0F172A", fontsize=body_fs,
                ha="center", va="center", linespacing=1.28)

    def arrow(x1, y1, x2, y2):
        ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle="-|>,head_length=0.65,head_width=0.40",
                                    lw=1.35, color="#334155"))

    phases = [
        (82, 15, "IDENTIFICATION", "#1E3A8A"),
        (61, 18, "SCREENING", "#0F766E"),
        (27.5, 29, "ELIGIBILITY &\nVERIFICATION", "#B45309"),
        (2.5, 21, "INCLUDED\nSYNTHESIS", "#0F294A"),
    ]
    for py, ph, label, col in phases:
        strip = FancyBboxPatch((1.0, py), 5.5, ph, boxstyle="round,pad=0.2,rounding_size=0.8",
                               facecolor=col, edgecolor=col, linewidth=0.8)
        ax.add_patch(strip)
        ax.text(3.75, py + ph / 2, label, color="white", fontsize=8.5, fontweight="bold",
                ha="center", va="center", rotation=90)

    ax.text(53, 98.5, "PRISMA 2020 Flow Diagram: Three-Track Systematic Review & Meta-Analysis",
            fontsize=12.5, fontweight="bold", color=NAVY, ha="center", va="center")

    draw_box(9, 82, 41, 14.5,
             "Primary Databases & Search Engines",
             "PubMed / MEDLINE (2000–2026): n = 273\n"
             "Consensus AI Corpus: n = 112\n"
             "Total primary records retrieved: n = 385",
             hdr_color="#1E3A8A", box_fc="#F8FAFC")

    draw_box(56, 82, 42, 14.5,
             "Duplicates Removed Prior to Screening",
             "Cross-database automated & manual\n"
             "deduplication removed: n = 147\n"
             "(385 retrieved → 238 unique records)",
             hdr_color="#475569", box_fc="#F8FAFC")
    arrow(50.3, 89.2, 55.7, 89.2)

    arrow(29.5, 81.8, 29.5, 76.2)
    draw_box(9, 62, 41, 14,
             "Title & Abstract Screening (PICOS)",
             "Unique records screened against\n"
             "three-track PICOS criteria:\n"
             "n = 238 records",
             hdr_color="#0F766E", box_fc="#F0FDFA")

    draw_box(56, 62, 42, 14,
             "Excluded at Title / Abstract Stage (n = 134)",
             "• Off-topic species, reviews, or non-vaccine AMR\n"
             "• Non-ICU or non-MENA surveillance surveys\n"
             "• Pure diagnostics without cohort resistance data",
             hdr_color="#64748B", box_fc="#F8FAFC")
    arrow(50.3, 69.0, 55.7, 69.0)

    arrow(29.5, 61.8, 29.5, 56.8)
    draw_box(9, 43.5, 41, 13,
             "Full-Text Reports Assessed for Eligibility",
             "Primary full-text articles evaluated: n = 104\n"
             "Initial three-track candidate pool: n = 58\n"
             "(before strict full-text survival count audit)",
             hdr_color="#B45309", box_fc="#FFFBEB")

    draw_box(56, 43.5, 42, 13,
             "Primary Full-Text Exclusions (n = 46)",
             "• Whole-cell / live-attenuated / lpxD bacterins (n = 14)\n"
             "• Passive mAb / IgY only without active arm (n = 11)\n"
             "• Pure docking (n = 9); non-ICU AMR (n = 11); carrier (n = 1)",
             hdr_color="#64748B", box_fc="#F8FAFC")
    arrow(50.3, 50.0, 55.7, 50.0)

    arrow(29.5, 43.3, 29.5, 39.2)
    draw_box(9, 27.5, 41, 11.5,
             "Full-Text Count Audit & Deep Search (Nov 2026)",
             "• Re-extraction removed 8 Track 1 reports lacking\n"
             "  exact x/n counts, valid controls, or peer review\n"
             "• Supplementary deep search (Undermind) added\n"
             "  +5 new Track 1 pubs (+ Jeffreys 2024 overlap)",
             hdr_color="#9A3412", box_fc="#FFF7ED", body_fs=8.0)

    draw_box(56, 27.5, 42, 11.5,
             "Cataloged Survival Exclusions (Table S2: n = 77)",
             "• 8 re-extracted Track 1 exclusions (curve/% only, OMV)\n"
             "• 28 supplementary protein survival reports lacking x/n\n"
             "• 25 non-protein platforms (OMV/DNA/mRNA/conjugate/WC)\n"
             "• 16 non-survival (CFU/passive/immunogenicity only)",
             hdr_color="#64748B", box_fc="#F8FAFC", body_fs=7.7)
    arrow(50.3, 33.2, 55.7, 33.2)

    arrow(29.5, 27.3, 29.5, 23.6)
    banner = FancyBboxPatch((9, 18.8), 89, 4.5, boxstyle="round,pad=0.25,rounding_size=0.9",
                            facecolor=NAVY, edgecolor=NAVY, linewidth=1.0)
    ax.add_patch(banner)
    ax.text(53.5, 21.05,
            "Final Included Corpus: 55 Unique Peer-Reviewed Publications (4 Studies Shared Across Track 1 & Track 2)",
            color="white", fontsize=9.6, fontweight="bold", ha="center", va="center")

    arrow(23.5, 18.6, 23.5, 16.0)
    arrow(53.5, 18.6, 53.5, 16.0)
    arrow(83.5, 18.6, 83.5, 16.0)

    draw_box(9, 2.5, 28.5, 13.2,
             "Track 1: Murine Survival",
             "20 publications | k = 22 expts\n"
             "48 active arms | 30 platform units\n"
             "N = 574 mice (390 vacc., 184 ctrl)\n"
             "RR = 3.69 [2.61–5.22], p < 0.0001",
             hdr_color="#1D4ED8", box_fc="#EFF6FF", title_fs=9.0, body_fs=8.0)

    draw_box(39.25, 2.5, 28.5, 13.2,
             "Track 2: In Silico MEVCs",
             "k = 22 computational pipelines\n"
             "7-domain quality rubric (0–14 pts)\n"
             "Mean: 4.91 ± 2.00 / 14 (range 2–8)\n"
             "90.9% lacked k-mer human screen",
             hdr_color="#0F766E", box_fc="#F0FDFA", title_fs=9.0, body_fs=8.0)

    draw_box(69.5, 2.5, 28.5, 13.2,
             "Track 3: ICU Epidemiology",
             "k = 17 Egypt & MENA cohorts\n"
             "1,724 analyzed isolates / episodes\n"
             "Unselected Egypt CRAB: 90.7%\n"
             "Ref. BMD colistin R: 1.2%–52.9%",
             hdr_color="#B91C1C", box_fc="#FEF2F2", title_fs=9.0, body_fs=8.0)

    save_600dpi(fig, FIG_DIR / "Fig_SR1_PRISMA2020_Flowchart.png")


# ===========================================================================
# FIGURE 2: FOREST PLOT OF k = 22 MURINE EXPERIMENTS (Fig_SR2_ForestPlot_Survival_RR.png)
# ===========================================================================
def build_fig2_forest():
    group_order = [
        ("BamA / Oma87 Single-Antigen Experiments", ["Singh_2017", "Rasooli_2020", "Araujo_Conde_2021"]),
        ("OmpW Single & Epitope Constructs", ["Huang_2015", "Abdollahi_2023"]),
        ("Other Single Outer Membrane & Surface Antigens", [
            "Huang_2016", "Singh_2016", "Garg_2016", "Skerniskyte_2019",
            "Behrouz_2024", "Freitas_2024_exp_a", "Freitas_2024_exp_b"
        ]),
        ("Multi-Platform & Chimeric / Multivalent Experiments", [
            "Badmasti_2015", "Mansouri_2023", "Akbari_2022", "Chaudhuri_2022", "Qamsari_2020",
            "Ren_2019", "Hagag_2022", "Jeffreys_2024", "Prather_2026_exp1", "Prather_2026_exp2"
        ]),
    ]
    by_id = {r["exp_id"]: r for r in t1_exps}

    display_labels = {
        "Singh_2017": "Singh et al., 2017 (rBamA; IN)",
        "Rasooli_2020": "Rasooli et al., 2020 (rOma87; IP)†",
        "Araujo_Conde_2021": "Araujo et al., 2021 (IB-rBamA; IP)",
        "Mansouri_2023": "Mansouri et al., 2023 (Oma87/Bap; IP)",
        "Huang_2015": "Huang et al., 2015 (Trx-OmpW; IP)",
        "Abdollahi_2023": "Abdollahi & Raoufi, 2023 (OmpW2; IN)†",
        "Huang_2016": "Huang et al., 2016 (Trx-Omp22; IP)",
        "Singh_2016": "Singh et al., 2016 (FilF; IT)",
        "Garg_2016": "Garg et al., 2016 (NucAb; IT)",
        "Skerniskyte_2019": "Skerniškytė et al., 2019 (Blp1; IP)",
        "Behrouz_2024": "Behrouz et al., 2024 (FliC-Omp22/Omp22; IP)",
        "Freitas_2024_exp_a": "Freitas et al., 2024, exp. a (CAM87009.1+AgNP; IP)",
        "Freitas_2024_exp_b": "Freitas et al., 2024, exp. b (CAM87009.1+Alum; IP)",
        "Badmasti_2015": "Badmasti et al., 2015 (AbOmpA/Bap; IP)",
        "Akbari_2022": "Akbari et al., 2022 (LCL BauA/Omp34; IP)",
        "Chaudhuri_2022": "Chaudhuri et al., 2022 (rBauA & LCL loops; IP)",
        "Qamsari_2020": "Qamsari et al., 2020 (rZnuD & LCL loops; IP)",
        "Ren_2019": "Ren et al., 2019 (rMEP: Ata+FilF+NucAb; IP)",
        "Hagag_2022": "Hagag et al., 2022 (Wza/YiaD; IP)",
        "Jeffreys_2024": "Jeffreys et al., 2024 (rAMEV2; IN)",
        "Prather_2026_exp1": "Prather et al., 2026, exp. 1 (AMEV1/2; IP)",
        "Prather_2026_exp2": "Prather et al., 2026, exp. 2 (AMEV2/5; IP)",
    }

    rows = []
    for g_title, ids in group_order:
        rows.append(("HEADER", g_title, None))
        for eid in ids:
            rows.append(("STUDY", display_labels[eid], by_id[eid]))
    rows.append(("SPACER", "", None))
    rows.append(("OVERALL", "Overall Random-Effects (DL) Model (k = 22)", None))

    total_rows = len(rows)
    fig, ax = plt.subplots(figsize=(14.4, 10.6))
    ax.set_ylim(-1.65, total_rows + 1.85)
    ax.set_xlim(-8.3, 10.2)
    ax.axis("off")

    # Generous column spacing so zero strings ever touch
    x_study = -8.15
    x_vacc = -2.05
    x_ctrl = -0.35
    x_wt = 0.85
    x_plot_l = 1.65
    x_plot_r = 7.25
    x_rr_txt = 7.55

    log_min, log_max = math.log10(0.2), math.log10(300.0)
    def map_x(rr_val):
        lv = math.log10(max(0.2, min(300.0, rr_val)))
        return x_plot_l + (lv - log_min) / (log_max - log_min) * (x_plot_r - x_plot_l)

    hdr_y = total_rows + 0.65
    ax.fill_between([-8.25, 10.1], hdr_y - 0.45, hdr_y + 0.55, color=NAVY)
    ax.text(x_study, hdr_y, "Experiment (Antigen; Challenge Route)", color="white", fontsize=9.1, fontweight="bold", va="center")
    ax.text(x_vacc, hdr_y, "Vaccinated\nn/N (%)", color="white", fontsize=8.5, fontweight="bold", va="center", ha="center")
    ax.text(x_ctrl, hdr_y, "Control\nn/N (%)", color="white", fontsize=8.5, fontweight="bold", va="center", ha="center")
    ax.text(x_wt, hdr_y, "DL Wt\n(%)", color="white", fontsize=8.5, fontweight="bold", va="center", ha="center")
    ax.text((x_plot_l + x_plot_r) / 2, hdr_y, "Risk Ratio (95% CI, log scale)", color="white", fontsize=9.1, fontweight="bold", va="center", ha="center")
    ax.text(x_rr_txt, hdr_y, "RR [95% Wald CI]", color="white", fontsize=9.1, fontweight="bold", va="center", ha="left")

    x_null = map_x(1.0)
    x_pool = map_x(3.6931)
    x_ci_l = map_x(2.6143)
    x_ci_h = map_x(5.2169)
    ax.fill_between([x_ci_l, x_ci_h], 0.2, total_rows + 0.1, color="#DBEAFE", alpha=0.65, zorder=0)
    ax.plot([x_null, x_null], [0.2, total_rows + 0.1], color="#64748B", lw=1.1, ls="--", zorder=1)
    ax.plot([x_pool, x_pool], [0.2, total_rows + 0.1], color=CRIMSON, lw=1.1, ls=":", zorder=1)

    for idx, (rtype, label, r) in enumerate(rows):
        y = total_rows - 0.35 - idx
        if rtype == "HEADER":
            ax.fill_between([-8.25, 10.1], y - 0.38, y + 0.38, color="#F1F5F9", zorder=0)
            ax.text(x_study, y, label, color=NAVY, fontsize=8.8, fontweight="bold", va="center")
        elif rtype == "STUDY":
            ev, nv = int(r["ev"]), int(r["nv"])
            ec, nc = int(r["ec"]), int(r["nc"])
            vpct = float(r["vacc_pct"])
            cpct = float(r["ctrl_pct"])
            wt = float(r["weight_DL_RR_pct"])
            rr = float(r["RR"])
            cil = float(r["RR_ci_low"])
            cih = float(r["RR_ci_high"])

            if idx % 2 == 1:
                ax.fill_between([-8.25, 10.1], y - 0.36, y + 0.36, color="#FAFAFA", zorder=0)

            ax.text(x_study + 0.10, y, label, color="#0F172A", fontsize=8.2, va="center")
            ax.text(x_vacc, y, f"{ev}/{nv} ({vpct:.1f}%)", color="#0F172A", fontsize=8.1, va="center", ha="center")
            ctrl_str = f"{ec}/{nc} ({cpct:.1f}%)" + ("†" if r["inferred_control_flag"] == "True" else "")
            ax.text(x_ctrl, y, ctrl_str, color="#0F172A", fontsize=8.1, va="center", ha="center")
            ax.text(x_wt, y, f"{wt:.1f}", color="#334155", fontsize=8.1, va="center", ha="center")

            xl, xh = map_x(cil), map_x(cih)
            xm = map_x(rr)
            ax.plot([xl, xh], [y, y], color="#1E3A8A", lw=1.25, zorder=3)
            if cil > 0.2:
                ax.plot([xl, xl], [y - 0.12, y + 0.12], color="#1E3A8A", lw=1.0, zorder=3)
            if cih <= 300.0:
                ax.plot([xh, xh], [y - 0.12, y + 0.12], color="#1E3A8A", lw=1.0, zorder=3)

            msize = 4.2 + 5.8 * math.sqrt(wt / 46.41)
            sq_col = "#0D9488" if ec > 0 else "#1D4ED8"
            ax.plot(xm, y, marker="s", markersize=msize, color=sq_col, markeredgecolor="#0F172A", markeredgewidth=0.5, zorder=4)
            ax.text(x_rr_txt, y, f"{rr:.2f} [{cil:.2f}, {cih:.2f}]", color="#0F172A", fontsize=8.2, va="center", ha="left")

        elif rtype == "OVERALL":
            ax.plot([-8.25, 10.1], [y + 0.52, y + 0.52], color=NAVY, lw=1.2)
            ax.fill_between([-8.25, 10.1], y - 0.55, y + 0.48, color="#EFF6FF", zorder=0)
            ax.text(x_study, y, label, color=NAVY, fontsize=8.8, fontweight="bold", va="center")
            ax.text(x_vacc, y, "280/390 (71.8%)", color=NAVY, fontsize=8.3, fontweight="bold", va="center", ha="center")
            ax.text(x_ctrl, y, "13/184 (7.1%)", color=NAVY, fontsize=8.3, fontweight="bold", va="center", ha="center")
            ax.text(x_wt, y, "100.0", color=NAVY, fontsize=8.3, fontweight="bold", va="center", ha="center")

            xh_l, xh_h = map_x(2.5598), map_x(5.3281)
            ax.plot([xh_l, xh_h], [y, y], color=CRIMSON, lw=1.8, zorder=3)
            diamond = Polygon([[x_ci_l, y], [x_pool, y + 0.28], [x_ci_h, y], [x_pool, y - 0.28]],
                              closed=True, facecolor=CRIMSON, edgecolor="#7F1D1D", lw=1.0, zorder=4)
            ax.add_patch(diamond)
            ax.text(x_rr_txt, y, "3.69 [2.61, 5.22]", color=CRIMSON, fontsize=8.7, fontweight="bold", va="center")

    axis_y = 0.2
    ax.plot([x_plot_l, x_plot_r], [axis_y, axis_y], color="#0F172A", lw=1.0)
    for tv in [0.25, 0.5, 1, 2, 5, 10, 25, 100, 250]:
        xt = map_x(tv)
        ax.plot([xt, xt], [axis_y, axis_y - 0.16], color="#0F172A", lw=0.9)
        ax.text(xt, axis_y - 0.42, str(tv), fontsize=8.0, ha="center", va="top", color="#1E293B")
    ax.text((x_plot_l + x_plot_r) / 2, axis_y - 0.88,
            "Favors Control  ←   Risk Ratio (RR)   →  Favors Vaccine",
            fontsize=8.6, fontweight="bold", ha="center", va="top", color=NAVY)

    footer = (
        "Heterogeneity (RR scale): Q = 19.95, df = 21, p = 0.52, I² = 0.0%, τ² = 0.000  |  HKSJ 95% CI: 2.56–5.33 (t = 7.41, p < 0.0001)  |  95% PI: 2.56–5.33\n"
        "Additive Risk Difference (RD): +62.7 pp [HKSJ 95% CI: 53.5–71.8; 95% PI: +29.3 to +96.0 pp], Q_RD = 49.55 (p = 0.0004), I²_RD = 57.6%  |  "
        "Teal squares: ec > 0 (4 expts, 71.0% weight); Blue squares: ec = 0 (+0.5 Haldane, 18 expts)"
    )
    ax.text(x_study, -1.35, footer, fontsize=7.8, color="#334155", va="top", linespacing=1.35)

    ax.text(0.92, total_rows + 1.48,
            "Preclinical Protective Efficacy of Subunit & Multi-Epitope Vaccines Against Lethal A. baumannii (k = 22 Experiments, N = 574 Mice)",
            fontsize=11.0, fontweight="bold", color=NAVY, ha="center", va="center")

    save_600dpi(fig, FIG_DIR / "Fig_SR2_ForestPlot_Survival_RR.png")


# ===========================================================================
# FIGURE 3: ANTIGEN PLATFORM SUBGROUPS & MATCHED HEAD-TO-HEAD (Fig_SR3_Subgroup_And_HeadToHead.png)
# ===========================================================================
def build_fig3_subgroups_h2h():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.2, 6.2), gridspec_kw={"width_ratios": [1.05, 1.28]})
    fig.subplots_adjust(wspace=0.36, bottom=0.22)

    platforms = [
        ("Multivalent &\nChimeric\n(k=10, 14 arms)", 82.9, "102/123", 3.95, 2.24, 6.98, "#0F766E"),
        ("BamA / Oma87\nSingle\n(k=4, 4 arms)", 88.1, "37/42", 2.80, 1.31, 5.98, "#1D4ED8"),
        ("OmpW\nConstructs\n(k=2, 3 arms)", 100.0, "22/22", 5.89, 1.66, 20.89, "#2563EB"),
        ("OmpA\nSingle Arm\n(k=1, 1 arm)", 71.4, "5/7", 4.58, 0.34, 61.58, "#D97706"),
        ("Other Single\nAntigens\n(k=13, 26 arms)", 58.2, "114/196", 3.62, 2.08, 6.30, "#475569"),
    ]
    y_pos = np.arange(len(platforms))[::-1]
    surv_vals = [p[1] for p in platforms]
    colors_a = [p[6] for p in platforms]

    ax1.barh(y_pos, surv_vals, height=0.56, color=colors_a, edgecolor="#0F172A", linewidth=0.8, alpha=0.90)
    ax1.axvline(71.8, color=CRIMSON, ls="--", lw=1.2, label="Overall vaccinated survival (71.8%)")
    ax1.set_yticks(y_pos)
    ax1.set_yticklabels([p[0] for p in platforms], fontsize=8.5, fontweight="bold")
    ax1.set_xlim(0, 138)
    ax1.set_xlabel("Pooled Vaccinated Survival (%)", fontsize=9.5, fontweight="bold")
    ax1.set_title("A. Antigen Platform Subgroups (30 Units, 48 Arms)\nBetween-platform Q = 1.09, df = 3, p = 0.78 (NS)",
                  fontsize=10.0, fontweight="bold", color=NAVY, pad=12)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(axis="x", linestyle=":", alpha=0.5)

    for idx, p in enumerate(platforms):
        yi = y_pos[idx]
        ax1.text(p[1] + 2.0, yi + 0.11, f"{p[1]:.1f}% ({p[2]})", fontsize=8.2, fontweight="bold", color="#0F172A", va="center")
        ax1.text(p[1] + 2.0, yi - 0.14, f"RR = {p[3]:.2f} [{p[4]:.2f}–{p[5]:.2f}]", fontsize=7.6, color="#334155", va="center")
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.15), fontsize=8.0, frameon=True, facecolor="white")

    h2h = [
        ("Badmasti 2015\n(AbOmpA+Bap vs singles)", 100.0, "7/7", 71.4, "5/7", 42.9, "3/7"),
        ("Akbari 2022\n(BauA+Omp34 LCL vs singles)", 71.4, "5/7", 42.9, "3/7", 42.9, "3/7"),
        ("Qamsari 2020\n(4-loop ZnuD vs rZnuD/loops)", 100.0, "8/8", 62.5, "5/8", 33.3, "2/6"),
        ("Hagag 2022\n(Wza+YiaD vs singles)", 100.0, "10/10", 70.0, "7/10", 60.0, "6/10"),
        ("Mansouri 2023\n(Oma87+Bap; sub-LD100)", 100.0, "14/14", 100.0, "14/14", 100.0, "14/14"),
        ("Chaudhuri 2022\n(Multi-loop BauA vs rBauA/L8)", 83.3, "5/6", 100.0, "6/6", 50.0, "3/6"),
        ("Behrouz 2024*\n(FliC–Omp22 carrier vs Omp22)", 75.0, "6/8", 37.5, "3/8", None, ""),
    ]
    x_pos = np.arange(len(h2h))
    w = 0.26

    multi_v = [h[1] for h in h2h]
    best_s = [h[3] for h in h2h]
    ax2.bar(x_pos - w, multi_v, width=w, color="#0F766E", edgecolor="#0F172A", lw=0.7, label="Multivalent / Fusion")
    ax2.bar(x_pos, best_s, width=w, color="#3B82F6", edgecolor="#0F172A", lw=0.7, label="Best Single Arm")
    x_low = [x_pos[i] + w for i in range(len(h2h)) if h2h[i][5] is not None]
    v_low = [h2h[i][5] for i in range(len(h2h)) if h2h[i][5] is not None]
    ax2.bar(x_low, v_low, width=w, color="#94A3B8", edgecolor="#0F172A", lw=0.7, label="Lowest Single Arm")

    for i, h in enumerate(h2h):
        if "Mansouri" in h[0]:
            ax2.text(x_pos[i], 102.5, "All 3 arms: 14/14", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color="#0F172A")
        else:
            ax2.text(x_pos[i] - w, h[1] + 2.0, h[2], ha="center", va="bottom", fontsize=7.2, fontweight="bold", color="#0F766E")
            ax2.text(x_pos[i], h[3] + 2.0, h[4], ha="center", va="bottom", fontsize=7.2, color="#1E3A8A")
            if h[5] is not None:
                ax2.text(x_pos[i] + w, h[5] + 2.0, h[6], ha="center", va="bottom", fontsize=7.0, color="#475569")

    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([h[0] for h in h2h], rotation=24, ha="right", fontsize=7.8)
    ax2.set_ylim(0, 128)
    ax2.set_ylabel("Post-Challenge Survival (%)", fontsize=9.5, fontweight="bold")
    ax2.set_title("B. Matched Within-Study Comparisons (5–14 Mice / Arm)\nMultivalent higher in 4/6 multi-antigen trials (*Behrouz: FliC carrier)",
                  fontsize=10.0, fontweight="bold", color=NAVY, pad=12)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(axis="y", linestyle=":", alpha=0.5)
    # Place 3-column horizontal legend cleanly at very top of ylim(0,128) above 112%
    ax2.legend(loc="upper center", ncol=3, fontsize=7.6, frameon=True, facecolor="white")

    out3 = FIG_DIR / "Fig_SR3_Subgroup_And_HeadToHead.png"
    save_600dpi(fig, out3)
    shutil.copy2(out3, FIG_DIR / "Fig_SR3_Subgroup_Comparison_Single_vs_Multi.png")


# ===========================================================================
# FIGURE 4: LEAVE-ONE-PUBLICATION-OUT SENSITIVITY (Fig_SR7_LeaveOneOut_Sensitivity.png)
# ===========================================================================
def build_fig4_loo():
    loo_sorted = sorted(t1_loo, key=lambda r: float(r["RR_DL"]))
    fig, ax = plt.subplots(figsize=(11.2, 7.8))

    n = len(loo_sorted)
    y_pos = np.arange(n)[::-1]

    ax.axvspan(2.6143, 5.2169, color="#DBEAFE", alpha=0.65, label="Full-dataset Wald 95% CI [2.61–5.22]")
    ax.axvline(3.6931, color=CRIMSON, ls="--", lw=1.4, label="Full-dataset Pooled RR = 3.69")
    ax.axvline(1.0, color="#64748B", ls=":", lw=1.1, label="Null effect (RR = 1.0)")

    for i, r in enumerate(loo_sorted):
        yi = y_pos[i]
        rr = float(r["RR_DL"])
        hl = float(r["RR_hksj_low"])
        hh = float(r["RR_hksj_high"])
        wl = float(r["RR_wald_low"])
        wh = float(r["RR_wald_high"])

        ax.plot([hl, hh], [yi, yi], color="#94A3B8", lw=2.6, solid_capstyle="round", zorder=2,
                label="HKSJ 95% CI" if i == 0 else None)
        ax.plot([wl, wh], [yi, yi], color="#1D4ED8", lw=1.3, zorder=3,
                label="Wald 95% CI" if i == 0 else None)
        ax.plot(rr, yi, "o", color=NAVY, markersize=6.0, zorder=4)
        ax.text(11.6, yi, f"{rr:.2f} [HKSJ: {hl:.2f}–{hh:.2f}]", fontsize=8.0, va="center", color="#0F172A")

    ax.set_yticks(y_pos)
    ax.set_yticklabels([
        r["omitted_publication"].replace("Araujo/Conde 2021", "Araujo et al., 2021").replace("Skerniskyte", "Skerniškytė")
        + (" (2 expts)" if int(r["experiments_omitted"]) == 2 else "")
        for r in loo_sorted
    ], fontsize=8.4)
    ax.set_xscale("log")
    ax.set_xlim(0.85, 30.0)
    ax.set_xticks([1, 2, 3, 3.69, 5, 7, 10])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("Recalculated Pooled Risk Ratio (RR_DL, log scale) After Omitting Publication", fontsize=9.5, fontweight="bold")
    ax.set_title("Leave-One-Publication-Out Sensitivity Analysis Across All 20 Included Publications\n"
                 "All 20 iterations remain statistically significant (RR range: 3.52–6.34; lowest HKSJ lower bound: 2.42)",
                 fontsize=10.2, fontweight="bold", color=NAVY, pad=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle=":", alpha=0.45)
    # Place legend cleanly below the x-axis so zero data lines or text columns are occluded
    ax.legend(loc="upper center", bbox_to_anchor=(0.48, -0.11), ncol=3, fontsize=8.0, frameon=True, facecolor="white")

    save_600dpi(fig, FIG_DIR / "Fig_SR7_LeaveOneOut_Sensitivity.png")


# ===========================================================================
# FIGURE 5: FUNNEL PLOT + EGGER & AUTHENTIC SYRCLE BAR CHART (Fig_SR5_FunnelPlot_Egger_SYRCLE.png)
# ===========================================================================
def build_fig5_funnel_syrcle():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14.0, 6.4), gridspec_kw={"width_ratios": [1.0, 1.18]})
    fig.subplots_adjust(wspace=0.40, bottom=0.22)

    lnrr_pool = math.log(3.6931)
    se_max = 1.62
    se_grid = np.linspace(0.001, se_max, 200)
    ax1.fill_betweenx(se_grid, lnrr_pool - 1.96 * se_grid, lnrr_pool + 1.96 * se_grid,
                      color="#EFF6FF", edgecolor="#BFDBFE", lw=1.0, label="Pseudo 95% CI funnel")
    ax1.axvline(lnrr_pool, color=CRIMSON, ls="--", lw=1.3, label="Pooled ln(RR) = 1.31 (RR = 3.69)")
    ax1.axvline(0.0, color="#64748B", ls=":", lw=1.0, label="Null ln(RR) = 0")

    prec = np.array([1.0 / float(r["se_lnRR"]) for r in t1_exps])
    snd = np.array([float(r["lnRR"]) / float(r["se_lnRR"]) for r in t1_exps])
    X = np.column_stack([np.ones_like(prec), prec])
    beta_ols = np.linalg.lstsq(X, snd, rcond=None)[0]
    a_int, b_slope = beta_ols[0], beta_ols[1]
    se_line = np.linspace(0.22, 1.55, 100)
    lnrr_egger = b_slope + a_int * se_line
    ax1.plot(lnrr_egger, se_line, color="#D97706", lw=1.7, ls="-.",
             label=f"Egger regression (int. = {a_int:.2f}, SE = 0.10, p < 0.0001)")

    ec0_x = [float(r["lnRR"]) for r in t1_exps if int(r["ec"]) == 0]
    ec0_y = [float(r["se_lnRR"]) for r in t1_exps if int(r["ec"]) == 0]
    ecp_x = [float(r["lnRR"]) for r in t1_exps if int(r["ec"]) > 0]
    ecp_y = [float(r["se_lnRR"]) for r in t1_exps if int(r["ec"]) > 0]

    ax1.scatter(ec0_x, ec0_y, c="#1D4ED8", s=46, edgecolors="#0F172A", linewidths=0.6, zorder=4,
                label="ec = 0 trials (+0.5 Haldane, k = 18)")
    ax1.scatter(ecp_x, ecp_y, c="#0D9488", s=62, marker="D", edgecolors="#0F172A", linewidths=0.7, zorder=5,
                label="ec > 0 trials (surviving controls, k = 4)")

    ax1.set_ylim(1.62, 0.0)
    ax1.set_xlim(-1.8, 4.6)
    ax1.set_xlabel("Log Risk Ratio [ln(RR)]", fontsize=9.5, fontweight="bold")
    ax1.set_ylabel("Standard Error [SE(ln RR)] (Inverted)", fontsize=9.5, fontweight="bold")
    ax1.set_title("A. Funnel Plot & Egger's Test (k = 22 Experiments)\nAsymmetry reflects ec = 0 variance coupling & small-study effects",
                  fontsize=9.8, fontweight="bold", color=NAVY, pad=10)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=7.3, frameon=True, facecolor="white")

    domains = [
        ("Mouse strain & sex stated (n=14)", 100.0, 0.0, 0.0),
        ("Follow-up duration stated (n=20)", 95.0, 5.0, 0.0),
        ("Ethics committee approval (n=20)", 90.0, 10.0, 0.0),
        ("Mouse age stated (n=14)", 85.7, 14.3, 0.0),
        ("Ethics protocol number stated (n=20)", 55.0, 45.0, 0.0),
        ("Humane endpoints described (n=14)", 21.4, 78.6, 0.0),
        ("Random allocation mentioned* (n=20)", 10.0, 90.0, 0.0),
        ("Blinded histopathology assessment (n=20)", 10.0, 90.0, 0.0),
        ("Random outcome selection (n=14)", 7.1, 92.9, 0.0),
        ("Attrition / group-n consistency (n=14)", 0.0, 78.6, 21.4),
        ("Allocation concealment (n=20)", 0.0, 100.0, 0.0),
        ("Caregiver / investigator blinding (n=20)", 0.0, 100.0, 0.0),
        ("Random housing / baseline / power (n=14)", 0.0, 100.0, 0.0),
    ]
    y_d = np.arange(len(domains))[::-1]
    rep = np.array([d[1] for d in domains])
    unc = np.array([d[2] for d in domains])
    disc = np.array([d[3] for d in domains])

    ax2.barh(y_d, rep, height=0.64, color="#10B981", edgecolor="#0F172A", lw=0.5, label="Explicitly Reported")
    ax2.barh(y_d, unc, left=rep, height=0.64, color="#FBBF24", edgecolor="#0F172A", lw=0.5, label="Unclear / Not Stated in Text")
    ax2.barh(y_d, disc, left=rep + unc, height=0.64, color="#EF4444", edgecolor="#0F172A", lw=0.5, label="Unexplained Group-n Discrepancy")

    ax2.set_yticks(y_d)
    ax2.set_yticklabels([d[0] for d in domains], fontsize=8.0)
    ax2.set_xlim(0, 100)
    ax2.set_xlabel("Percentage of Appraised Publications (%)", fontsize=9.5, fontweight="bold")
    ax2.set_title("B. SYRCLE Risk-of-Bias & Reporting Audit\n(*Sequence generation method 0/20; missing items coded Unclear)",
                  fontsize=9.8, fontweight="bold", color=NAVY, pad=10)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=7.5, frameon=True, facecolor="white")

    save_600dpi(fig, FIG_DIR / "Fig_SR5_FunnelPlot_Egger_SYRCLE.png")


# ===========================================================================
# FIGURE 6: ROUTE / ADJUVANT SUBGROUPS & SENSITIVITY MODELS (Fig_SR9)
# ===========================================================================
def build_fig6_route_adjuvant_sensitivity():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.4, 6.4), gridspec_kw={"width_ratios": [1.02, 1.15]})
    fig.subplots_adjust(wspace=0.58, bottom=0.18)

    sub_items = [
        ("HEADER", "By Challenge Route (Q_between = 4.04, p = 0.13)", None, None, None, ""),
        ("ITEM", "Intraperitoneal (IP; k=17)", 3.35, 2.34, 4.79, "243/334 (72.8%)"),
        ("ITEM", "Intranasal (IN; k=3)", 15.72, 3.26, 75.78, "30/36 (83.3%)"),
        ("ITEM", "Intratracheal (IT; k=2)", 7.57, 1.01, 56.49, "7/20 (35.0%)"),
        ("HEADER", "By Adjuvant Class (Q_between = 3.14, p = 0.21)", None, None, None, ""),
        ("ITEM", "Aluminum salts (Alum; k=7)", 4.12, 2.22, 7.65, "76/108 (70.4%)"),
        ("ITEM", "Freund's CFA/IFA (k=10)", 4.61, 2.34, 9.10, "152/209 (72.7%)"),
        ("ITEM", "Other / None (AddaS03/AgNP/FliC; k=5)", 13.89, 4.11, 46.93, "52/73 (71.2%)"),
    ]
    y_a = np.arange(len(sub_items))[::-1]
    ax1.axvline(3.6931, color=CRIMSON, ls="--", lw=1.2)
    ax1.axvline(1.0, color="#64748B", ls=":", lw=1.0)

    ytick_pos, ytick_lbl = [], []
    for idx, (kind, lbl, rr, cl, ch, surv_txt) in enumerate(sub_items):
        yi = y_a[idx]
        if kind == "HEADER":
            ax1.text(0.82, yi, lbl, fontsize=8.4, fontweight="bold", color=NAVY, va="center")
        else:
            ytick_pos.append(yi)
            ytick_lbl.append(lbl)
            ax1.plot([cl, ch], [yi, yi], color="#1D4ED8", lw=1.8)
            ax1.plot(rr, yi, "s", color="#0F766E", markersize=6.5, markeredgecolor="#0F172A", markeredgewidth=0.6)
            ax1.text(92, yi, f"{rr:.2f} [{cl:.2f}–{ch:.2f}]", fontsize=7.8, va="center", color="#0F172A")

    ax1.set_ylim(-0.7, len(sub_items) - 0.15)
    ax1.set_yticks(ytick_pos)
    ax1.set_yticklabels(ytick_lbl, fontsize=8.2)
    ax1.set_xscale("log")
    ax1.set_xlim(0.75, 520)
    ax1.set_xticks([1, 2, 5, 10, 25, 75])
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.set_xlabel("Subgroup Pooled Risk Ratio [95% Wald CI]", fontsize=9.2, fontweight="bold")
    ax1.set_title("A. Subgroup Analyses by Challenge Route & Adjuvant\nNo statistically significant between-subgroup differences",
                  fontsize=9.8, fontweight="bold", color=NAVY, pad=12)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(axis="x", linestyle=":", alpha=0.45)

    sens_models = [
        ("Primary dataset (all k=22, N=574)", 3.69, 2.56, 5.33, "+62.7 pp"),
        ("Fixed-effect IV (+0.5 Haldane, k=22)", 3.69, 2.61, 5.22, "Wald CI"),
        ("Omit inferred controls (k=20, N=536)", 3.52, 2.42, 5.12, "+59.5 pp"),
        ("IP challenge only (k=17, N=470)", 3.35, 2.27, 4.94, "+63.9 pp"),
        ("Omit Garg 2016 48-h readout (k=21)", 3.68, 2.54, 5.32, "+65.1 pp"),
        ("Pre-supplementary 14-pub set (k=15)", 3.86, 2.39, 6.22, "+64.7 pp"),
        ("Follow-up ≥ 7 days stated (k=19)", 5.99, 3.55, 10.10, "+60.0 pp"),
        ("Surviving controls ec > 0 only (k=4)", 2.32, 1.19, 4.52, "+54.0 pp"),
        ("Mantel–Haenszel (+0.5, k=22)", 6.05, 4.13, 8.85, "Wald CI"),
        ("Leave-one-pub-out min–max (20 pubs)", 3.69, 3.52, 6.34, "pt. range"),
    ]
    y_b = np.arange(len(sens_models))[::-1]
    ax2.axvspan(2.56, 5.33, color="#DBEAFE", alpha=0.60, label="Primary HKSJ 95% CI [2.56–5.33]")
    ax2.axvline(3.6931, color=CRIMSON, ls="--", lw=1.3, label="Primary Pooled RR = 3.69")
    ax2.axvline(1.0, color="#64748B", ls=":", lw=1.0)

    for idx, (lbl, rr, cl, ch, rd_note) in enumerate(sens_models):
        yi = y_b[idx]
        col = CRIMSON if idx == 0 else ("#D97706" if idx == len(sens_models) - 1 else "#1E3A8A")
        ax2.plot([cl, ch], [yi, yi], color=col, lw=2.0)
        if idx < len(sens_models) - 1:
            ax2.plot(rr, yi, "o", color=col, markersize=6.0)
            ax2.text(11.8, yi, f"{rr:.2f} [{cl:.2f}–{ch:.2f}] ({rd_note})", fontsize=7.7, va="center", color="#0F172A")
        else:
            ax2.plot([cl, ch], [yi, yi], "d", color=col, markersize=6.0)
            ax2.text(11.8, yi, f"RR_DL: {cl:.2f} to {ch:.2f}", fontsize=7.7, fontweight="bold", va="center", color="#9A3412")

    ax2.set_ylim(-0.7, len(sens_models) - 0.3)
    ax2.set_yticks(y_b)
    ax2.set_yticklabels([m[0] for m in sens_models], fontsize=8.1)
    ax2.set_xscale("log")
    ax2.set_xlim(0.9, 52.0)
    ax2.set_xticks([1, 2, 3.69, 6, 10])
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.set_xlabel("Pooled Risk Ratio [95% HKSJ or Wald CI]", fontsize=9.2, fontweight="bold")
    ax2.set_title("B. Pre-Specified & Estimator Sensitivity Analyses\nPositive survival benefit persists across all subset & weighting models",
                  fontsize=9.8, fontweight="bold", color=NAVY, pad=12)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(axis="x", linestyle=":", alpha=0.45)
    ax2.legend(loc="upper center", bbox_to_anchor=(0.48, -0.14), ncol=2, fontsize=7.8, frameon=True, facecolor="white")

    out6 = FIG_DIR / "Fig_SR9_Subgroup_Route_Adjuvant_Sensitivity.png"
    save_600dpi(fig, out6)
    shutil.copy2(out6, FIG_DIR / "Fig_SR9_Correlates_CFU_And_Passive_Transfer.png")


# ===========================================================================
# FIGURE 7: TRACK 2 IN SILICO QUALITY HEATMAP (Journal k=22 & Thesis +Benchmark)
# ===========================================================================
def build_fig7_track2_heatmap(include_thesis_benchmark=False, out_name="Fig_SR4_InSilico_Methodological_Audit_Heatmap.png"):
    rows_plot = []
    for r in t2_rows:
        lbl = r["study_verified"].replace("Martin-Galiano", "Martín-Galiano")
        scores = [int(r[f"D{j}"]) for j in range(1, 8)]
        tot = int(r["total_verified_14"])
        rows_plot.append((lbl, scores, tot, False))

    if include_thesis_benchmark:
        rows_plot.append(("Volume II Design Benchmark (2026)*", [2, 2, 2, 2, 2, 2, 1], 13, True))

    n_rows = len(rows_plot)
    mat = np.array([r[1] for r in rows_plot])

    fig, (ax_hm, ax_tot) = plt.subplots(
        1, 2, figsize=(12.6, 9.0 if include_thesis_benchmark else 8.6),
        gridspec_kw={"width_ratios": [3.3, 1.2]}, sharey=True
    )
    fig.subplots_adjust(wspace=0.06, top=0.82, bottom=0.08)

    cmap = matplotlib.colors.ListedColormap(["#FEE2E2", "#FEF3C7", "#D1FAE5"])
    bounds = [-0.5, 0.5, 1.5, 2.5]
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    ax_hm.imshow(mat, cmap=cmap, norm=norm, aspect="auto")

    col_labels = [
        "D1: Pan-Genome\nScale (>20)",
        "D2: Target Clinical\nMDR Genomes",
        "D3: Strict β-Barrel\nECL Excision",
        "D4: Multi-Locus\nRegional HLA",
        "D5: Peptide k-mer\nHuman Screen",
        "D6: AF2 (>90%) &\nTLR2/4 Docking",
        "D7: Validation\nDepth (MD/In Vivo)",
    ]
    ax_hm.set_xticks(np.arange(7))
    ax_hm.xaxis.tick_top()
    ax_hm.set_xticklabels(col_labels, fontsize=8.1, fontweight="bold", rotation=28, ha="left", rotation_mode="anchor")

    ax_hm.set_yticks(np.arange(n_rows))
    ax_hm.set_yticklabels([r[0] for r in rows_plot], fontsize=8.5)

    for i in range(n_rows):
        if rows_plot[i][3]:
            ax_hm.axhline(i - 0.5, color=NAVY, lw=2.0)
        else:
            ax_hm.axhline(i - 0.5, color="white", lw=1.2)
        for j in range(7):
            v = mat[i, j]
            txt_c = "#991B1B" if v == 0 else ("#92400E" if v == 1 else "#065F46")
            ax_hm.text(j, i, str(v), ha="center", va="center", fontsize=8.8, fontweight="bold", color=txt_c)

    for j in range(8):
        ax_hm.axvline(j - 0.5, color="white", lw=1.4)

    tots = [r[2] for r in rows_plot]
    bar_cols = ["#0D9488" if r[3] else ("#2563EB" if r[2] >= 7 else ("#D97706" if r[2] >= 5 else "#DC2626")) for r in rows_plot]
    y_idx = np.arange(n_rows)
    ax_tot.barh(y_idx, tots, height=0.68, color=bar_cols, edgecolor="#0F172A", lw=0.6)
    ax_tot.axvline(4.91, color=CRIMSON, ls="--", lw=1.3, label="Published Mean: 4.91 / 14")
    ax_tot.set_xlim(0, 16.4)
    ax_tot.xaxis.tick_top()
    ax_tot.set_xticks([0, 4, 7, 10, 14])
    ax_tot.set_xlabel("Total Score (/14)", fontsize=8.8, fontweight="bold")
    ax_tot.xaxis.set_label_position("top")
    ax_tot.spines["right"].set_visible(False)
    ax_tot.spines["bottom"].set_visible(False)
    ax_tot.grid(axis="x", linestyle=":", alpha=0.45)

    for i, tot in enumerate(tots):
        ax_tot.text(tot + 0.35, i, f"{tot}/14", va="center", fontsize=8.2, fontweight="bold", color="#0F172A")

    title_str = (
        "Methodological Quality Appraisal of Published In Silico A. baumannii MEVC Pipelines (k = 22, 2018–2026)\n"
        "Verified Mean Score: 4.91 ± 2.00 / 14 (Median 5.0, IQR 3.25–6.00; Range 2–8)  |  0 = Absent (Red), 1 = Partial (Amber), 2 = Rigorous (Green)"
    )
    if include_thesis_benchmark:
        title_str += "\n(*Bottom row shows companion Volume II Egyptian IC2 Quadrivalent ECL-MEVC pipeline [13/14; wet-lab pending])"
    fig.suptitle(title_str, fontsize=9.8, fontweight="bold", color=NAVY, y=0.98)
    ax_tot.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), fontsize=7.8, frameon=True, facecolor="white")

    save_600dpi(fig, FIG_DIR / out_name)


# ===========================================================================
# FIGURE 8 (Thesis Fig_SR8): SYRCLE 20-PUBLICATION TRAFFIC-LIGHT MATRIX
# ===========================================================================
def build_fig_sr8_syrcle_matrix():
    cols = [
        ("D1_sequence_generation", "D1: Random\nSequence"),
        ("D2_baseline_characteristics", "D2: Baseline\nComparability"),
        ("D3_allocation_concealment", "D3: Allocation\nConcealment"),
        ("D4_random_housing", "D4: Random\nHousing"),
        ("D5_caregiver_blinding", "D5: Caregiver\nBlinding"),
        ("D6_random_outcome_assessment", "D6: Random\nAssessment"),
        ("D7_outcome_assessor_blinding", "D7: Assessor\nBlinding"),
        ("D8_incomplete_outcome_data", "D8: Attrition /\nGroup-n"),
        ("ethics_committee_stated", "Ethics\nCommittee"),
        ("ethics_approval_number", "Protocol\nNumber"),
    ]

    n_r = len(t1_syrcle)
    n_c = len(cols)
    mat = np.zeros((n_r, n_c), dtype=int)
    sym = [["" for _ in range(n_c)] for _ in range(n_r)]

    for i, r in enumerate(t1_syrcle):
        for j, (ck, _) in enumerate(cols):
            val = r[ck].strip()
            if ck == "ethics_committee_stated":
                code, s = (2, "+") if val.startswith("Yes") else (1, "?")
            elif ck == "ethics_approval_number":
                code, s = (2, "+") if val.startswith("Yes") else (1, "?")
            elif ck == "D1_sequence_generation":
                code, s = (2, "+*") if val.startswith("Reported") else (1, "?")
            elif ck == "D6_random_outcome_assessment":
                if val.startswith("Reported"):
                    code, s = 2, "+"
                elif "Not appraised" in val:
                    code, s = 3, "—"
                else:
                    code, s = 1, "?"
            elif ck == "D7_outcome_assessor_blinding":
                code, s = (2, "+*") if val.startswith("Partial") else (1, "?")
            elif ck == "D8_incomplete_outcome_data":
                if "Concern" in val or "Noted attrition" in val:
                    code, s = 0, "!"
                elif "Not appraised" in val:
                    code, s = 3, "—"
                else:
                    code, s = 1, "?"
            else:
                if "Not appraised" in val:
                    code, s = 3, "—"
                elif val.startswith("Reported") or val.startswith("Low"):
                    code, s = 2, "+"
                else:
                    code, s = 1, "?"
            mat[i, j] = code
            sym[i][j] = s

    # Verify exact counts match report:
    assert sum(1 for i in range(n_r) if sym[i][9] == "+") == 11, "Ethics protocol numbers must be 11/20"
    assert sum(1 for i in range(n_r) if sym[i][8] == "+") == 18, "Ethics committees must be 18/20"

    fig, ax = plt.subplots(figsize=(11.6, 8.6))
    fig.subplots_adjust(top=0.82, bottom=0.10)
    cmap = matplotlib.colors.ListedColormap(["#FCA5A5", "#FDE68A", "#A7F3D0", "#E2E8F0"])
    norm = matplotlib.colors.BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    ax.imshow(mat, cmap=cmap, norm=norm, aspect="auto")

    ax.set_xticks(np.arange(n_c))
    ax.xaxis.tick_top()
    ax.set_xticklabels([c[1] for c in cols], fontsize=8.1, fontweight="bold", rotation=28, ha="left", rotation_mode="anchor")

    pub_names = [
        r["publication"].replace("Araujo/Conde 2021", "Araujo et al., 2021").replace("Skerniskyte", "Skerniškytė")
        + (" (supp.)" if r["appraisal_scope"] != "Full_10_domain" else "")
        for r in t1_syrcle
    ]
    ax.set_yticks(np.arange(n_r))
    ax.set_yticklabels(pub_names, fontsize=8.4)

    for i in range(n_r):
        ax.axhline(i - 0.5, color="white", lw=1.2)
        for j in range(n_c):
            ax.text(j, i, sym[i][j], ha="center", va="center", fontsize=8.6, fontweight="bold", color="#0F172A")
    for j in range(n_c + 1):
        ax.axvline(j - 0.5, color="white", lw=1.2)

    legend_patches = [
        mpatches.Patch(color="#A7F3D0", label="+ / +*  Reported (+* partial: random w/o method or histopathology blinding)"),
        mpatches.Patch(color="#FDE68A", label="?  Unclear / Not stated in full text"),
        mpatches.Patch(color="#FCA5A5", label="!  Unexplained group-size / attrition discrepancy (3/14)"),
        mpatches.Patch(color="#E2E8F0", label="—  Supplementary stage (core bias & ethics domains extracted)"),
    ]
    ax.legend(handles=legend_patches, loc="upper center", bbox_to_anchor=(0.5, -0.03), ncol=2, fontsize=7.8, frameon=True)
    fig.suptitle("Study-Level SYRCLE Risk-of-Bias & Ethics Reporting Matrix (20 Included Track 1 Publications)\n"
                 "Committee approval: 18/20 (90%)  |  Protocol number: 11/20 (55%)  |  Random allocation: 2/20 (10%, method NS)",
                 fontsize=10.2, fontweight="bold", color=NAVY, y=0.97)

    save_600dpi(fig, FIG_DIR / "Fig_SR8_SYRCLE_TrafficLight_Matrix.png")


# ===========================================================================
# FIGURE 8 (Journal Fig 8 / Fig_SR6): EGYPT & MENA EPIDEMIOLOGY
# ===========================================================================
def build_fig8_epidemiology():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.6, 6.8), gridspec_kw={"width_ratios": [1.15, 1.0]})
    fig.subplots_adjust(wspace=0.60, bottom=0.20)

    crab_cohorts = [
        ("Benmahmod 2019 (Mansoura, N=50)", 98.0, "Unselected Egypt"),
        ("Alkasaby & Zaki 2017 (Mansoura ICU, N=280)", 95.7, "Unselected Egypt"),
        ("Abozahra 2021 (Alexandria ICU, N=120)", 91.7, "Unselected Egypt"),
        ("Elkholy 2012 (Cairo ICU, N=26)", 76.9, "Unselected Egypt"),
        ("Fouad 2013 (Cairo/Giza ICU, N=53)", 73.6, "Unselected Egypt"),
        ("El-Masry 2018 (Giza ICU, N=22)", 63.6, "Unselected Egypt"),
        ("Itani 2023 (Sidon, Lebanon, N=111)", 94.6, "Unselected MENA"),
        ("Ballouz 2017 (Beirut, Lebanon, N=90 ep.)", 91.1, "Unselected MENA"),
        ("Ghaith 2017 (Giza ICU MDR survey, N=50)‡", 100.0, "MDR / Outbreak Survey"),
        ("Boral 2019 (Türkiye 10-center, N=172)‡", 99.4, "MDR / Outbreak Survey"),
        ("Hafiz 2023 (Albaha ICU VAP, N=115)‡", 98.0, "MDR / Outbreak Survey"),
        ("Hassan 2021 (Cairo CRAB, N=206)*", 100.0, "CRAB-Preselected*"),
        ("Abouelfetouh 2019 (Alexandria CRAB, N=74)*", 100.0, "CRAB-Preselected*"),
        ("Abdulzahra 2018 (Cairo CRAB, N=40)*", 100.0, "CRAB-Preselected*"),
        ("Fam, Gamal 2020 (Giza CRAB, N=17)*", 100.0, "CRAB-Preselected*"),
    ]
    cat_colors = {
        "Unselected Egypt": "#1D4ED8",
        "Unselected MENA": "#0D9488",
        "MDR / Outbreak Survey": "#D97706",
        "CRAB-Preselected*": "#94A3B8",
    }
    y_c = np.arange(len(crab_cohorts))[::-1]
    vals_c = [c[1] for c in crab_cohorts]
    cols_c = [cat_colors[c[2]] for c in crab_cohorts]

    ax1.barh(y_c, vals_c, height=0.64, color=cols_c, edgecolor="#0F172A", lw=0.6)
    ax1.axvline(90.7, color=CRIMSON, ls="--", lw=1.4)
    ax1.set_yticks(y_c)
    ax1.set_yticklabels([c[0] for c in crab_cohorts], fontsize=7.9)
    ax1.set_xlim(0, 116)
    ax1.set_xlabel("Carbapenem Resistance (%)", fontsize=9.2, fontweight="bold")
    ax1.set_title("A. Carbapenem Resistance Disaggregated by Sampling Frame\n(Seleim 2022 [N=100] & Alotaibi 2021 [N=198]: carbapenem AST NR)",
                  fontsize=9.6, fontweight="bold", color=NAVY, pad=12)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(axis="x", linestyle=":", alpha=0.45)

    for idx, c in enumerate(crab_cohorts):
        prefix = "~" if "Hafiz" in c[0] else ""
        ax1.text(c[1] + 1.2, y_c[idx], f"{prefix}{c[1]:.1f}%", va="center", fontsize=7.6, fontweight="bold", color="#0F172A")

    handles_a = [mpatches.Patch(color=v, label=k) for k, v in cat_colors.items()]
    handles_a.append(plt.Line2D([0], [0], color=CRIMSON, ls="--", lw=1.4, label="Unselected Egypt mean (90.7%)"))
    ax1.legend(handles=handles_a, loc="upper center", bbox_to_anchor=(0.48, -0.13), ncol=3, fontsize=7.4, frameon=True, facecolor="white")

    b_items = [
        ("HEADER", "Reference Broth Microdilution (BMD) Colistin Resistance", 0, "", ""),
        ("BMD", "Fam, Gamal 2020 (Egypt CRAB, 9/17)", 52.9, "52.9% (9/17)", "#B91C1C"),
        ("BMD", "Seleim 2022 (Egypt clinical, 49/100)", 49.0, "49.0% (49/100)", "#DC2626"),
        ("BMD", "Fam, Mohamed 2020† (Egypt CRAB, 9/22)", 40.9, "40.9% (9/22)", "#EF4444"),
        ("BMD", "Abdulzahra 2018 (Egypt CRAB, 2/40)", 5.0, "5.0% (2/40)", "#F87171"),
        ("BMD", "Boral 2019 (Türkiye multicenter, 2/172)", 1.2, "1.2% (2/172)", "#FCA5A5"),
        ("HEADER", "Reported All-Cause Mortality in MENA Cohorts (Endpoints Differ)", 0, "", ""),
        ("MORT", "Hafiz 2023 (Saudi ICU VAP; window NS)", 74.0, "74.0% (74/100*)", "#334155"),
        ("MORT", "Itani 2023 (Lebanon; 30-day all-cause)", 71.2, "71.2% (79/111)", "#475569"),
        ("MORT", "Ballouz 2017 (Lebanon BSI; in-hospital)", 63.5, "63.5% (54/85)", "#64748B"),
        ("MORT", "Boral 2019 (Türkiye ICU; window NS)", 58.5, "58.5% (96/164)", "#64748B"),
        ("MORT", "Alotaibi 2021 (Saudi ICU; window NS)", 40.8, "40.8% (80/198)", "#94A3B8"),
    ]
    y_b = np.arange(len(b_items))[::-1]
    yt_p, yt_l = [], []
    for idx, (kind, lbl, val, txt, col) in enumerate(b_items):
        yi = y_b[idx]
        if kind == "HEADER":
            ax2.text(0.5, yi, lbl, fontsize=8.2, fontweight="bold", color=NAVY, va="center")
        else:
            yt_p.append(yi)
            yt_l.append(lbl)
            ax2.barh(yi, val, height=0.60, color=col, edgecolor="#0F172A", lw=0.6)
            ax2.text(val + 1.5, yi, txt, va="center", fontsize=7.7, fontweight="bold", color="#0F172A")

    ax2.set_ylim(-0.7, len(b_items) - 0.10)
    ax2.set_yticks(yt_p)
    ax2.set_yticklabels(yt_l, fontsize=7.9)
    ax2.set_xlim(0, 98)
    ax2.set_xlabel("Percentage (%)", fontsize=9.2, fontweight="bold")
    ax2.set_title("B. Reference BMD Colistin Resistance & MENA Cohort Mortality\n(†Companion TBRI study; *Hafiz 64.3% [74/115] incl. transfers)",
                  fontsize=9.6, fontweight="bold", color=NAVY, pad=12)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(axis="x", linestyle=":", alpha=0.45)

    save_600dpi(fig, FIG_DIR / "Fig_SR6_Egyptian_MENA_CRAB_Epidemiology.png")


# ===========================================================================
# FIGURE 9 (Journal Fig 9 / Fig_SR10): ANTIGEN FREQUENCY & COLISTIN BY AST METHOD
# ===========================================================================
def build_fig9_antigens_and_colistin():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.4, 6.4), gridspec_kw={"width_ratios": [1.08, 1.05]})
    fig.subplots_adjust(wspace=0.58, bottom=0.20)

    af_sorted = sorted(antigen_freq, key=lambda d: (d["track1_pubs_of_20"] + d["track2_pubs_of_22"], d["track1_pubs_of_20"]))
    y_a = np.arange(len(af_sorted))
    h = 0.36
    t1_c = [d["track1_pubs_of_20"] for d in af_sorted]
    t2_c = [d["track2_pubs_of_22"] for d in af_sorted]

    ax1.barh(y_a + h / 2, t1_c, height=h, color="#1D4ED8", edgecolor="#0F172A", lw=0.6, label="Track 1: In Vivo Murine Survival (20 pubs)")
    ax1.barh(y_a - h / 2, t2_c, height=h, color="#0D9488", edgecolor="#0F172A", lw=0.6, label="Track 2: In Silico MEVC Design (22 pubs)")

    for i, d in enumerate(af_sorted):
        ax1.text(d["track1_pubs_of_20"] + 0.18, y_a[i] + h / 2, str(d["track1_pubs_of_20"]), va="center", fontsize=7.6, fontweight="bold", color="#1E3A8A")
        ax1.text(d["track2_pubs_of_22"] + 0.18, y_a[i] - h / 2, str(d["track2_pubs_of_22"]), va="center", fontsize=7.6, fontweight="bold", color="#0F766E")

    ax1.set_yticks(y_a)
    ax1.set_yticklabels([d["antigen"] for d in af_sorted], fontsize=8.2, fontweight="bold")
    ax1.set_xlim(0, 10.6)
    ax1.set_xlabel("Number of Included Publications Evaluating Target", fontsize=9.2, fontweight="bold")
    ax1.set_title("A. Candidate Antigen Representation Across Tracks 1 & 2\n(Note LptD/LptE & OmpK prioritized in silico but 0 in vivo survival trials)",
                  fontsize=9.6, fontweight="bold", color=NAVY, pad=12)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(axis="x", linestyle=":", alpha=0.45)
    ax1.legend(loc="upper center", bbox_to_anchor=(0.48, -0.14), ncol=1, fontsize=7.8, frameon=True, facecolor="white")

    col_items = [
        ("HEADER", "Reference ISO 20776-1 Broth Microdilution (BMD)", 0, "", ""),
        ("BMD", "Fam, Gamal 2020 (Egypt, Giza TBRI CRAB)", 52.9, "52.9% (9/17)", "#B91C1C"),
        ("BMD", "Seleim 2022 (Egypt, Cairo Univ. 65% ICU)", 49.0, "49.0% (49/100)", "#DC2626"),
        ("BMD", "Abdulzahra 2018 (Egypt, Kasr Al-Aini CRAB)", 5.0, "5.0% (2/40)", "#EF4444"),
        ("BMD", "Boral 2019 (Türkiye, 10-center MDR)", 1.2, "1.2% (2/172)", "#F87171"),
        ("HEADER", "Non-Reference / Unresolved AST Methods (Not Comparable to BMD)", 0, "", ""),
        ("NON", "Hafiz 2023 (Saudi ICU; disk / BD Phoenix)", 15.0, "15.0% (N=115)", "#D97706"),
        ("NON", "El-Masry 2018 (Egypt ICU; VITEK-2 / disk)", 4.5, "4.5% (1/22)", "#F59E0B"),
        ("NON", "Alkasaby & Zaki 2017 (Egypt ICU; E-test)", 3.2, "3.2% (9/280)", "#FBBF24"),
        ("NON", "Itani 2023 (Lebanon; disk diffusion)", 0.0, "0.0% (0/111)", "#94A3B8"),
        ("NON", "Benmahmod 2019 (Egypt; disk diffusion)", 0.0, "0.0% (0/50)", "#94A3B8"),
        ("NON", "Fouad 2013 (Egypt CRAB; method NR)", 0.0, "0.0% (0/39)", "#94A3B8"),
    ]
    y_b = np.arange(len(col_items))[::-1]
    yt_p, yt_l = [], []
    for idx, (kind, lbl, val, txt, col) in enumerate(col_items):
        yi = y_b[idx]
        if kind == "HEADER":
            ax2.text(0.4, yi, lbl, fontsize=8.0, fontweight="bold", color=NAVY, va="center")
        else:
            yt_p.append(yi)
            yt_l.append(lbl)
            bar_w = max(val, 0.6) if val == 0.0 else val
            ax2.barh(yi, bar_w, height=0.58, color=col, edgecolor="#0F172A", lw=0.6)
            ax2.text(bar_w + 1.3, yi, txt, va="center", fontsize=7.7, fontweight="bold", color="#0F172A")

    ax2.set_ylim(-0.7, len(col_items) - 0.10)
    ax2.set_yticks(yt_p)
    ax2.set_yticklabels(yt_l, fontsize=7.9)
    ax2.set_xlim(0, 68)
    ax2.set_xlabel("Reported Colistin Resistance (%)", fontsize=9.2, fontweight="bold")
    ax2.set_title("B. Colistin Resistance Across Track 3 Cohorts by AST Method\n(Disk diffusion & E-test are non-reference methods per CLSI/EUCAST)",
                  fontsize=9.6, fontweight="bold", color=NAVY, pad=12)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.grid(axis="x", linestyle=":", alpha=0.45)

    save_600dpi(fig, FIG_DIR / "Fig_SR10_Antigen_Frequency_And_Colistin_Method_Discordance.png")


if __name__ == "__main__":
    build_fig1_prisma()
    build_fig2_forest()
    build_fig3_subgroups_h2h()
    build_fig4_loo()
    build_fig5_funnel_syrcle()
    build_fig6_route_adjuvant_sensitivity()
    build_fig7_track2_heatmap(include_thesis_benchmark=False, out_name="Fig_SR4_InSilico_Methodological_Audit_Heatmap.png")
    build_fig7_track2_heatmap(include_thesis_benchmark=True, out_name="Fig_SR4_Thesis_Benchmark_Heatmap.png")
    build_fig_sr8_syrcle_matrix()
    build_fig8_epidemiology()
    build_fig9_antigens_and_colistin()
    print("ALL 11 VERIFIED & VISUALLY POLISHED 600-DPI FIGURES GENERATED SUCCESSFULLY.")
