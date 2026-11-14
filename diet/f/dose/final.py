#!/usr/bin/env python3
"""Assemble DOSE manuscript figures from computed aggregates and native figures.

Usage: python diet/f/dose/final.py --output /path/to/analysis/diet
No participant records are required. No estimates are copied from the manuscript.
Native PDF panels retain their vector content when PyMuPDF composes the pages.
"""
from __future__ import annotations

import argparse
import csv
import math
import shutil
from pathlib import Path

import fitz
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def number(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None


def count(value):
    x = number(value)
    return f"{int(round(x)):,}" if x is not None else "not available"


def extension(stem: Path, suffix: str):
    """Append a format extension without dropping dotted scientific figure names."""
    return Path(str(stem) + suffix)


def choose(stem: Path):
    # Vector source first; no screenshots or crops of an original panel.
    for suffix in (".pdf", ".png"):
        path = extension(stem, suffix)
        if path.is_file() and path.stat().st_size:
            return path
    return None


def clear_generated(stem: Path):
    for suffix in (".pdf", ".png"):
        path = extension(stem, suffix)
        if path.is_file():
            path.unlink()


def render_pdf(path: Path):
    with fitz.open(path) as document:
        document[0].get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False).save(path.with_suffix(".png"))


def add_asset(page, rect, path: Path):
    if path.suffix.lower() == ".pdf":
        with fitz.open(path) as source:
            if not source.page_count:
                raise ValueError(f"Empty PDF panel: {path}")
            page.show_pdf_page(rect, source, 0, keep_proportion=True)
    else:
        page.insert_image(rect, filename=str(path), keep_proportion=True)


def aspect(path: Path):
    if path.suffix.lower() == ".pdf":
        with fitz.open(path) as document:
            rect = document[0].rect
            return rect.height / rect.width
    pixmap = fitz.Pixmap(str(path))
    return pixmap.height / pixmap.width


def save_native(source, stem: Path, figure: str):
    if source is None:
        clear_generated(stem)
        return status(figure, "input_missing", stem, ["native figure"], [])
    if source.suffix == ".pdf":
        shutil.copy2(source, extension(stem, ".pdf"))
    else:
        document = fitz.open()
        width = 900
        page = document.new_page(width=width, height=width * aspect(source))
        add_asset(page, page.rect, source)
        document.save(extension(stem, ".pdf"), deflate=True)
        document.close()
    render_pdf(extension(stem, ".pdf"))
    return status(figure, "complete", stem, [], [source])


def status(figure, state, stem, missing, sources):
    return dict(figure=figure, status=state,
                pdf=str(extension(stem, ".pdf")) if extension(stem, ".pdf").is_file() else "",
                png=str(extension(stem, ".png")) if extension(stem, ".png").is_file() else "",
                missing_panels="; ".join(missing), sources="; ".join(map(str, sources)))


def compose(stem, figure, title, panels, rows, notes=()):
    """Rows contain panel keys and width proportions; all content is fitted whole."""
    available = [p["source"] for p in panels.values() if p["source"]]
    missing = [key for key, p in panels.items() if p["source"] is None]
    missing.extend(notes)
    if not available:
        clear_generated(stem)
        return status(figure, "input_missing", stem, missing, [])
    width, margin, gap, label_height = 1100.0, 25.0, 14.0, 25.0
    usable = width - 2 * margin
    row_heights = []
    for row in rows:
        widths = [(usable - gap * (len(row) - 1)) * fraction for _, fraction in row]
        heights = []
        for (key, _), cell_width in zip(row, widths):
            source = panels[key]["source"]
            heights.append(cell_width * aspect(source) if source else cell_width * .62)
        row_heights.append(max(heights) + label_height)
    title_height = 77 if missing else 58
    height = title_height + sum(row_heights) + gap * (len(rows) - 1) + margin
    document = fitz.open()
    page = document.new_page(width=width, height=height)
    page.insert_text((margin, 27), title, fontsize=18, color=(.10, .18, .25))
    if missing:
        page.insert_textbox(fitz.Rect(margin, 38, width - margin, 72),
                            "PARTIAL BUILD - unavailable panels or estimates: " + ", ".join(
                                key for key, p in panels.items() if p["source"] is None) +
                            (". Some panel rows are unavailable; see composite_status.tsv." if notes else ""),
                            fontsize=10, color=(.62, .22, .12))
    y = title_height
    for row, row_height in zip(rows, row_heights):
        x = margin
        for key, fraction in row:
            cell_width = (usable - gap * (len(row) - 1)) * fraction
            panel = panels[key]
            page.insert_text((x + 3, y + 15), f"{key}  {panel['label']}", fontsize=12,
                             color=(.10, .18, .25))
            rect = fitz.Rect(x, y + label_height, x + cell_width, y + row_height)
            if panel["source"]:
                panel_height = cell_width * aspect(panel["source"])
                asset_rect = fitz.Rect(x, y + label_height, x + cell_width, y + label_height + panel_height)
                add_asset(page, asset_rect, panel["source"])
            else:
                page.draw_rect(rect, color=(.83, .85, .88), fill=(.96, .97, .98))
                page.insert_textbox(fitz.Rect(x + 15, y + label_height + 30,
                                              x + cell_width - 15, y + row_height - 15),
                                    "Analysis output has not been generated.\n" + panel["label"],
                                    fontsize=12, color=(.38, .40, .43))
            x += cell_width + gap
        y += row_height + gap
    document.set_metadata({"title": title, "subject": "Computed DOSE results; whole native panels, no estimates imputed"})
    document.save(extension(stem, ".pdf"), deflate=True)
    document.close()
    render_pdf(extension(stem, ".pdf"))
    return status(figure, "partial" if missing else "complete", stem, missing, available)


def stroke_comparison(base: Path, final: Path):
    ukb_path = base / "ukb/04score/overall_stroke_comparison/UKB_overall_stroke_M2.csv"
    nh_path = base / "nhanes/NHANES_M2_results.csv"
    ukb = read_csv(ukb_path)
    nh = read_csv(nh_path)
    stem = final / "Fig3f.DOSE.stroke_comparison"
    if not ukb or not nh:
        clear_generated(stem)
        return None, ["3f UKB/NHANES aggregate comparisons"]
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.8), sharey=True)
    labels = ["UKB Q1", "UKB Q2", "UKB Q3", "UKB Q4", "UKB per SD", "",
              "NHANES T1", "NHANES T2", "NHANES T3", "NHANES per SD"]
    missing = []
    for ax, score in zip(axes, ("dose", "dash")):
        u = {r.get("Contrast"): r for r in ukb if r.get("Diet") == score and r.get("Model") == "M2"}
        h = {r.get("Contrast"): r for r in nh if r.get("Score") == score and
             r.get("Analysis") == "two_day_2003_2018" and r.get("Model") == "M2"}
        bounds = [1.0]
        for key, y in zip(("Q1", "Q2", "Q3", "Q4", "Per SD"), (0, 1, 2, 3, 4)):
            row = u.get(key)
            if not row:
                missing.append(f"3f {score} UKB {key}"); continue
            value, lower, upper = (number(row.get(v)) for v in ("Estimate", "LCI", "UCI"))
            if None in (value, lower, upper) or min(value, lower, upper) <= 0 or not lower <= value <= upper:
                missing.append(f"3f {score} UKB {key} estimate"); continue
            ax.errorbar(value, y, xerr=[[value - lower], [upper - value]], fmt="s",
                        color="#2F5D8A", capsize=2.8, ms=5.5, lw=1.25)
            bounds.extend((lower, upper))
        for key, y in zip(("T1", "T2", "T3", "Per SD"), (6, 7, 8, 9)):
            row = h.get(key)
            if not row:
                missing.append(f"3f {score} NHANES {key}"); continue
            value, lower, upper = (number(row.get(v)) for v in ("OR", "LCI", "UCI"))
            if None in (value, lower, upper) or min(value, lower, upper) <= 0 or not lower <= value <= upper:
                missing.append(f"3f {score} NHANES {key} estimate"); continue
            ax.errorbar(value, y, xerr=[[value - lower], [upper - value]], fmt="o",
                        color="#B96537", capsize=2.8, ms=5.5, lw=1.25)
            bounds.extend((lower, upper))
        ax.axvline(1, ls="--", lw=.8, color=".5")
        ax.axhline(5, lw=.6, color=".88")
        ax.set_xscale("log")
        ax.set_xlim(min(bounds) / 1.12, max(bounds) * 1.12)
        ax.set_yticks(range(len(labels)), labels, fontsize=10)
        ax.set_ylim(9.7, -.7)
        ax.set_xlabel("Hazard ratio (UKB) / odds ratio (NHANES), 95% CI", fontsize=9)
        ax.set_title(score.upper(), loc="left", fontsize=13, weight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        u_main, h_main = u.get("Per SD", {}), h.get("Per SD", {})
        ax.text(.01, -0.17,
                f"UKB: {count(u_main.get('N'))} participants, {count(u_main.get('Events'))} events\n"
                f"NHANES: {count(h_main.get('N'))} participants, {count(h_main.get('Cases'))} cases",
                transform=ax.transAxes, fontsize=8.4, va="top")
    fig.suptitle("Overall stroke: incident disease in UKB and prevalent disease in NHANES", fontsize=14, x=.06, ha="left")
    fig.text(.06, .025, "Model 2. UKB: score quartiles; NHANES 2003-2018: weighted tertiles from two valid recalls.\n"
             "Per-SD estimates use each cohort's own score SD. The estimates are displayed separately; no pooled effect is calculated.",
             fontsize=8.3, va="bottom")
    if missing:
        fig.text(.06, .915, "PARTIAL BUILD: some comparison rows are unavailable", color="#9E3A20", fontsize=9)
    fig.subplots_adjust(top=.86, bottom=.25, left=.14, right=.97, wspace=.20)
    fig.savefig(extension(stem, ".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(extension(stem, ".png"), dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return extension(stem, ".pdf"), missing


def overview(base: Path, final: Path):
    stem = final / "Fig1.DOSE.overview"
    event_path = base / "ukb/00data/tab_event.csv"
    flow_path = base / "ukb/00data/flow.csv"
    events, flow = read_csv(event_path), read_csv(flow_path)
    primary = next((r for r in events if r.get("outcome", "").lower() == "ischemic stroke"), None)
    if primary is None:
        clear_generated(stem)
        return status("Fig1", "input_missing", stem, ["computed UKB cohort counts"], [])
    nh_path = base / "nhanes/NHANES_M2_results.csv"
    nh = next((r for r in read_csv(nh_path) if r.get("Analysis") == "two_day_2003_2018" and
               r.get("Score") == "dose" and r.get("Contrast") == "Per SD"), None)
    map_path = base / "ukb/02foodwas/food_group_map.csv"
    group_count = len({r.get("Food.Group") for r in read_csv(map_path) if r.get("Food.Group")})
    n, cases = count(primary.get("total")), count(primary.get("case"))
    model_path = base / "ukb/04score/pattern/pattern_continuous_10patterns.csv"
    if not model_path.is_file(): model_path = base / "ukb/04score/pattern/pattern_continuous.csv"
    fitted = next((r for r in read_csv(model_path) if r.get("diet") == "dose" and
                   r.get("outcome") == "cvd_stroke_i" and r.get("model") == "M2"), None)
    if fitted:
        primary_text = (f"Incident ischemic stroke: M2\n"
                        f"n = {count(fitted.get('N'))}; {count(fitted.get('Events'))} events\n"
                        "Cox models and sensitivity analyses")
    else:
        primary_text = (f"Eligible ischemic stroke cohort\nn = {n}; {cases} events\n"
                        "Fitted model counts pending")
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set(xlim=(0, 1), ylim=(0, 1)); ax.axis("off")
    ax.text(.025, .97, "DOSE study overview", fontsize=22, weight="bold", va="top", color="#19364B")
    def box(x, y, w, h, text, fill="#EDF3F8", fontsize=11):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.008,rounding_size=0.009",
                                   fc=fill, ec="#8A9EAC", lw=1))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize,
                color="#162F3F", linespacing=1.45)
    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14,
                                     lw=1.3, color="#7890A0"))
    flow_labels = []
    if flow:
        names = {"Any WebQ record": "Any Oxford WebQ record",
                 "Baseline diet not often varied": "Baseline diet not often varied",
                 "At least one valid recall": "At least one valid dietary recall",
                 "No prevalent stroke and valid follow-up": "No prevalent stroke; valid follow-up"}
        for row in flow[:4]:
            flow_labels.append(f"{names.get(row.get('stage'), row.get('stage', 'Cohort'))}\nn = {count(row.get('n'))}")
    else:
        flow_labels = [f"UKB analytical cohort\nn = {n}"]
    yy = [.73, .54, .35, .16] if len(flow_labels) > 1 else [.53]
    for i, (text, y) in enumerate(zip(flow_labels, yy)):
        box(.035, y, .28, .125, text, fontsize=10.5)
        if i: arrow(.175, yy[i-1], .175, y+.125)
    box(.39, .69, .26, .16, "Food-wide associations\n" +
        (f"{group_count} mapped food groups\n" if group_count else "Food groups from foods.xlsx\n") +
        "Machine-learning-assisted selection", fontsize=11)
    box(.39, .43, .26, .17, "Frozen DOSE9half scoring\n11 food criteria; 9 domains\n0-9 points in 0.5-point steps", fill="#E6EFD9")
    box(.39, .17, .26, .16, primary_text, fontsize=11)
    arrow(.315, yy[-1]+.0625, .39, .77)
    arrow(.52, .69, .52, .60); arrow(.52, .43, .52, .33)
    nh_text = (f"External NHANES 2003-2018\nn = {count(nh.get('N'))}; {count(nh.get('Cases'))} cases\nSurvey-weighted prevalent stroke"
               if nh else "External NHANES validation\nAggregate results not generated")
    box(.725, .68, .25, .17, nh_text, fill="#F5EADB", fontsize=10.5)
    box(.725, .42, .25, .17, "Outcome-wide associations\nOutcome-specific exclusions\nMultiple-testing adjustment", fontsize=10.5)
    box(.725, .16, .25, .17, "Biological correlates\nMetabolites, inflammation, proteins\nBrain structure and white matter", fontsize=10.5)
    for y in (.765, .505, .245): arrow(.65, .515, .725, y)
    missing = [] if nh else ["NHANES participant and case counts"]
    if not flow: missing.append("initial cohort screening counts")
    if not fitted: missing.append("fitted primary M2 participant and event counts")
    footer = "Screening and fitted-model counts are read separately from the current run. UKB development samples overlap; NHANES is cross-sectional."
    if not flow: footer += " Initial screening counts were not exported in this run."
    if missing: footer = "PARTIAL BUILD: some cohort counts are unavailable. " + footer
    ax.text(.035, .055, footer, fontsize=9, color="#566979", va="top", wrap=True)
    fig.savefig(extension(stem, ".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(extension(stem, ".png"), dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return status("Fig1", "partial" if missing else "complete", stem, missing,
                  [p for p in (event_path, flow_path, model_path, map_path, nh_path) if p.is_file()])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="analysis/diet root")
    args = parser.parse_args(argv)
    output = args.output.resolve()
    base, final = output / "dose", output / "final/dose"
    final.mkdir(parents=True, exist_ok=True)
    rows = [overview(base, final)]
    rows.append(save_native(choose(base / "ukb/02foodwas/group/figures/food_group_MYTH_circular_ref_first_ordered"),
                            final / "Fig2.DOSE.foodwas", "Fig2"))
    comparison, comparison_missing = stroke_comparison(base, final)
    panels = {
        "a": dict(label="Food prioritization and cumulative association", source=choose(base / "ukb/03ml/figures/Fig_ML_pAB_Cox")),
        "b": dict(label="DOSE scoring criteria", source=choose(base / "ukb/04score/figures/Fig_DOSE9HALF")),
        "c": dict(label="Food intake and score distribution", source=choose(base / "ukb/04score/figures/Fig_DOSE_adherence")),
        "d": dict(label="Ten dietary patterns and ischemic stroke", source=choose(base / "ukb/04score/pattern/figures/Fig_Pattern_Stroke_Forest_M2")),
        "e": dict(label="DOSE quartiles and spline association", source=choose(base / "ukb/04score/pattern/figures/Fig_DOSE_Quartiles_RCS_M2")),
        "f": dict(label="UKB and NHANES overall stroke", source=comparison),
    }
    rows.append(compose(final / "Fig3.DOSE.score", "Fig3", "DOSE construction and stroke associations", panels,
                        [[("a", .43), ("b", .57)], [("c", .5), ("d", .5)], [("e", .5), ("f", .5)]],
                        comparison_missing if comparison else ()))
    rows.append(save_native(choose(base / "ukb/04score/otherOutcomes/figures/Fig_DOSE_OutcomeWide_M1_ChapterForest_P10P90"),
                            final / "Fig4.DOSE.outcomes", "Fig4"))
    mechanism = {
        "a-j": dict(label="Metabolomic, inflammatory and proteomic correlates", source=choose(base / "ukb/05mechanism/summary_figures/Fig_s4_mechanism_summary")),
        "k": dict(label="Brain structural correlates", source=choose(base / "ukb/05mechanism/04mri/figures/main/Figure4_MRI")),
    }
    rows.append(compose(final / "Fig5.DOSE.mechanisms", "Fig5", "Biological correlates of DOSE", mechanism,
                        [[("a-j", 1.)], [("k", 1.)]]))
    path = final / "composite_status.tsv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["figure", "status", "pdf", "png", "missing_panels", "sources"], delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    for row in rows:
        print(f"{row['figure']}: {row['status']}" + (f" ({row['missing_panels']})" if row['missing_panels'] else ""))
    print(f"Saved figure status: {path}")


if __name__ == "__main__":
    main()
