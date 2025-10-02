# coding: utf-8
"""Proteome Task 2D analysis script."""
import os
import sys
import logging
import traceback
import math
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any

OUTDIR = r"C:/Users/23102/Desktop/蛋白质组学output"
RESULT_DIR = os.path.join(OUTDIR, "2D")
LOG_FILE = os.path.join(RESULT_DIR, "2D_task.log")

REQUIRED_PACKAGES = [
    ("pandas", None),
    ("numpy", None),
    ("scipy", None),
    ("seaborn", None),
    ("matplotlib", None),
    ("statsmodels", None),
    ("pingouin", None),
    ("openpyxl", None),
    ("xlsxwriter", None),
]

KEY_PROTEINS = [
    "ENPP2",
    "PLPP3",
    "RAF1",
    "MAP2K1",
    "MAPK1",
    "MAPK3",
    "RPS6KA1",
    "CCL2",
    "LYZ",
    "S100A8",
    "S100A9",
    "ITGAM",
    "NFKBIA",
    "IKBKB",
]

KEY_PATHWAYS = [
    "CUSTOM_ERK_DOWNSTREAM",
    "GO_ERK1_AND_ERK2_CASCADE",
    "HALLMARK_KRAS_SIGNALING_UP",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_COMPLEMENT",
]

os.makedirs(RESULT_DIR, exist_ok=True)


def setup_logging() -> logging.Logger:
    logger = logging.getLogger("proteome_task2D")
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    if not logger.handlers:
        file_handler = logging.FileHandler(LOG_FILE, mode="a", encoding="utf-8-sig")
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    return logger


def ensure_packages(logger: logging.Logger):
    import importlib
    import subprocess

    for pkg, import_name in REQUIRED_PACKAGES:
        module_name = import_name or pkg
        try:
            importlib.import_module(module_name)
            logger.info("Package '%s' is available.", pkg)
        except ImportError:
            logger.warning("Package '%s' not found. Attempting installation.", pkg)
            try:
                subprocess.check_call([sys.executable, "-m", "pip", "install", pkg])
                importlib.import_module(module_name)
                logger.info("Package '%s' installed successfully.", pkg)
            except Exception as exc:
                logger.error("Failed to install package '%s': %s", pkg, exc)
                logger.error(traceback.format_exc())


def norm_name(value: Any) -> str:
    if value is None:
        value = ""
    if isinstance(value, float) and math.isnan(value):
        value = ""
    value = str(value).strip().lower()
    for ch in [" ", ".", "-", "|", "/", "\\", ":", ";", "[", "]", "(", ")", "{", "}"]:
        value = value.replace(ch, "_")
    while "__" in value:
        value = value.replace("__", "_")
    value = value.strip("_")
    if not value:
        value = "unnamed"
    return value


def normalize_axis(values: List[Any]) -> Tuple[List[str], Dict[str, str]]:
    normalized = []
    mapping: Dict[str, str] = {}
    counts: Dict[str, int] = {}
    for val in values:
        original = str(val)
        key = norm_name(original)
        if key in counts:
            counts[key] += 1
            key = f"{key}_{counts[key]}"
        else:
            counts[key] = 0
        normalized.append(key)
        mapping[key] = original
    return normalized, mapping


def determine_orientation(df1, df2, logger: logging.Logger) -> Tuple[str, str]:
    def norm_set(values):
        return set(norm_name(v) for v in values)

    options1 = {"columns": norm_set(df1.columns), "rows": norm_set(df1.index)}
    options2 = {"columns": norm_set(df2.columns), "rows": norm_set(df2.index)}
    best_pair = ("columns", "columns")
    best_overlap = -1
    for o1, set1 in options1.items():
        for o2, set2 in options2.items():
            overlap = len(set1 & set2)
            if overlap > best_overlap:
                best_pair = (o1, o2)
                best_overlap = overlap
    logger.info("Orientation chosen: df1 -> %s, df2 -> %s (overlap %d)", best_pair[0], best_pair[1], best_overlap)
    return best_pair


def apply_orientation(df, orientation: str, label: str, logger: logging.Logger):
    try:
        if orientation == "rows":
            logger.info("Transposing %s to align samples as columns.", label)
            df = df.transpose()
    except Exception as exc:
        logger.error("Failed to transpose %s: %s", label, exc)
        logger.error(traceback.format_exc())
    return df


def maybe_transpose(df, label: str, orientation: str, logger: logging.Logger):
    df = apply_orientation(df, orientation, label, logger)
    try:
        normalized_columns, col_map = normalize_axis(list(df.columns))
        df.columns = normalized_columns
        normalized_index, idx_map = normalize_axis(list(df.index))
        df.index = normalized_index
        logger.info("Normalized axis for %s (rows=%d, cols=%d)", label, len(df.index), len(df.columns))
        return df, col_map, idx_map
    except Exception as exc:
        logger.error("Failed to normalize axis for %s: %s", label, exc)
        logger.error(traceback.format_exc())
        return df, {}, {}


def load_csv_safe(path: str, logger: logging.Logger):
    import pandas as pd

    try:
        df = pd.read_csv(path, encoding="utf-8-sig")
        logger.info("Loaded %s with shape %s", path, df.shape)
        return df
    except UnicodeDecodeError:
        df = pd.read_csv(path, encoding="utf-8")
        logger.info("Loaded %s with utf-8 encoding", path)
        return df
    except Exception as exc:
        logger.error("Failed to load %s: %s", path, exc)
        logger.error(traceback.format_exc())
        raise


def fuzzy_map(target_list: List[str], column_names: List[str], logger: logging.Logger) -> Dict[str, Optional[str]]:
    from difflib import SequenceMatcher

    normalized_candidates = {}
    for name in column_names:
        normalized_candidates[norm_name(name).replace("_", "")] = name

    mapping: Dict[str, Optional[str]] = {}
    for target in target_list:
        target_norm = norm_name(target).replace("_", "")
        if target_norm in normalized_candidates:
            mapping[target] = normalized_candidates[target_norm]
            continue
        best_match = None
        best_ratio = 0.0
        for cand_norm, original in normalized_candidates.items():
            ratio = SequenceMatcher(None, target_norm, cand_norm).ratio()
            if ratio > best_ratio:
                best_ratio = ratio
                best_match = original
        if best_match and best_ratio >= 0.7:
            mapping[target] = best_match
            logger.info("Fuzzy matched %s -> %s (%.2f)", target, best_match, best_ratio)
        else:
            mapping[target] = None
            logger.warning("Failed to match %s (best score %.2f)", target, best_ratio)
    return mapping


def starify(q: Optional[float]) -> str:
    if q is None:
        return ""
    if isinstance(q, float) and math.isnan(q):
        return ""
    if q < 0.05:
        return "**"
    if q < 0.1:
        return "*"
    return ""


def adjust_fdr(pvalues: List[float]) -> List[float]:
    from statsmodels.stats.multitest import multipletests

    if not pvalues:
        return []
    safe = [1.0 if (p is None or (isinstance(p, float) and math.isnan(p))) else p for p in pvalues]
    try:
        _, qvals, _, _ = multipletests(safe, method="fdr_bh")
        return list(qvals)
    except Exception as exc:
        logging.getLogger("proteome_task2D").error("FDR adjustment failed: %s", exc)
        logging.getLogger("proteome_task2D").error(traceback.format_exc())
        return safe


def build_covars(metadata, columns: List[str]):
    import pandas as pd

    covar_frames = []
    covar_cols: List[str] = []
    for col in columns:
        if col in metadata.columns:
            dummy = pd.get_dummies(metadata[col].astype(str), prefix=col, drop_first=True)
            if not dummy.empty:
                covar_frames.append(dummy)
                covar_cols.extend(list(dummy.columns))
    if covar_frames:
        covars = pd.concat(covar_frames, axis=1)
        return covar_cols, covars
    return [], None


def write_csv(df, path: str, logger: logging.Logger):
    try:
        df.to_csv(path, encoding="utf-8-sig")
        logger.info("Saved CSV %s", path)
    except Exception as exc:
        logger.error("Failed to write CSV %s: %s", path, exc)
        logger.error(traceback.format_exc())


def save_excel(dfs: Dict[str, 'pd.DataFrame'], path: str, logger: logging.Logger):
    import pandas as pd

    try:
        with pd.ExcelWriter(path, engine="xlsxwriter") as writer:
            for sheet, data in dfs.items():
                data.to_excel(writer, sheet_name=sheet[:31], index=False)
        logger.info("Saved Excel %s", path)
    except Exception as exc:
        logger.error("Failed to write Excel %s: %s", path, exc)
        logger.error(traceback.format_exc())


def configure_matplotlib(logger: logging.Logger):
    try:
        import matplotlib
        matplotlib.rcParams["axes.unicode_minus"] = False
        for font in ["SimHei", "Microsoft YaHei", "Arial Unicode MS"]:
            try:
                matplotlib.rcParams["font.sans-serif"] = [font]
                break
            except Exception:
                continue
        import matplotlib.pyplot as plt
        return plt
    except Exception as exc:
        logger.error("Failed to configure matplotlib: %s", exc)
        logger.error(traceback.format_exc())
        return None


def compute_spearman_groupwise(proteins_mat, scores_mat, metadata, matched_proteins, matched_pathways, logger: logging.Logger):
    from scipy.stats import spearmanr
    import pandas as pd

    records = []
    for group in ["BA", "Ctrl"]:
        group_samples = metadata.index[metadata["group"] == group].tolist()
        if not group_samples:
            logger.warning("Group %s has no samples", group)
            continue
        for protein_label, protein_col in matched_proteins.items():
            if protein_col is None:
                continue
            for pathway_label, pathway_col in matched_pathways.items():
                if pathway_col is None:
                    continue
                try:
                    df = pd.DataFrame(
                        {
                            "x": proteins_mat.loc[group_samples, protein_col],
                            "y": scores_mat.loc[group_samples, pathway_col],
                        }
                    ).dropna()
                    n = len(df)
                    if n < 3:
                        logger.info("Insufficient data for %s-%s in %s (n=%d)", protein_col, pathway_col, group, n)
                        continue
                    rho, pval = spearmanr(df["x"], df["y"])
                    records.append(
                        {
                            "protein": protein_col,
                            "pathway": pathway_col,
                            "group": group,
                            "n": n,
                            "rho": rho,
                            "p": pval,
                        }
                    )
                except Exception as exc:
                    logger.error("Spearman failed for %s-%s (%s): %s", protein_col, pathway_col, group, exc)
                    logger.error(traceback.format_exc())
    import pandas as pd

    df = pd.DataFrame(records)
    if not df.empty:
        for group, subset in df.groupby("group"):
            qvals = adjust_fdr(subset["p"].tolist())
            df.loc[subset.index, "q"] = qvals
    return df


def create_heatmap(data_df, value_field: str, title: str, output_path: str, cmap, logger: logging.Logger):
    if data_df is None or data_df.empty:
        logger.warning("No data for heatmap %s", title)
        return
    import pandas as pd

    pivot = data_df.pivot_table(index="protein", columns="pathway", values=value_field)
    q_pivot = data_df.pivot_table(index="protein", columns="pathway", values="q")
    plt = configure_matplotlib(logger)
    if plt is None:
        return
    import seaborn as sns

    try:
        sns.set(style="whitegrid")
        fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(pivot.columns)), max(4, 0.5 * len(pivot.index))))
        sns.heatmap(
            pivot,
            cmap=cmap if cmap is not None else "coolwarm",
            vmin=-1,
            vmax=1,
            center=0,
            ax=ax,
            linewidths=0.5,
            linecolor="white",
            cbar_kws={"label": value_field},
        )
        ax.set_title(title)
        for i, protein in enumerate(pivot.index):
            for j, pathway in enumerate(pivot.columns):
                try:
                    q_val = q_pivot.loc[protein, pathway]
                except KeyError:
                    q_val = None
                stars = starify(q_val)
                if stars:
                    ax.text(j + 0.5, i + 0.5, stars, ha="center", va="center", color="black", fontsize=12)
        plt.tight_layout()
        fig.savefig(output_path, dpi=300)
        plt.close(fig)
        logger.info("Saved heatmap %s", output_path)
    except Exception as exc:
        logger.error("Failed to create heatmap %s: %s", title, exc)
        logger.error(traceback.format_exc())


def compute_partial_correlations(proteins_mat, scores_mat, metadata, matched_proteins, matched_pathways, logger: logging.Logger):
    import pandas as pd
    import pingouin as pg

    configs = [
        ("all_partial", metadata.index.tolist(), ["group", "tissue"]),
        ("BA_partial", metadata.index[metadata["group"] == "BA"].tolist(), ["tissue"]),
        ("Ctrl_partial", metadata.index[metadata["group"] == "Ctrl"].tolist(), ["tissue"]),
    ]
    result_dfs: Dict[str, pd.DataFrame] = {}
    for name, samples, covar_cols in configs:
        if len(samples) < 6:
            logger.warning("Skipping %s due to insufficient samples (%d)", name, len(samples))
            continue
        records = []
        meta_sub = metadata.loc[samples]
        covar_names, covars = build_covars(meta_sub, covar_cols)
        for protein_label, protein_col in matched_proteins.items():
            if protein_col is None:
                continue
            for pathway_label, pathway_col in matched_pathways.items():
                if pathway_col is None:
                    continue
                try:
                    df = pd.DataFrame(
                        {
                            "x": proteins_mat.loc[samples, protein_col],
                            "y": scores_mat.loc[samples, pathway_col],
                        }
                    )
                    if covars is not None:
                        df = pd.concat([df, covars], axis=1)
                    df = df.dropna()
                    if len(df) < 3:
                        logger.info("Insufficient data for partial corr %s-%s in %s", protein_col, pathway_col, name)
                        continue
                    covar_arg = covar_names if covars is not None and covar_names else None
                    res = pg.partial_corr(data=df, x="x", y="y", covar=covar_arg, method="spearman")
                    rho = res.loc[0, "r"]
                    pval = res.loc[0, "p-val"]
                    records.append(
                        {
                            "protein": protein_col,
                            "pathway": pathway_col,
                            "analysis": name,
                            "n": res.loc[0, "n"] if "n" in res.columns else len(df),
                            "rho": rho,
                            "p": pval,
                        }
                    )
                except Exception as exc:
                    logger.error("Partial corr failed for %s-%s in %s: %s", protein_col, pathway_col, name, exc)
                    logger.error(traceback.format_exc())
        df_out = pd.DataFrame(records)
        if not df_out.empty:
            df_out["q"] = adjust_fdr(df_out["p"].tolist())
        result_dfs[name] = df_out
    return result_dfs


def compute_linear_models(proteins_mat, scores_mat, metadata, matched_proteins, matched_pathways, logger: logging.Logger):
    import pandas as pd
    import statsmodels.formula.api as smf

    records = []
    for protein_label, protein_col in matched_proteins.items():
        if protein_col is None:
            continue
        for pathway_label, pathway_col in matched_pathways.items():
            if pathway_col is None:
                continue
            try:
                df = pd.DataFrame(
                    {
                        "score": scores_mat[pathway_col],
                        "protein_raw": proteins_mat[protein_col],
                        "group": metadata["group"],
                        "tissue": metadata["tissue"],
                    }
                ).dropna()
                if len(df) < 4:
                    logger.info("Skipping linear model %s-%s due to n=%d", protein_col, pathway_col, len(df))
                    continue
                df["zProtein"] = (df["protein_raw"] - df["protein_raw"].mean()) / (df["protein_raw"].std(ddof=0) + 1e-9)
                model = smf.ols("score ~ zProtein + C(group) + C(tissue)", data=df).fit()
                robust = model.get_robustcov_results(cov_type="HC3")
                beta = robust.params.get("zProtein", float("nan"))
                tval = robust.tvalues.get("zProtein", float("nan"))
                pval = robust.pvalues.get("zProtein", float("nan"))
                records.append(
                    {
                        "protein": protein_col,
                        "pathway": pathway_col,
                        "n": len(df),
                        "beta_protein": beta,
                        "t": tval,
                        "p": pval,
                        "adj_R2": model.rsquared_adj,
                    }
                )
            except Exception as exc:
                logger.error("Linear model failed for %s-%s: %s", protein_col, pathway_col, exc)
                logger.error(traceback.format_exc())
    import pandas as pd

    lm_df = pd.DataFrame(records)
    if not lm_df.empty:
        for pathway, subset in lm_df.groupby("pathway"):
            qvals = adjust_fdr(subset["p"].tolist())
            lm_df.loc[subset.index, "q"] = qvals
    return lm_df


def create_scatter_plots(top_pairs, proteins_mat, scores_mat, metadata, output_pdf, logger: logging.Logger):
    if not top_pairs:
        logger.warning("No pairs for scatter plots.")
        return
    plt = configure_matplotlib(logger)
    if plt is None:
        return
    import seaborn as sns
    from matplotlib.backends.backend_pdf import PdfPages
    import numpy as np
    from scipy.stats import spearmanr

    try:
        sns.set(style="whitegrid")
        with PdfPages(output_pdf) as pdf:
            for entry in top_pairs:
                protein = entry["protein"]
                pathway = entry["pathway"]
                fig, ax = plt.subplots(figsize=(6, 5))
                df = proteins_mat[[protein]].join(scores_mat[[pathway]])
                df = df.join(metadata[["group"]])
                df.columns = ["protein", "score", "group"]
                df = df.dropna()
                sns.scatterplot(data=df, x="protein", y="score", hue="group", ax=ax, palette="Set1")
                if len(df) >= 2:
                    slope, intercept = np.polyfit(df["protein"], df["score"], 1)
                    xs = np.linspace(df["protein"].min(), df["protein"].max(), 100)
                    ax.plot(xs, slope * xs + intercept, linestyle="--", color="black", label="Fit")
                    rho, pval = spearmanr(df["protein"], df["score"])
                    ax.text(0.05, 0.95, f"Spearman rho={rho:.3f}\np={pval:.3g}", transform=ax.transAxes, va="top")
                ax.set_title(f"{protein} vs {pathway}")
                ax.legend()
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
        logger.info("Saved scatter plots to %s", output_pdf)
    except Exception as exc:
        logger.error("Failed to create scatter plots: %s", exc)
        logger.error(traceback.format_exc())


def integrate_evidence(group_df, partial_dfs, lm_df, logger: logging.Logger):
    import pandas as pd

    records: Dict[Tuple[str, str], Dict[str, Any]] = {}

    def ensure_record(protein, pathway):
        key = (protein, pathway)
        if key not in records:
            records[key] = {
                "protein": protein,
                "pathway": pathway,
                "flag_group_corr": 0,
                "flag_partial": 0,
                "flag_lm": 0,
                "best_group_rho": None,
                "best_group_q": None,
                "best_partial_rho": None,
                "partial_q": None,
                "beta_lm": None,
                "lm_q": None,
            }
        return records[key]

    if group_df is not None and not group_df.empty:
        sorted_group = group_df.sort_values("q")
        for _, row in sorted_group.iterrows():
            rec = ensure_record(row["protein"], row["pathway"])
            qval = row.get("q")
            if qval is not None and not (isinstance(qval, float) and math.isnan(qval)) and qval < 0.1:
                rec["flag_group_corr"] = 1
            rho = row.get("rho")
            if rho is not None:
                current = rec.get("best_group_rho")
                if current is None or abs(rho) > abs(current):
                    rec["best_group_rho"] = rho
                    rec["best_group_q"] = qval
    if partial_dfs:
        df_all = partial_dfs.get("all_partial")
        if df_all is not None and not df_all.empty:
            for _, row in df_all.sort_values("q").iterrows():
                rec = ensure_record(row["protein"], row["pathway"])
                qval = row.get("q")
                if qval is not None and not (isinstance(qval, float) and math.isnan(qval)) and qval < 0.1:
                    rec["flag_partial"] = 1
                rho = row.get("rho")
                if rho is not None:
                    current = rec.get("best_partial_rho")
                    if current is None or abs(rho) > abs(current):
                        rec["best_partial_rho"] = rho
                        rec["partial_q"] = qval
    if lm_df is not None and not lm_df.empty:
        for _, row in lm_df.iterrows():
            rec = ensure_record(row["protein"], row["pathway"])
            qval = row.get("q")
            if qval is not None and not (isinstance(qval, float) and math.isnan(qval)) and qval < 0.1:
                rec["flag_lm"] = 1
            rec["beta_lm"] = row.get("beta_protein")
            rec["lm_q"] = qval
    data = []
    for rec in records.values():
        effects = []
        for val in [rec.get("best_group_rho"), rec.get("best_partial_rho"), rec.get("beta_lm")]:
            if val is not None and not (isinstance(val, float) and math.isnan(val)):
                effects.append(abs(val))
        rec["best_abs_effect"] = max(effects) if effects else None
        rec["evidence_score"] = 2 * rec["flag_lm"] + rec["flag_partial"] + rec["flag_group_corr"]
        data.append(rec)
    df = pd.DataFrame(data)
    if not df.empty:
        df.sort_values(by=["evidence_score", "best_abs_effect"], ascending=[False, False], inplace=True)
    return df


def write_readme(path: str, info: Dict[str, Any], logger: logging.Logger):
    lines = []
    lines.append("Proteome 2D Analysis Summary")
    lines.append("=" * 40)
    lines.append("")
    lines.append("Input files used:")
    for item in info.get("inputs", []):
        lines.append(f"- {item}")
    lines.append("")
    lines.append(f"Total samples included: {info.get('n_samples', 'NA')}")
    lines.append("")
    lines.append("Matched proteins:")
    if info.get("matched_proteins"):
        for name in info.get("matched_proteins"):
            lines.append(f"- {name}")
    else:
        lines.append("- None")
    lines.append("")
    lines.append("Matched pathways:")
    if info.get("matched_pathways"):
        for name in info.get("matched_pathways"):
            lines.append(f"- {name}")
    else:
        lines.append("- None")
    lines.append("")
    lines.append("Significance coding:")
    lines.append("- ** : q < 0.05")
    lines.append("- *  : q < 0.10")
    lines.append("")
    lines.append("Top findings:")
    if info.get("top_findings"):
        for item in info.get("top_findings"):
            lines.append(f"- {item}")
    else:
        lines.append("- None")
    lines.append("")
    lines.append(f"Metadata inferred: {info.get('metadata_inferred', False)}")
    content = "\n".join(lines)
    try:
        with open(path, "w", encoding="utf-8-sig") as f:
            f.write(content)
        logger.info("Saved README %s", path)
    except Exception as exc:
        logger.error("Failed to write README %s: %s", path, exc)
        logger.error(traceback.format_exc())


def write_session_info(path: str, logger: logging.Logger):
    import platform
    import importlib

    lines = []
    lines.append("Python session information")
    lines.append("=" * 40)
    lines.append(f"Python version: {sys.version}")
    lines.append(f"Platform: {platform.platform()}")
    lines.append("")
    lines.append("Packages:")
    for pkg in ["pandas", "numpy", "scipy", "statsmodels", "pingouin", "seaborn", "matplotlib", "openpyxl", "xlsxwriter"]:
        try:
            module = importlib.import_module(pkg)
            version = getattr(module, "__version__", "unknown")
        except Exception:
            version = "not available"
        lines.append(f"- {pkg}: {version}")
    content = "\n".join(lines)
    try:
        with open(path, "w", encoding="utf-8-sig") as f:
            f.write(content)
        logger.info("Saved session info %s", path)
    except Exception as exc:
        logger.error("Failed to write session info %s: %s", path, exc)
        logger.error(traceback.format_exc())


def main():
    logger = setup_logging()
    logger.info("Starting proteome 2D analysis")
    ensure_packages(logger)
    try:
        import pandas as pd
    except Exception as exc:
        logger.error("Critical import failure: %s", exc)
        logger.error(traceback.format_exc())
        return

    configure_matplotlib(logger)

    scores_path = os.path.join(OUTDIR, "ssgsea_scores.csv")
    proteins_path = os.path.join(OUTDIR, "protein_matrix_wide.csv")
    metadata_path = os.path.join(OUTDIR, "sample_metadata.csv")

    try:
        ssgsea_raw = load_csv_safe(scores_path, logger)
        proteins_raw = load_csv_safe(proteins_path, logger)
    except Exception as exc:
        logger.error("Failed to load core matrices: %s", exc)
        return

    orientation_scores, orientation_proteins = determine_orientation(ssgsea_raw, proteins_raw, logger)
    ssgsea_df, scores_col_map, scores_idx_map = maybe_transpose(ssgsea_raw, "ssgsea", orientation_scores, logger)
    proteins_df, proteins_col_map, proteins_idx_map = maybe_transpose(proteins_raw, "proteins", orientation_proteins, logger)

    common_samples = sorted(set(ssgsea_df.columns) & set(proteins_df.columns))
    logger.info("Common samples between matrices: %d", len(common_samples))

    metadata_inferred = False
    if os.path.exists(metadata_path):
        try:
            metadata_raw = load_csv_safe(metadata_path, logger)
        except Exception:
            metadata_raw = None
    else:
        metadata_raw = None

    import pandas as pd

    if metadata_raw is None or metadata_raw.empty:
        metadata_inferred = True
        logger.info("Inferring metadata from sample names.")
        records = []
        for sample in common_samples:
            name_lower = sample.lower()
            if any(tag in name_lower for tag in ["ba", "case", "patient"]):
                group = "BA"
            elif any(tag in name_lower for tag in ["ctrl", "control", "healthy"]):
                group = "Ctrl"
            else:
                group = "Unknown"
            if any(tag in name_lower for tag in ["chol", "bile", "duct"]):
                tissue = "Cholangio"
            elif any(tag in name_lower for tag in ["liv", "hep"]):
                tissue = "Liver"
            else:
                tissue = "Unknown"
            records.append({"sample_id": sample, "group": group, "tissue": tissue})
        metadata = pd.DataFrame(records)
        inferred_path = os.path.join(RESULT_DIR, "sample_metadata_inferred.csv")
        write_csv(metadata, inferred_path, logger)
    else:
        metadata_raw.columns = [norm_name(c) for c in metadata_raw.columns]
        sample_col = None
        for col in metadata_raw.columns:
            if "sample" in col:
                sample_col = col
                break
        if sample_col is None:
            metadata_raw["sample_id"] = [norm_name(idx) for idx in metadata_raw.index]
        else:
            metadata_raw["sample_id"] = metadata_raw[sample_col].apply(norm_name)
        if "group" not in metadata_raw.columns:
            metadata_raw["group"] = "Unknown"
        if "tissue" not in metadata_raw.columns:
            metadata_raw["tissue"] = "Unknown"
        metadata = metadata_raw[["sample_id", "group", "tissue"]].copy()
        metadata["sample_id"] = metadata["sample_id"].apply(norm_name)

    metadata.index = metadata["sample_id"].apply(norm_name)
    metadata["group"] = metadata["group"].fillna("Unknown")
    metadata["tissue"] = metadata["tissue"].fillna("Unknown")
    metadata = metadata.loc[~metadata.index.duplicated(keep="first")]

    sample_intersection = sorted(set(common_samples) & set(metadata.index))
    removed_samples = sorted(set(ssgsea_df.columns) | set(proteins_df.columns) | set(metadata.index) - set(sample_intersection))
    if removed_samples:
        logger.info("Samples removed after alignment: %s", ", ".join(removed_samples))
    if not sample_intersection:
        logger.error("No overlapping samples after alignment. Exiting.")
        return

    ssgsea_aligned = ssgsea_df[sample_intersection].transpose()
    proteins_aligned = proteins_df[sample_intersection].transpose()
    metadata_aligned = metadata.loc[sample_intersection]

    write_csv(ssgsea_aligned, os.path.join(RESULT_DIR, "scores_mat_2D.csv"), logger)
    write_csv(proteins_aligned, os.path.join(RESULT_DIR, "proteins_mat_2D.csv"), logger)

    matched_pathways = fuzzy_map(KEY_PATHWAYS, list(ssgsea_aligned.columns), logger)
    matched_proteins = fuzzy_map(KEY_PROTEINS, list(proteins_aligned.columns), logger)

    matched_pathway_cols = {k: v for k, v in matched_pathways.items() if v is not None}
    matched_protein_cols = {k: v for k, v in matched_proteins.items() if v is not None}

    print(
        f"找到的样本数量: {len(sample_intersection)}，成功匹配的蛋白数量: {len(matched_protein_cols)}, 成功匹配的路径数量: {len(matched_pathway_cols)}"
    )

    try:
        import seaborn as sns
        cmap = sns.color_palette("coolwarm", as_cmap=True)
    except Exception:
        cmap = None

    group_corr_df = compute_spearman_groupwise(
        proteins_aligned,
        ssgsea_aligned,
        metadata_aligned,
        matched_proteins,
        matched_pathways,
        logger,
    )
    if not group_corr_df.empty:
        grouped = {group: sub for group, sub in group_corr_df.groupby("group")}
        save_excel(grouped, os.path.join(RESULT_DIR, "2D_groupwise_correlations.xlsx"), logger)
        for group, df in grouped.items():
            create_heatmap(
                df,
                "rho",
                f"Group Spearman ({group})",
                os.path.join(RESULT_DIR, f"Figure_2D_GroupCorr_Heatmap_{group}.png"),
                cmap,
                logger,
            )
    else:
        logger.warning("Group-wise correlation results are empty.")

    partial_corrs = compute_partial_correlations(
        proteins_aligned,
        ssgsea_aligned,
        metadata_aligned,
        matched_proteins,
        matched_pathways,
        logger,
    )
    if partial_corrs:
        save_excel(partial_corrs, os.path.join(RESULT_DIR, "2D_partial_correlations.xlsx"), logger)
        if "all_partial" in partial_corrs and not partial_corrs["all_partial"].empty:
            create_heatmap(
                partial_corrs["all_partial"],
                "rho",
                "Partial Spearman (All samples)",
                os.path.join(RESULT_DIR, "Figure_2D_PartialCorr_Heatmap.png"),
                cmap,
                logger,
            )
    else:
        logger.warning("Partial correlation results are empty.")

    lm_df = compute_linear_models(
        proteins_aligned,
        ssgsea_aligned,
        metadata_aligned,
        matched_proteins,
        matched_pathways,
        logger,
    )
    if not lm_df.empty:
        save_excel({"linear_models": lm_df}, os.path.join(RESULT_DIR, "2D_linear_models.xlsx"), logger)
        significant = lm_df[lm_df["q"].fillna(1) < 0.1].copy()
        if significant.empty:
            temp = lm_df.copy()
            temp["abs_beta"] = temp["beta_protein"].abs()
            top_pairs = temp.sort_values("abs_beta", ascending=False).head(6).to_dict("records")
        else:
            significant["abs_beta"] = significant["beta_protein"].abs()
            top_pairs = significant.sort_values("abs_beta", ascending=False).head(6).to_dict("records")
        create_scatter_plots(
            top_pairs,
            proteins_aligned,
            ssgsea_aligned,
            metadata_aligned,
            os.path.join(RESULT_DIR, "Figure_2D_LM_Top6_scatter.pdf"),
            logger,
        )
    else:
        logger.warning("Linear model dataframe is empty.")

    evidence_df = integrate_evidence(group_corr_df, partial_corrs, lm_df, logger)
    if evidence_df is not None and not evidence_df.empty:
        write_csv(evidence_df, os.path.join(RESULT_DIR, "evidence_rank.csv"), logger)
    else:
        logger.warning("Evidence dataframe is empty.")

    matched_protein_names = []
    for target, actual in matched_proteins.items():
        if actual is not None:
            original = proteins_idx_map.get(actual, actual)
            matched_protein_names.append(f"{target} -> {actual}")
        else:
            matched_protein_names.append(f"{target} -> None")
    matched_pathway_names = []
    for target, actual in matched_pathways.items():
        if actual is not None:
            original = scores_idx_map.get(actual, actual)
            matched_pathway_names.append(f"{target} -> {actual}")
        else:
            matched_pathway_names.append(f"{target} -> None")

    top_findings = []
    if evidence_df is not None and not evidence_df.empty:
        top_rows = evidence_df.head(10)
        for _, row in top_rows.iterrows():
            protein = row["protein"]
            pathway = row["pathway"]
            segments = [f"{protein} 与 {pathway}"]
            if row.get("best_group_rho") is not None:
                segments.append(f"组内ρ={row['best_group_rho']:.3f} (q={row.get('best_group_q', float('nan')):.3g})")
            if row.get("best_partial_rho") is not None:
                segments.append(f"偏相关ρ={row['best_partial_rho']:.3f} (q={row.get('partial_q', float('nan')):.3g})")
            if row.get("beta_lm") is not None:
                segments.append(f"线性模型β={row['beta_lm']:.3f} (q={row.get('lm_q', float('nan')):.3g})")
            top_findings.append("，".join(segments))

    readme_info = {
        "inputs": [
            f"ssgsea_scores.csv ({ssgsea_raw.shape[0]} x {ssgsea_raw.shape[1]})",
            f"protein_matrix_wide.csv ({proteins_raw.shape[0]} x {proteins_raw.shape[1]})",
            f"sample_metadata.csv (exists={os.path.exists(metadata_path)})",
        ],
        "n_samples": len(sample_intersection),
        "matched_proteins": matched_protein_names,
        "matched_pathways": matched_pathway_names,
        "top_findings": top_findings,
        "metadata_inferred": metadata_inferred,
    }

    write_readme(os.path.join(RESULT_DIR, "README_2D.txt"), readme_info, logger)
    write_session_info(os.path.join(RESULT_DIR, "sessionInfo_2D.txt"), logger)

    print(f"2D 分析完成，输出位于 {RESULT_DIR}")
    logger.info("Proteome 2D analysis finished.")


if __name__ == "__main__":
    main()
