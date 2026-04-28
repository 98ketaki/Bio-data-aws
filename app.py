"""Bouncer — Streamlit frontend."""

import os
import sys
import tempfile
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)

# ensure local bouncer package is importable
sys.path.insert(0, str(Path(__file__).parent))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import streamlit as st

import bouncer
from bouncer.qc.finding import Severity

# ── page config ──────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="Bouncer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

SCHEMA      = Path(__file__).parent / "schema/basic-schema.yaml"
QC_CONTRACT = Path(__file__).parent / "schema/basic-qc.yaml"
DB_PATH     = Path(__file__).parent / "bouncer_store.duckdb"

# ── helpers ───────────────────────────────────────────────────────────────────

SEV_COLOUR = {
    Severity.HARD:    ("#fce8e8", "#b22222", "✖ HARD"),
    Severity.SOFT:    ("#fff8e1", "#b8860b", "⚠ SOFT"),
    Severity.WARNING: ("#eaf4fb", "#1a6fa0", "ℹ ADVISORY"),
}

SECTION_LABELS = {
    "ct_export":     "CT Export Integrity",
    "ntc":           "NTC Contamination",
    "reference_genes": "Reference Genes",
    "samplesheet":   "Samplesheet Integrity",
    "cross_file":    "Cross-File Concordance",
    "design":        "Experimental Design",
}

sns.set_theme(style="whitegrid", font_scale=1.05)


def save_upload(uploaded) -> Path:
    suffix = Path(uploaded.name).suffix
    tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
    tmp.write(uploaded.read())
    tmp.flush()
    return Path(tmp.name)


def severity_badge(sev: Severity) -> str:
    bg, fg, label = SEV_COLOUR[sev]
    return (
        f'<span style="background:{bg};color:{fg};'
        f'font-weight:600;padding:2px 8px;border-radius:4px;'
        f'font-size:0.78rem;border:1px solid {fg}33">{label}</span>'
    )


def render_findings(findings):
    if not findings:
        st.success("All checks passed for this section.")
        return
    for f in findings:
        bg, fg, _ = SEV_COLOUR[f.severity]
        badge = severity_badge(f.severity)
        with st.container():
            st.markdown(
                f"""
                <div style="background:{bg};border-left:4px solid {fg};
                            border-radius:4px;padding:10px 14px;margin-bottom:8px">
                  <div style="display:flex;align-items:center;gap:8px;margin-bottom:4px">
                    {badge}
                    <code style="font-size:0.85rem">{f.check}</code>
                  </div>
                  <div style="font-size:0.9rem;color:#333">{f.message}</div>
                  {"<div style='font-size:0.8rem;color:#555;margin-top:4px'>Wells: " + ", ".join(str(w) for w in f.wells[:20]) + ("…" if len(f.wells) > 20 else "") + "</div>" if f.wells else ""}
                  {"<div style='font-size:0.8rem;color:#555;margin-top:2px'>Samples: " + ", ".join(str(s) for s in f.samples[:10]) + ("…" if len(f.samples) > 10 else "") + "</div>" if f.samples else ""}
                </div>
                """,
                unsafe_allow_html=True,
            )


# ── sidebar ───────────────────────────────────────────────────────────────────

with st.sidebar:
    st.image("https://em-content.zobj.net/source/apple/354/dna_1f9ec.png", width=48)
    st.title("Bouncer")
    st.caption(f"v{bouncer.__version__}  ·  qPCR QC & Feature Store")
    st.divider()

    page = st.radio(
        "Navigate",
        ["QC Check", "Feature Store", "Visualise"],
        label_visibility="collapsed",
    )
    st.divider()

    st.subheader("Upload files")
    ct_file  = st.file_uploader("Ct export CSV",   type=["csv"], key="ct")
    ss_file  = st.file_uploader("Samplesheet CSV", type=["csv"], key="ss")
    pdf_file = st.file_uploader("Protocol PDF (optional)", type=["pdf"], key="pdf")

    use_demo = st.checkbox("Use demo dataset (correct)", value=True)
    st.caption("Schema and QC contract are loaded from `schema/`.")

# ── resolve paths ─────────────────────────────────────────────────────────────

def resolve_paths():
    """Return (ct_path, ss_path, pdf_path) or None if nothing available."""
    if ct_file and ss_file:
        ct_path  = save_upload(ct_file)
        ss_path  = save_upload(ss_file)
        pdf_path = save_upload(pdf_file) if pdf_file else None
        return ct_path, ss_path, pdf_path
    if use_demo:
        base = Path(__file__).parent / "qPCR" / "correct"
        return (
            base / "Ct.csv",
            base / "samplesheet-qpcr - correct.csv",
            base / "protocol.pdf",
        )
    return None


# ── PAGE: QC Check ────────────────────────────────────────────────────────────

if page == "QC Check":
    st.title("QC Check")
    st.markdown(
        "Upload your files in the sidebar (or use the demo dataset) then click **Run QC**."
    )

    paths = resolve_paths()
    if paths is None:
        st.info("Upload a Ct export CSV and samplesheet CSV to get started.")
        st.stop()

    ct_path, ss_path, pdf_path = paths

    if st.button("▶ Run QC", type="primary", use_container_width=True):
        with st.spinner("Running checks…"):
            result = bouncer.check(
                ct_export=ct_path,
                samplesheet=ss_path,
                schema=SCHEMA,
                qc_contract=QC_CONTRACT,
                protocol=pdf_path,
                verbose=False,
            )
        st.session_state["result"] = result

    result = st.session_state.get("result")
    if result is None:
        st.stop()

    # ── status banner
    col1, col2, col3, col4 = st.columns(4)
    status_colour = "green" if result.passed else "red"
    status_label  = "PASS ✓" if result.passed else "FAIL ✗"
    col1.metric("Status", status_label)
    col2.metric("Hard failures", len(result.hard_failures))
    col3.metric("Soft failures", len(result.soft_failures))
    col4.metric("Advisories",    len(result.warnings))

    if result.ct:
        st.markdown(
            f"**Experiment** `{result.ct.experiment_name}` &nbsp;·&nbsp; "
            f"**Plate** `{result.ct.plate_barcode}` &nbsp;·&nbsp; "
            f"**Run date** `{result.ct.run_date}` &nbsp;·&nbsp; "
            f"**Operator** `{result.ct.operator}`"
        )

    st.divider()

    # ── per-section findings
    sections_with_findings: dict[str, list] = {}
    for f in result.findings:
        sections_with_findings.setdefault(f.section, []).append(f)

    for key, label in SECTION_LABELS.items():
        section_findings = sections_with_findings.get(key, [])
        n_hard = sum(1 for f in section_findings if f.severity == Severity.HARD)
        n_soft = sum(1 for f in section_findings if f.severity == Severity.SOFT)
        n_warn = sum(1 for f in section_findings if f.severity == Severity.WARNING)

        icon = "✅" if not section_findings else ("❌" if n_hard else ("⚠️" if n_soft else "ℹ️"))
        with st.expander(f"{icon} {label}  ({n_hard} hard · {n_soft} soft · {n_warn} advisory)",
                         expanded=bool(n_hard or n_soft)):
            render_findings(section_findings)

    st.divider()

    # ── findings table
    with st.expander("All findings as table"):
        rows = [
            dict(severity=f.severity.value, section=f.section,
                 check=f.check, message=f.message)
            for f in result.findings
        ]
        if rows:
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
        else:
            st.info("No findings.")

    # ── register button
    st.divider()
    if result.passed:
        st.success("All hard checks passed — data is ready to register.")
        if st.button("📥 Register to feature store", type="primary"):
            with st.spinner("Building AnnData and writing to DuckDB…"):
                adata = bouncer.register(result, db_path=DB_PATH)
            st.session_state["adata"] = adata
            st.success(
                f"Registered **{len(adata.obs_names)} samples × {len(adata.var_names)} targets** "
                f"to the feature store."
            )
    else:
        st.error(
            f"{len(result.hard_failures)} hard failure(s) block registration. "
            "Fix the issues highlighted above and re-run."
        )


# ── PAGE: Feature Store ───────────────────────────────────────────────────────

elif page == "Feature Store":
    st.title("Feature Store")

    # ── auto-populate store with demo data if empty
    if not DB_PATH.exists():
        paths = resolve_paths()
        if paths:
            ct_path, ss_path, pdf_path = paths
            result = bouncer.check(
                ct_export=ct_path, samplesheet=ss_path,
                schema=SCHEMA, qc_contract=QC_CONTRACT,
                protocol=pdf_path, verbose=False,
            )
            if result.passed:
                bouncer.register(result, db_path=DB_PATH)

    experiments = bouncer.list_experiments(db_path=DB_PATH)
    if experiments.empty:
        st.info("No experiments registered yet. Run QC and register from the **QC Check** page.")
        st.stop()

    st.subheader("Registered experiments")
    st.dataframe(experiments, use_container_width=True, hide_index=True)

    st.divider()
    st.subheader("Pull data")

    # ── filter widgets
    all_df = bouncer.pull_data(db_path=DB_PATH)

    col1, col2, col3 = st.columns(3)
    with col1:
        assay_opts = sorted(all_df["assay_type"].dropna().unique())
        sel_assay = st.multiselect("Assay type", assay_opts, default=assay_opts)
    with col2:
        cond_opts = sorted(all_df["condition"].dropna().unique())
        sel_cond  = st.multiselect("Condition", cond_opts, default=cond_opts)
    with col3:
        ttype_opts = sorted(all_df["target_type"].dropna().unique())
        sel_ttype  = st.multiselect("Target type", ttype_opts, default=ttype_opts)

    col4, col5 = st.columns(2)
    with col4:
        tgt_opts = sorted(all_df["target_name"].dropna().unique())
        sel_tgt  = st.multiselect("Target name", tgt_opts, default=tgt_opts)
    with col5:
        org_opts = sorted(all_df["organism"].dropna().unique())
        sel_org  = st.multiselect("Organism", org_opts, default=org_opts)

    df = bouncer.pull_data(
        db_path=DB_PATH,
        assay=sel_assay      or None,
        condition=sel_cond   or None,
        target_type=sel_ttype or None,
        target_name=sel_tgt  or None,
        organism=sel_org     or None,
    )

    st.caption(f"{len(df)} rows matched")

    show_cols = ["sample_id", "target_name", "target_type", "ct_mean", "ct_sd",
                 "condition", "organism", "sex", "batch", "experiment_id"]
    st.dataframe(
        df[[c for c in show_cols if c in df.columns]],
        use_container_width=True, hide_index=True,
    )

    csv = df.to_csv(index=False).encode()
    st.download_button("⬇ Download CSV", csv, "bouncer_query.csv", "text/csv")


# ── PAGE: Visualise ───────────────────────────────────────────────────────────

elif page == "Visualise":
    st.title("Visualise")

    # ── load data
    if not DB_PATH.exists():
        paths = resolve_paths()
        if paths:
            ct_path, ss_path, pdf_path = paths
            result = bouncer.check(
                ct_export=ct_path, samplesheet=ss_path,
                schema=SCHEMA, qc_contract=QC_CONTRACT,
                protocol=pdf_path, verbose=False,
            )
            if result.passed:
                bouncer.register(result, db_path=DB_PATH)

    df_all = bouncer.pull_data(db_path=DB_PATH)
    if df_all.empty:
        st.info("No data in the feature store yet. Register an experiment first.")
        st.stop()

    adata = bouncer.pull_data(db_path=DB_PATH, as_anndata=True)

    tab1, tab2, tab3, tab4 = st.tabs([
        "CT heatmap", "Replicate SD", "Condition comparison", "Reference stability"
    ])

    # ── CT heatmap ────────────────────────────────────────────────────────────
    with tab1:
        st.subheader("Mean Ct — sample × target")
        ct_df = pd.DataFrame(
            adata.X,
            index=[f"{s} ({c})" for s, c in zip(adata.obs_names, adata.obs["condition"])],
            columns=adata.var_names,
        )
        fig, ax = plt.subplots(figsize=(8, max(3, len(ct_df) * 0.7 + 1)))
        sns.heatmap(
            ct_df, annot=True, fmt=".1f", cmap="YlOrRd_r",
            linewidths=0.5, cbar_kws={"label": "Mean Ct"}, ax=ax,
        )
        ax.set_xlabel("Target"); ax.set_ylabel("Sample")
        plt.tight_layout()
        st.pyplot(fig, use_container_width=True)
        plt.close(fig)

    # ── Replicate SD ──────────────────────────────────────────────────────────
    with tab2:
        st.subheader("Technical replicate SD")
        st.caption("Soft threshold 0.5 · Hard threshold 1.0")
        sd_df = pd.DataFrame(
            adata.layers["ct_sd"],
            index=[f"{s} ({c})" for s, c in zip(adata.obs_names, adata.obs["condition"])],
            columns=adata.var_names,
        )
        fig, ax = plt.subplots(figsize=(8, max(3, len(sd_df) * 0.7 + 1)))
        sns.heatmap(
            sd_df, annot=True, fmt=".3f", cmap="Blues",
            linewidths=0.5, vmin=0, vmax=0.5,
            cbar_kws={"label": "CT SD"}, ax=ax,
        )
        ax.set_xlabel("Target"); ax.set_ylabel("Sample")
        plt.tight_layout()
        st.pyplot(fig, use_container_width=True)
        plt.close(fig)

    # ── Condition comparison ──────────────────────────────────────────────────
    with tab3:
        st.subheader("Ct by condition — target genes")
        df_tgt = df_all[df_all["target_type"] == "target"].copy()
        if df_tgt.empty:
            st.info("No target-type genes in the store.")
        else:
            palette = dict(zip(
                sorted(df_tgt["condition"].unique()),
                sns.color_palette("tab10", n_colors=df_tgt["condition"].nunique()),
            ))
            fig, ax = plt.subplots(figsize=(8, 4))
            sns.boxplot(data=df_tgt, x="target_name", y="ct_mean",
                        hue="condition", palette=palette,
                        width=0.5, linewidth=1.2, fliersize=0, ax=ax)
            sns.stripplot(data=df_tgt, x="target_name", y="ct_mean",
                          hue="condition", palette=palette,
                          dodge=True, size=7, linewidth=0.8, alpha=0.85, ax=ax)
            handles, labels = ax.get_legend_handles_labels()
            n = df_tgt["condition"].nunique()
            ax.legend(handles[:n], labels[:n], title="Condition")
            ax.set_xlabel("Target gene"); ax.set_ylabel("Mean Ct")
            ax.invert_yaxis()
            ax.set_title("Lower Ct = higher expression")
            plt.tight_layout()
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)

    # ── Reference stability ───────────────────────────────────────────────────
    with tab4:
        st.subheader("Reference gene stability (CV across samples)")
        df_ref = df_all[df_all["target_type"] == "reference"].copy()
        if df_ref.empty:
            st.info("No reference-type genes in the store.")
        else:
            cv_data = (
                df_ref.groupby("target_name")["ct_mean"]
                .agg(mean="mean", std="std")
                .assign(cv=lambda d: d["std"] / d["mean"] * 100)
                .reset_index()
            )

            colours = ["#2ca02c" if v < 15 else "#d62728" for v in cv_data["cv"]]
            fig, ax = plt.subplots(figsize=(6, 3.5))
            ax.bar(cv_data["target_name"], cv_data["cv"], color=colours,
                   edgecolor="white", linewidth=0.8)
            ax.axhline(15, color="orange", linestyle="--", linewidth=1.2,
                       label="Soft (15 %)")
            ax.axhline(25, color="red", linestyle="--", linewidth=1.2,
                       label="Hard (25 %)")
            ax.set_ylabel("CV (%)"); ax.set_xlabel("Reference gene")
            ax.legend(fontsize=9)
            plt.tight_layout()
            st.pyplot(fig, use_container_width=True)
            plt.close(fig)

            st.dataframe(
                cv_data.rename(columns={"mean": "mean_ct", "std": "std_ct",
                                        "cv": "cv_%"})
                       .round(3),
                use_container_width=True, hide_index=True,
            )
