import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os


st.set_page_config(page_title="Phenotypic Screen", page_icon="🔭", layout="wide")

TEAL   = "#0D869B"
ORANGE = "#E66F02"
GREEN  = "#3fb950"
PURPLE = "#9b5de5"

st.markdown(f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@300;400;500;600&family=DM+Serif+Display&display=swap');
html, body, [class*="css"] {{ font-family: 'DM Sans', sans-serif; }}
.section-title {{
    font-family: 'DM Serif Display', serif;
    font-size: 1.4rem;
    margin-bottom: 0.5rem;
    margin-top: 1.5rem;
    color: {TEAL};
}}
.info-card {{
    border-radius: 12px;
    padding: 20px 24px;
    border: 1px solid #88888833;
    margin-bottom: 12px;
}}
.info-row {{
    display: flex;
    justify-content: space-between;
    padding: 8px 0;
    border-bottom: 1px solid #88888822;
    font-size: 0.9rem;
}}
.info-label {{ font-weight: 500; opacity: 0.6; min-width: 220px; }}
.info-value {{ font-weight: 600; }}
.badge {{
    display: inline-block;
    border-radius: 6px;
    padding: 3px 10px;
    font-size: 0.78rem;
    font-weight: 600;
    margin: 2px 3px;
}}
.badge-teal   {{ background: {TEAL}22;   color: {TEAL};   border: 1px solid {TEAL}55; }}
.badge-orange {{ background: {ORANGE}22; color: {ORANGE}; border: 1px solid {ORANGE}55; }}
.badge-green  {{ background: {GREEN}22;  color: {GREEN};  border: 1px solid {GREEN}55; }}
.badge-purple {{ background: {PURPLE}22; color: {PURPLE}; border: 1px solid {PURPLE}55; }}
.badge-grey   {{ background: #88888822;  color: #888;     border: 1px solid #88888844; }}
</style>
""", unsafe_allow_html=True)

# ── Load data ─────────────────────────────────────────────────────────────────
@st.cache_data
def load_data():
    inf_raw = pd.read_csv("data/infezione.csv", sep=None, engine="python")
    inf_raw.columns = inf_raw.columns.str.strip()
    inf_raw["Gene"] = inf_raw["WELL NAME"].astype(str).str.extract(r"-\s*(.+)$")[0].str.strip().str.upper()

    num_cols = ["INFECTION 1","INFECTION 2","INFECTION 3","INFECTION 4",
                "REINFECTION 1","REINFECTION 2","REINFECTION 3",
                "REINFECTION 4","REINFECTION 5","REINFECTION 6"]
    for col in num_cols:
        inf_raw[col] = inf_raw[col].astype(str).str.replace(",",".").str.strip()
        inf_raw[col] = pd.to_numeric(inf_raw[col], errors="coerce")

    inf_raw["Mean_Infection"]   = inf_raw[["INFECTION 1","INFECTION 2","INFECTION 3","INFECTION 4"]].mean(axis=1)
    inf_raw["Mean_Reinfection"] = inf_raw[["REINFECTION 1","REINFECTION 2","REINFECTION 3","REINFECTION 4","REINFECTION 5","REINFECTION 6"]].mean(axis=1)
    inf_raw["Std_Infection"]    = inf_raw[["INFECTION 1","INFECTION 2","INFECTION 3","INFECTION 4"]].std(axis=1)
    inf_raw["Std_Reinfection"]  = inf_raw[["REINFECTION 1","REINFECTION 2","REINFECTION 3","REINFECTION 4","REINFECTION 5","REINFECTION 6"]].std(axis=1)

    ctrl_inf    = inf_raw[inf_raw["Gene"] == "INF"]["Mean_Infection"].mean()
    ctrl_mock   = inf_raw[inf_raw["Gene"] == "MOCK"]["Mean_Infection"].mean()
    ctrl_inf_r  = inf_raw[inf_raw["Gene"] == "INF"]["Mean_Reinfection"].mean()
    ctrl_mock_r = inf_raw[inf_raw["Gene"] == "MOCK"]["Mean_Reinfection"].mean()

    controls = ["MOCK", "INF", "PLK", "VMP1", "EMPTY"]
    gene_df = inf_raw[~inf_raw["Gene"].isin(controls)].groupby("Gene").agg(
        Mean_Infection   = ("Mean_Infection",   "mean"),
        Mean_Reinfection = ("Mean_Reinfection", "mean"),
        Std_Infection    = ("Std_Infection",    "mean"),
        Std_Reinfection  = ("Std_Reinfection",  "mean"),
    ).reset_index()

    gene_df["Infection_norm"]   = gene_df["Mean_Infection"]   / ctrl_inf
    gene_df["Reinfection_norm"] = gene_df["Mean_Reinfection"] / ctrl_inf_r

    crit = pd.read_csv("data/criteri.csv", sep=None, engine="python")
    crit.columns = crit.columns.str.strip()
    col_map = {
        crit.columns[0]: "Known Host Factor",
        crit.columns[1]: "Interactor of Known HF",
        crit.columns[2]: "Phosphosite SARS-CoV-2",
        crit.columns[3]: "Glycosite SARS-CoV-2",
        crit.columns[4]: "Enzymatic Activity",
        crit.columns[5]: "Druggable",
        crit.columns[6]: "New Interactions Post-Translocation",
    }
    crit = crit.rename(columns=col_map)
    criteria_cols = list(col_map.values())

    gene_criteria = {}
    for col in criteria_cols:
        for gene in crit[col].dropna():
            g = str(gene).strip().upper()
            if g not in gene_criteria:
                gene_criteria[g] = []
            gene_criteria[g].append(col)

    for col in criteria_cols:
        gene_df[col] = gene_df["Gene"].map(lambda g, c=col: c in gene_criteria.get(g, []))

    return gene_df, criteria_cols, ctrl_inf, ctrl_mock, ctrl_inf_r, ctrl_mock_r, inf_raw


@st.cache_data
def load_arbo():
    arbo = pd.read_csv("data/arbo.csv", sep=None, engine="python")
    arbo.columns = arbo.columns.str.strip()
    arbo = arbo.rename(columns={"T: GENE NAME": "Gene"})
    arbo["Gene"] = arbo["Gene"].str.strip().str.upper()

    viruses = {
        "RVFV": ("Rift Valley Fever", ["RVFV 1","RVFV 2","RVFV 3","RVFV 4"]),
        "WNV":  ("West Nile Virus",   ["WNV 1 ","WNV 2","WNV 3","WNV 4"]),
        "MAYV": ("Mayaro Virus",      ["MAYV 1 ","MAYV 2","MAYV 3","MAYV 4"]),
        "ZIKV": ("Zika Virus",        ["ZIKV 1 ","ZIKV 2","ZIKV 3","ZIKV 4"]),
        "SFV":  ("Semliki Forest",    ["SFV 1","SFV 2","SFV 3","SFV 4"]),
    }

    for code, (name, cols) in viruses.items():
        for c in cols:
            if c in arbo.columns:
                arbo[c] = arbo[c].astype(str).str.replace(",",".").str.strip()
                arbo[c] = pd.to_numeric(arbo[c], errors="coerce")
        valid = [c for c in cols if c in arbo.columns]
        arbo[f"{code}_mean"] = arbo[valid].mean(axis=1)
        arbo[f"{code}_std"]  = arbo[valid].std(axis=1)

    return arbo, viruses


gene_df, criteria_cols, ctrl_inf, ctrl_mock, ctrl_inf_r, ctrl_mock_r, inf_raw = load_data()
arbo_df, VIRUSES = load_arbo()

badge_colors = ["teal","orange","green","purple","teal","orange","green"]
VIRUS_COLORS = {
    "RVFV": "#f85149",
    "WNV":  "#58a6ff",
    "MAYV": "#E66F02",
    "ZIKV": "#3fb950",
    "SFV":  "#9b5de5",
}

# ── Header ────────────────────────────────────────────────────────────────────
st.title("🔭 Phenotypic Screen")
st.markdown("""
From the **773 proteins** that relocalize upon HCoV-OC43 infection, we selected
**166 candidates** based on functional criteria and performed a targeted siRNA screen
using SARS-CoV-2 infection as readout (High-Content Screening).
Infection was quantified by automated imaging and normalized to the infected control.
A subset of candidates was also tested against a panel of **arboviruses** to assess broad-spectrum antiviral activity.
""")
st.divider()

# ── Stats ─────────────────────────────────────────────────────────────────────
crit_counts = {col: int(gene_df[col].sum()) for col in criteria_cols}
c1, c2, c3, c4 = st.columns(4)
c1.metric("Candidates screened", len(gene_df))
c2.metric("Druggable targets", crit_counts["Druggable"])
c3.metric("New interactors post-infection", crit_counts["New Interactions Post-Translocation"])
c4.metric("Tested against arboviruses", len(arbo_df))

badges = "".join([f'<span class="badge badge-{color}">{col} ({crit_counts[col]})</span>'
                  for col, color in zip(criteria_cols, badge_colors)])
st.markdown(f"**Selection criteria:** {badges}", unsafe_allow_html=True)
st.divider()

# ── Find your protein ─────────────────────────────────────────────────────────
st.markdown('<div class="section-title">🔍 Find your protein</div>', unsafe_allow_html=True)

search = st.text_input(
    "Search by gene name (e.g. ACE2, NAPA, RND3...)",
    placeholder="Type a gene name..."
).strip().upper()

readout = st.radio("SARS-CoV-2 readout", ["Infection", "Reinfection"], horizontal=True)

if search:
    match = gene_df[gene_df["Gene"] == search]

    if match.empty:
        partial = gene_df[gene_df["Gene"].str.contains(search, na=False)]
        if not partial.empty:
            st.info("Exact match not found. Did you mean:")
            st.dataframe(partial[["Gene"] + criteria_cols].head(10), use_container_width=True, hide_index=True)
        else:
            st.warning(f"Gene **{search}** not found in the screen dataset.")
    else:
        row = match.iloc[0]

        # Check HCS images
        hcs_images = []
        hcs_folder = f"data/images_screen/{row['Gene']}"
        if os.path.isdir(hcs_folder):
            for f in sorted(os.listdir(hcs_folder)):
                if f.lower().endswith((".png",".jpg",".jpeg",".tif",".bmp")):
                    hcs_images.append(os.path.join(hcs_folder, f))
        else:
            for ext in [".png",".jpg",".jpeg",".tif",".bmp"]:
                p = f"data/images_screen/{row['Gene']}{ext}"
                if os.path.exists(p):
                    hcs_images = [p]
                    break
        has_hcs = len(hcs_images) > 0

        arbo_row = arbo_df[arbo_df["Gene"] == search]
        has_arbo = not arbo_row.empty

        col_info, col_right = st.columns([3, 2] if has_hcs else [1, 0.001])

        with col_info:
            gene_badges = "".join([
                f'<span class="badge badge-{color}">{col}</span>'
                for col, color in zip(criteria_cols, badge_colors) if row[col]
            ])
            st.markdown(f"""
            <div style="margin-bottom: 12px;">
                <span style="font-family: DM Serif Display, serif; font-size: 1.6rem; font-weight: bold;">{row['Gene']}</span>
            </div>
            {gene_badges or '<span style="opacity:0.5; font-size:0.85rem;">No specific criteria</span>'}
            """, unsafe_allow_html=True)
            st.markdown("")

            if readout == "Infection":
                val, norm, std = row["Mean_Infection"], row["Infection_norm"], row["Std_Infection"]
                ctrl, ctrl_m = ctrl_inf, ctrl_mock
                rep_cols = ["INFECTION 1","INFECTION 2","INFECTION 3","INFECTION 4"]
            else:
                val, norm, std = row["Mean_Reinfection"], row["Reinfection_norm"], row["Std_Reinfection"]
                ctrl, ctrl_m = ctrl_inf_r, ctrl_mock_r
                rep_cols = ["REINFECTION 1","REINFECTION 2","REINFECTION 3","REINFECTION 4","REINFECTION 5","REINFECTION 6"]

            if norm < 0.5:
                effect, effect_color = "🔵 Dependency factor", TEAL
            elif norm > 1.5:
                effect, effect_color = "🔴 Restriction factor", ORANGE
            else:
                effect, effect_color = "⚪ No significant effect", "#888"

            st.markdown(f"""
            <div style="margin: 12px 0; padding: 12px 16px; border-radius: 10px;
                        background: {effect_color}15; border: 1px solid {effect_color}44;">
                <span style="font-size: 1rem; font-weight: 600; color: {effect_color};">{effect}</span>
            </div>
            """, unsafe_allow_html=True)

            st.markdown('<div class="info-card">', unsafe_allow_html=True)
            for label_text, value in [
                (f"Mean {readout} (raw)", f"{val:.3f}"),
                ("Normalized to infected ctrl", f"{norm:.3f}"),
                ("Std deviation", f"{std:.3f}"),
                ("Infected control (ref = 1.0)", f"{ctrl:.3f}"),
                ("Mock control", f"{ctrl_m:.3f}"),
            ]:
                st.markdown(
                    f'<div class="info-row"><span class="info-label">{label_text}</span>'
                    f'<span class="info-value">{value}</span></div>',
                    unsafe_allow_html=True
                )
            st.markdown('</div>', unsafe_allow_html=True)

            # SARS-CoV-2 replicates bar chart
            gene_rows = inf_raw[inf_raw["Gene"] == search]
            rep_values = gene_rows[rep_cols].values.flatten()
            rep_values = rep_values[~np.isnan(rep_values)]

            fig_bar = go.Figure()
            fig_bar.add_trace(go.Bar(
                x=[f"Rep {i+1}" for i in range(len(rep_values))],
                y=rep_values,
                marker_color=TEAL, marker_line=dict(width=0),
                hovertemplate="%{y:.3f}<extra></extra>",
            ))
            fig_bar.add_hline(y=ctrl, line_dash="dash", line_color=ORANGE,
                              annotation_text="Infected ctrl", annotation_position="right")
            fig_bar.add_hline(y=ctrl_m, line_dash="dash", line_color="#888",
                              annotation_text="Mock ctrl", annotation_position="right")
            fig_bar.update_layout(
                title=dict(text=f"SARS-CoV-2 {readout} — {row['Gene']}", font=dict(size=12)),
                plot_bgcolor="#0d1117", paper_bgcolor="#0d1117",
                font=dict(color="#e6edf3"), height=260,
                margin=dict(l=40, r=80, t=40, b=30),
                yaxis=dict(title="Normalized infection", gridcolor="#21262d", zerolinecolor="#484f58"),
                xaxis=dict(gridcolor="#21262d"), showlegend=False,
            )
            st.plotly_chart(fig_bar, use_container_width=True)

            # Arbovirus bar chart per protein
            st.markdown('<div class="section-title">🌍 Arbovirus screen</div>', unsafe_allow_html=True)
            if has_arbo:
                st.caption("Normalized infection upon siRNA knockdown across a panel of arboviruses (4 replicates each).")
                ar = arbo_row.iloc[0]
                virus_names = [VIRUSES[v][0] for v in VIRUSES]
                virus_means = [ar[f"{v}_mean"] for v in VIRUSES]
                virus_stds  = [ar[f"{v}_std"]  for v in VIRUSES]
                arbo_colors = [VIRUS_COLORS[v] for v in VIRUSES]

                fig_arbo = go.Figure()
                fig_arbo.add_trace(go.Bar(
                    x=virus_names, y=virus_means,
                    error_y=dict(type="data", array=virus_stds, visible=True,
                                 color="rgba(255,255,255,0.35)"),
                    marker_color=arbo_colors, marker_line=dict(width=0),
                    hovertemplate="<b>%{x}</b><br>Mean: %{y:.3f}<extra></extra>",
                ))
                fig_arbo.add_hline(y=1.0, line_dash="dash", line_color=ORANGE,
                                   annotation_text="Infected ctrl", annotation_position="right")
                fig_arbo.add_hline(y=0.5, line_dash="dot", line_color=TEAL,
                                   annotation_text="0.5×", annotation_position="right")
                fig_arbo.update_layout(
                    plot_bgcolor="#0d1117", paper_bgcolor="#0d1117",
                    font=dict(color="#e6edf3"), height=280,
                    margin=dict(l=40, r=80, t=20, b=40),
                    yaxis=dict(title="Normalized infection", gridcolor="#21262d", zerolinecolor="#484f58"),
                    xaxis=dict(gridcolor="#21262d"), showlegend=False,
                )
                st.plotly_chart(fig_arbo, use_container_width=True)
            else:
                st.caption("This protein was not included in the arbovirus screen subset.")

        # Right column: HCS channel viewer
        with col_right:
            if has_hcs:
                st.markdown('<div class="section-title">🔬 HCS Immunofluorescence</div>', unsafe_allow_html=True)
                st.caption(f"{row['Gene']} · siRNA knockdown · SARS-CoV-2 infection")

                gene_folder = f"data/images_screen/{row['Gene']}"
                CHANNELS = {
                    "merge": ("Merge",    "#ffffff"),
                    "DAPI":  ("DAPI",     "#4488ff"),
                    "N":     ("N (virus)","#ff4444"),
                }

                available = {}
                for key, (label, color) in CHANNELS.items():
                    p = os.path.join(gene_folder, f"{row['Gene']}_{key}.bmp")
                    if not os.path.exists(p):
                        for ext in [".png",".jpg",".jpeg",".tif",".bmp"]:
                            p2 = os.path.join(gene_folder, f"{row['Gene']}_{key}{ext}")
                            if os.path.exists(p2):
                                p = p2
                                break
                    if os.path.exists(p):
                        available[key] = (label, color, p)

                if available:
                    st.markdown("**Select channels:**")
                    selected = {}
                    cols_ch = st.columns(len(available))
                    for i, (key, (label, color, path)) in enumerate(available.items()):
                        with cols_ch[i]:
                            default = key == "merge"
                            selected[key] = st.checkbox(label, value=default, key=f"ch_{row['Gene']}_{key}")

                    active = {k: v for k, v in available.items() if selected.get(k, False)}

                    if not active:
                        st.caption("Select at least one channel.")
                    elif len(active) == 1:
                        key, (label, color, path) = list(active.items())[0]
                        st.image(path, caption=label)
                    else:
                        from PIL import Image as PILImage
                        imgs = []
                        for key, (label, color, path) in active.items():
                            try:
                                img = PILImage.open(path).convert("RGB")
                                imgs.append(np.array(img).astype(np.float32))
                            except Exception:
                                pass

                        if imgs:
                            h = min(i.shape[0] for i in imgs)
                            w = min(i.shape[1] for i in imgs)
                            imgs = [i[:h,:w] for i in imgs]
                            result = np.ones((h, w, 3), dtype=np.float32)
                            for img_arr in imgs:
                                result *= (1.0 - img_arr / 255.0)
                            blended = ((1.0 - result) * 255).clip(0, 255).astype(np.uint8)
                            blended_img = PILImage.fromarray(blended)
                            labels = " + ".join([v[0] for v in active.values()])
                            st.image(blended_img, caption=labels)
                else:
                    for img_path in hcs_images:
                        label = os.path.splitext(os.path.basename(img_path))[0]
                        st.image(img_path, caption=label)

else:
    # ── Overview barplot ──────────────────────────────────────────────────────
    st.markdown('<div class="section-title">📊 Overview — all screened proteins</div>', unsafe_allow_html=True)

    y_col   = "Infection_norm" if readout == "Infection" else "Reinfection_norm"
    y_label = "Normalized Infection" if readout == "Infection" else "Normalized Reinfection"
    plot_df = gene_df.sort_values(y_col)
    colors  = [TEAL if v < 0.5 else ORANGE if v > 1.5 else "#484f58" for v in plot_df[y_col]]

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=plot_df["Gene"], y=plot_df[y_col],
        marker_color=colors, marker_line=dict(width=0),
        hovertemplate="<b>%{x}</b><br>Normalized: %{y:.3f}<extra></extra>",
    ))
    fig.add_hline(y=1.0, line_dash="dash", line_color=ORANGE,
                  annotation_text="Infected ctrl", annotation_position="right")
    fig.add_hline(y=0.5, line_dash="dot", line_color=TEAL, annotation_text="0.5×", annotation_position="right")
    fig.add_hline(y=1.5, line_dash="dot", line_color=ORANGE, annotation_text="1.5×", annotation_position="right")
    fig.update_layout(
        xaxis_title="Gene", yaxis_title=y_label,
        plot_bgcolor="#0d1117", paper_bgcolor="#0d1117",
        font=dict(color="#e6edf3"), height=500,
        margin=dict(l=60, r=80, t=30, b=80),
        xaxis=dict(gridcolor="#21262d", tickangle=-45, tickfont=dict(size=8)),
        yaxis=dict(gridcolor="#21262d", zerolinecolor="#484f58"),
        showlegend=False,
    )
    st.plotly_chart(fig, use_container_width=True)

    st.markdown(f"""
    <span class="badge badge-teal">🔵 Dependency factor (&lt;0.5)</span>
    <span class="badge badge-orange">🔴 Restriction factor (&gt;1.5)</span>
    <span class="badge badge-grey">⚪ No significant effect</span>
    """, unsafe_allow_html=True)

    with st.expander("📋 Show full table"):
        show_cols = ["Gene","Mean_Infection","Infection_norm","Mean_Reinfection","Reinfection_norm"] + criteria_cols
        st.dataframe(plot_df[[c for c in show_cols if c in gene_df.columns]].reset_index(drop=True),
                     use_container_width=True, hide_index=True)

    # ── Arbovirus heatmap ─────────────────────────────────────────────────────
    st.divider()
    st.markdown('<div class="section-title">🌍 Arbovirus screen overview</div>', unsafe_allow_html=True)
    st.caption(f"{len(arbo_df)} proteins tested across 5 arboviruses · values normalized to infected control")

    # Build matrix: genes x viruses (+ SARS-CoV-2)
    heatmap_df = arbo_df[["Gene"]].copy()
    sars_map = gene_df.set_index("Gene")["Infection_norm"].to_dict()
    heatmap_df["SARS-CoV-2"] = heatmap_df["Gene"].map(sars_map)
    for code, (name, _) in VIRUSES.items():
        heatmap_df[name] = arbo_df[f"{code}_mean"].values
    heatmap_df = heatmap_df.set_index("Gene")

    # Sort by mean across all viruses
    heatmap_df = heatmap_df.loc[heatmap_df.mean(axis=1).sort_values().index]

    genes          = heatmap_df.index.tolist()
    virus_cols     = heatmap_df.columns.tolist()
    z              = heatmap_df.values.tolist()

    fig_heat = go.Figure(go.Heatmap(
        z=z,
        x=virus_cols,
        y=genes,
        colorscale=[
            [0.0,  "#0D869B"],
            [0.35, "#5bb8c9"],
            [0.5,  "#f0f0f0"],
            [0.65, "#f5a96a"],
            [1.0,  "#E66F02"],
        ],
        zmid=1.0,
        zmin=0.0,
        zmax=2.0,
        colorbar=dict(
            title=dict(text="Normalized<br>infection", font=dict(size=10, color="#e6edf3")),
            tickfont=dict(color="#e6edf3", size=9),
            thickness=12,
            len=0.6,
            tickvals=[0, 0.5, 1.0, 1.5, 2.0],
            ticktext=["0", "0.5", "1.0", "1.5", "≥2.0"],
        ),
        hovertemplate="<b>%{y}</b> — %{x}<br>Value: %{z:.3f}<extra></extra>",
        xgap=2,
        ygap=1,
    ))

    fig_heat.update_layout(
        plot_bgcolor="#0d1117",
        paper_bgcolor="#0d1117",
        font=dict(color="#e6edf3"),
        height=max(400, len(genes) * 22 + 100),
        margin=dict(l=80, r=100, t=30, b=60),
        xaxis=dict(side="top", tickfont=dict(size=10), tickangle=-30),
        yaxis=dict(tickfont=dict(size=9), autorange="reversed"),
    )

    st.plotly_chart(fig_heat, use_container_width=True)