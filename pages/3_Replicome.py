import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os

st.set_page_config(page_title="Replicome", page_icon="🧫", layout="wide")

TEAL   = "#0D869B"
ORANGE = "#E66F02"

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
.badge-rep {{
    background: {TEAL}22; color: {TEAL};
    border: 1px solid {TEAL}55; border-radius: 6px;
    padding: 3px 10px; font-size: 0.8rem; font-weight: 600;
}}
.badge-er {{
    background: {ORANGE}22; color: {ORANGE};
    border: 1px solid {ORANGE}55; border-radius: 6px;
    padding: 3px 10px; font-size: 0.8rem; font-weight: 600;
}}
.badge-yes {{
    background: #3fb95022; color: #3fb950;
    border: 1px solid #3fb95055; border-radius: 6px;
    padding: 3px 10px; font-size: 0.8rem; font-weight: 600;
}}
.badge-no {{
    background: #88888822; color: #888;
    border: 1px solid #88888855; border-radius: 6px;
    padding: 3px 10px; font-size: 0.8rem; font-weight: 600;
}}
</style>
""", unsafe_allow_html=True)

# ── Load data ─────────────────────────────────────────────────────────────────
@st.cache_data
def load_data():
    df = pd.read_csv("data/replicome.csv", sep=None, engine="python")
    df.columns = df.columns.str.strip()
    df = df.rename(columns={
        "Gene Name":                         "Gene",
        "Protein group":                     "Protein Group",
        "Euclidean Distance from REPLICOME": "Distance from Replicome",
        "Q-VALUE from REPLICOME":            "Q-value Replicome",
        "Euclidean Distance from ER":        "Distance from ER",
        "Q-VALUE from ER":                   "Q-value ER",
        "RISULTATO":                         "Closest Organelle",
        "DELTA DISTANZA":                    "Delta Distance",
        "Protein description":               "Protein Description",
    })
    df["Gene"] = df["Gene"].str.strip()
    df["Closest Organelle"] = df["Closest Organelle"].replace({
        "Vicino all'ER":       "Close to ER",
        "Vicino al REPLICOMA": "Close to Replicome",
        "Vicino al Replicoma": "Close to Replicome",
    })
    for col in ["Distance from Replicome","Q-value Replicome","Distance from ER","Q-value ER","Delta Distance"]:
        if df[col].dtype == object:
            df[col] = df[col].astype(str).str.replace(",",".").str.strip()
            df[col] = pd.to_numeric(df[col], errors="coerce")

    rep = df[df["Closest Organelle"] == "Close to Replicome"].copy()
    if "Quantile" not in df.columns:
        rep_sorted = rep.sort_values("Delta Distance", ascending=False)
        rep_sorted["Quantile"] = pd.qcut(rep_sorted["Delta Distance"], q=4, labels=["Q1","Q2","Q3","Q4"])
        df = df.merge(rep_sorted[["Gene","Quantile"]], on="Gene", how="left")
    df["Final Selected"] = df["Quantile"].isin(["Q2","Q3","Q4"])
    return df

df = load_data()

# ── IF channel viewer helper ──────────────────────────────────────────────────
CHANNELS = {
    "merge":   "Merge",
    "DAPI":    "DAPI",
    "dsRNA":   "dsRNA",
    "protein": "Protein",
}

def find_channel(gene, condition, channel):
    """Find image file for gene/condition/channel, any extension."""
    folder = f"data/images/{gene}"
    base = f"{gene}_{condition}_{channel}"
    for ext in [".bmp",".png",".jpg",".jpeg",".tif",".tiff"]:
        p = os.path.join(folder, base + ext)
        if os.path.exists(p):
            return p
    return None

def blend_images(paths):
    """Screen blend multiple images (works great on black background)."""
    from PIL import Image as PILImage
    imgs = []
    for p in paths:
        try:
            arr = np.array(PILImage.open(p).convert("RGB")).astype(np.float32)
            imgs.append(arr)
        except Exception:
            pass
    if not imgs:
        return None
    h = min(i.shape[0] for i in imgs)
    w = min(i.shape[1] for i in imgs)
    imgs = [i[:h,:w] for i in imgs]
    result = np.ones((h, w, 3), dtype=np.float32)
    for img_arr in imgs:
        result *= (1.0 - img_arr / 255.0)
    blended = ((1.0 - result) * 255).clip(0,255).astype(np.uint8)
    from PIL import Image as PILImage
    return PILImage.fromarray(blended)

def if_viewer(gene, condition_label, condition_key):
    """Render channel selector + blended image for one condition."""
    available = {}
    for ch_key, ch_label in CHANNELS.items():
        p = find_channel(gene, condition_key, ch_key)
        if p:
            available[ch_key] = (ch_label, p)

    if not available:
        st.caption(f"No images available for {condition_label}.")
        return

    st.markdown(f"**{condition_label}**")
    selected = {}
    cols_ch = st.columns(len(available))
    for i, (ch_key, (ch_label, _)) in enumerate(available.items()):
        with cols_ch[i]:
            selected[ch_key] = st.checkbox(
                ch_label, value=(ch_key == "merge"),
                key=f"ch_{gene}_{condition_key}_{ch_key}"
            )

    active_paths  = [available[k][1] for k in available if selected.get(k, False)]
    active_labels = [available[k][0] for k in available if selected.get(k, False)]

    if not active_paths:
        st.caption("Select at least one channel.")
    elif len(active_paths) == 1:
        st.image(active_paths[0], caption=active_labels[0])
    else:
        blended = blend_images(active_paths)
        if blended:
            st.image(blended, caption=" + ".join(active_labels))

# ── UniProt info ──────────────────────────────────────────────────────────────
import requests

@st.cache_data
def get_uniprot_info(gene_name):
    try:
        url = (f"https://rest.uniprot.org/uniprotkb/search"
               f"?query=gene:{gene_name}+AND+organism_id:9606+AND+reviewed:true"
               f"&fields=protein_name,cc_function&format=json&size=1")
        r = requests.get(url, timeout=5)
        data = r.json()
        if data["results"]:
            result = data["results"][0]
            name = (result.get("proteinDescription",{})
                    .get("recommendedName",{})
                    .get("fullName",{})
                    .get("value",""))
            comments = result.get("comments",[])
            function = next(
                (c["texts"][0]["value"] for c in comments if c["commentType"] == "FUNCTION"),
                None
            )
            return name, function
    except Exception:
        pass
    return None, None

# ── Header ────────────────────────────────────────────────────────────────────
st.title("🧫 Replicome")
st.markdown("""
During coronavirus infection, non-structural viral proteins accumulate in a specific
suborganellar ER domain that we term the **Replicome** — a proxy for the viral replication
organelle (RO). By comparing the spatial profiles of host proteins to this viral signature,
we identified host factors enriched near the site of viral replication.

Here you can explore the full dataset of **6,149 detected proteins** and their proximity
to the Replicome and the ER.
""")
st.divider()

# ── Stats ─────────────────────────────────────────────────────────────────────
c1, c2, c3, c4 = st.columns(4)
c1.metric("Total proteins detected", "6,149")
c2.metric("Close to Replicome", "280")
c3.metric("Final selected hits", "210")
st.divider()

# ── Find your protein ─────────────────────────────────────────────────────────
st.markdown('<div class="section-title">🔍 Find your protein</div>', unsafe_allow_html=True)

search = st.text_input(
    "Search by gene name (e.g. DHCR24, TMX2, EGFR...)",
    placeholder="Type a gene name...",
).strip().upper()

if search:
    match = df[df["Gene"].str.upper() == search]

    if match.empty:
        partial = df[df["Gene"].str.upper().str.contains(search, na=False)]
        if not partial.empty:
            st.info("Exact match not found. Did you mean one of these?")
            st.dataframe(
                partial[["Gene","Protein Description","Closest Organelle"]].head(10),
                use_container_width=True, hide_index=True
            )
        else:
            st.warning(f"Protein **{search}** not found in the dataset.")
    else:
        row = match.iloc[0]

        # Check if IF images exist
        has_if = os.path.isdir(f"data/images/{row['Gene']}")

        col_info, col_img = st.columns([3, 2] if has_if else [1, 0.001])

        with col_info:
            is_rep      = row["Closest Organelle"] == "Close to Replicome"
            is_selected = row["Final Selected"]
            badge_org   = f'<span class="badge-rep">🔵 Close to Replicome</span>' if is_rep else f'<span class="badge-er">🟠 Close to ER</span>'
            badge_sel   = f'<span class="badge-yes">✅ Final selected hit</span>' if is_selected else f'<span class="badge-no">➖ Not selected</span>'

            st.markdown(f"""
            <div style="margin-bottom: 16px;">
                <span style="font-family: DM Serif Display, serif; font-size: 1.6rem; font-weight: bold;">{row['Gene']}</span>
                &nbsp;&nbsp;{badge_org}&nbsp;{badge_sel}
            </div>
            """, unsafe_allow_html=True)

            # UniProt info
            uniprot_name, uniprot_function = get_uniprot_info(row["Gene"])
            if uniprot_function:
                with st.expander("📖 Protein function (UniProt)"):
                    if uniprot_name:
                        st.caption(f"**{uniprot_name}**")
                    st.write(uniprot_function)

            if pd.notna(row.get("Protein Description")):
                st.caption(f"📖 {row['Protein Description']}")

            st.markdown('<div class="info-card">', unsafe_allow_html=True)

            def info_row(label, value):
                st.markdown(f"""
                <div class="info-row">
                    <span class="info-label">{label}</span>
                    <span class="info-value">{value}</span>
                </div>
                """, unsafe_allow_html=True)

            info_row("Protein Group", row.get("Protein Group","—"))
            info_row("Organism", row.get("Organism","—"))
            info_row("Distance from Replicome", f"{row['Distance from Replicome']:.4f}" if pd.notna(row['Distance from Replicome']) else "—")
            info_row("Q-value (Replicome)", f"{row['Q-value Replicome']:.2e}" if pd.notna(row['Q-value Replicome']) else "—")
            info_row("Distance from ER", f"{row['Distance from ER']:.4f}" if pd.notna(row['Distance from ER']) else "—")
            info_row("Q-value (ER)", f"{row['Q-value ER']:.2e}" if pd.notna(row['Q-value ER']) else "—")
            info_row("Closest Organelle", row["Closest Organelle"])
            if is_rep and pd.notna(row.get("Delta Distance")):
                info_row("Delta Distance (Rep vs ER)", f"{row['Delta Distance']:.4f}")
            if pd.notna(row.get("Quantile")):
                info_row("Quantile", str(row["Quantile"]))

            st.markdown('</div>', unsafe_allow_html=True)

            # Scatter position
            st.markdown('<div class="section-title">📍 Position in the dataset</div>', unsafe_allow_html=True)
            fig = go.Figure()
            for label, color, opacity, size in [
                ("Close to ER",        "#8b949e", 0.3, 3),
                ("Close to Replicome", TEAL,      0.6, 5),
            ]:
                sub = df[df["Closest Organelle"] == label]
                fig.add_trace(go.Scatter(
                    x=sub["Distance from Replicome"], y=sub["Distance from ER"],
                    mode="markers", name=label,
                    marker=dict(color=color, size=size, opacity=opacity),
                    text=sub["Gene"],
                    hovertemplate="<b>%{text}</b><br>Rep: %{x:.3f}<br>ER: %{y:.3f}<extra></extra>",
                ))
            fig.add_trace(go.Scatter(
                x=[row["Distance from Replicome"]], y=[row["Distance from ER"]],
                mode="markers+text", name=row["Gene"],
                marker=dict(color=ORANGE, size=14, symbol="star", line=dict(width=1.5, color="white")),
                text=[row["Gene"]], textposition="top center",
                textfont=dict(size=11, color=ORANGE), hoverinfo="skip",
            ))
            max_val = max(df["Distance from Replicome"].max(), df["Distance from ER"].max()) + 0.1
            fig.add_shape(type="line", x0=0, y0=0, x1=max_val, y1=max_val,
                          line=dict(color="#3fb950", dash="dash", width=1.5))
            fig.update_layout(
                xaxis_title="Distance from Replicome", yaxis_title="Distance from ER",
                plot_bgcolor="#0d1117", paper_bgcolor="#0d1117",
                font=dict(color="#e6edf3"), height=380,
                margin=dict(l=60, r=30, t=30, b=60),
                xaxis=dict(gridcolor="#21262d", zerolinecolor="#484f58"),
                yaxis=dict(gridcolor="#21262d", zerolinecolor="#484f58"),
                legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            )
            st.plotly_chart(fig, use_container_width=True)

        # Right column: IF viewer with mock vs infected
        if has_if:
            with col_img:
                st.markdown('<div class="section-title">🔬 Immunofluorescence</div>', unsafe_allow_html=True)

                tab_mock, tab_inf = st.tabs(["🔘 Mock", "🦠 Infected (SARS-CoV-2)"])
                with tab_mock:
                    if_viewer(row["Gene"], "Mock", "mock")
                with tab_inf:
                    if_viewer(row["Gene"], "Infected", "inf")

else:
    # Browse all proteins
    st.markdown('<div class="section-title">📊 Browse all proteins</div>', unsafe_allow_html=True)

    col1, col2 = st.columns(2)
    with col1:
        filter_org = st.selectbox("Filter by closest organelle", ["All","Close to Replicome","Close to ER"])
    with col2:
        filter_sel = st.selectbox("Filter by selection", ["All","Final selected hits only"])

    display = df.copy()
    if filter_org != "All":
        display = display[display["Closest Organelle"] == filter_org]
    if filter_sel == "Final selected hits only":
        display = display[display["Final Selected"] == True]

    show_cols = ["Gene","Protein Description","Distance from Replicome",
                 "Q-value Replicome","Distance from ER","Q-value ER",
                 "Closest Organelle","Delta Distance","Quantile","Final Selected"]
    show_cols = [c for c in show_cols if c in display.columns]

    st.dataframe(display[show_cols].reset_index(drop=True), use_container_width=True, hide_index=True)
    st.caption(f"Showing {len(display):,} proteins")