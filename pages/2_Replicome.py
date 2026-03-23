import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import os

st.set_page_config(page_title="Replicome", page_icon="🧫", layout="wide")

TEAL   = "#0D869B"
ORANGE = "#E66F02"

# ── Custom CSS ────────────────────────────────────────────────────────────────
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
.info-label {{
    font-weight: 500;
    opacity: 0.6;
    min-width: 220px;
}}
.info-value {{ font-weight: 600; }}

.badge-rep {{
    background: {TEAL}22;
    color: {TEAL};
    border: 1px solid {TEAL}55;
    border-radius: 6px;
    padding: 3px 10px;
    font-size: 0.8rem;
    font-weight: 600;
}}
.badge-er {{
    background: {ORANGE}22;
    color: {ORANGE};
    border: 1px solid {ORANGE}55;
    border-radius: 6px;
    padding: 3px 10px;
    font-size: 0.8rem;
    font-weight: 600;
}}
.badge-yes {{
    background: #3fb95022;
    color: #3fb950;
    border: 1px solid #3fb95055;
    border-radius: 6px;
    padding: 3px 10px;
    font-size: 0.8rem;
    font-weight: 600;
}}
.badge-no {{
    background: #88888822;
    color: #888;
    border: 1px solid #88888855;
    border-radius: 6px;
    padding: 3px 10px;
    font-size: 0.8rem;
    font-weight: 600;
}}
</style>
""", unsafe_allow_html=True)

# ── Load data ─────────────────────────────────────────────────────────────────
@st.cache_data
def load_data():
    df = pd.read_csv("data/replicome.csv", sep=None, engine="python")
    df.columns = df.columns.str.strip()
    df = df.rename(columns={
        "Gene Name":                        "Gene",
        "Protein group":                    "Protein Group",
        "Euclidean Distance from REPLICOME":"Distance from Replicome",
        "Q-VALUE from REPLICOME":           "Q-value Replicome",
        "Euclidean Distance from ER":       "Distance from ER",
        "Q-VALUE from ER":                  "Q-value ER",
        "RISULTATO":                        "Closest Organelle",
        "DELTA DISTANZA":                   "Delta Distance",
        "Protein description":              "Protein Description",
    })
    df["Gene"] = df["Gene"].str.strip()
    df["Closest Organelle"] = df["Closest Organelle"].replace({
        "Vicino all'ER":       "Close to ER",
        "Vicino al REPLICOMA": "Close to Replicome",
        "Vicino al Replicoma": "Close to Replicome",
    })

    # Converti colonne numeriche (virgola → punto)
    for col in ["Distance from Replicome", "Q-value Replicome",
                "Distance from ER", "Q-value ER", "Delta Distance"]:
        if df[col].dtype == object:
            df[col] = df[col].astype(str).str.replace(",", ".").str.strip()
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # Final selected: proteine nel quantile 2-3-4
    rep = df[df["Closest Organelle"] == "Close to Replicome"].copy()
    if "Quantile" not in df.columns:
        rep_sorted = rep.sort_values("Delta Distance", ascending=False)
        rep_sorted["Quantile"] = pd.qcut(
            rep_sorted["Delta Distance"], q=4,
            labels=["Q1", "Q2", "Q3", "Q4"]
        )
        df = df.merge(rep_sorted[["Gene", "Quantile"]], on="Gene", how="left")
    df["Final Selected"] = df["Quantile"].isin(["Q2", "Q3", "Q4"])

    return df

df = load_data()

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
close_rep = df[df["Closest Organelle"] == "Close to Replicome"]
selected  = df[df["Final Selected"] == True]

c1, c2, c3, c4 = st.columns(4)
c1.metric("Total proteins detected", "6,149")
c2.metric("Close to Replicome", "280")
c3.metric("Final selected hits", "210")
c4.metric("Proteins with IF validation", "5")

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
        # Cerca parziale
        partial = df[df["Gene"].str.upper().str.contains(search, na=False)]
        if not partial.empty:
            st.info(f"Exact match not found. Did you mean one of these?")
            st.dataframe(
                partial[["Gene", "Protein Description", "Closest Organelle"]].head(10),
                use_container_width=True, hide_index=True
            )
        else:
            st.warning(f"Protein **{search}** not found in the dataset.")
    else:
        row = match.iloc[0]

        # Layout: info a sinistra, IF a destra se disponibile
        has_image = False
        img_extensions = [".png", ".jpg", ".jpeg", ".tif", ".tiff"]
        img_path = None
        for ext in img_extensions:
            p = f"data/images/{row['Gene']}{ext}"
            if os.path.exists(p):
                has_image = True
                img_path = p
                break

        col_info, col_img = st.columns([3, 2] if has_image else [1, 0.001])

        with col_info:
            # Badge risultato
            is_rep = row["Closest Organelle"] == "Close to Replicome"
            is_selected = row["Final Selected"]

            badge_org = f'<span class="badge-rep">🔵 Close to Replicome</span>' if is_rep else f'<span class="badge-er">🟠 Close to ER</span>'
            badge_sel = f'<span class="badge-yes">✅ Final selected hit</span>' if is_selected else f'<span class="badge-no">➖ Not selected</span>'

            st.markdown(f"""
            <div style="margin-bottom: 16px;">
                <span style="font-family: DM Serif Display, serif; font-size: 1.6rem; font-weight: bold;">{row['Gene']}</span>
                &nbsp;&nbsp;{badge_org}&nbsp;{badge_sel}
            </div>
            """, unsafe_allow_html=True)

            # Descrizione proteina
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

            info_row("Protein Group", row.get("Protein Group", "—"))
            info_row("Organism", row.get("Organism", "—"))
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

        # IF image
        if has_image:
            with col_img:
                st.markdown(f'<div class="section-title">🔬 Immunofluorescence</div>', unsafe_allow_html=True)
                st.image(img_path, caption=f"{row['Gene']} — IF during HCoV-OC43 infection", use_column_width=True)

        # Mini scatter: dove si posiziona questa proteina rispetto alle altre
        st.markdown('<div class="section-title">📍 Position in the dataset</div>', unsafe_allow_html=True)

        fig = go.Figure()

        # Tutti i punti
        for label, color, opacity, size in [
            ("Close to ER",        "#8b949e", 0.3, 3),
            ("Close to Replicome", TEAL,      0.6, 5),
        ]:
            sub = df[df["Closest Organelle"] == label]
            fig.add_trace(go.Scatter(
                x=sub["Distance from Replicome"],
                y=sub["Distance from ER"],
                mode="markers",
                name=label,
                marker=dict(color=color, size=size, opacity=opacity),
                text=sub["Gene"],
                hovertemplate="<b>%{text}</b><br>Rep: %{x:.3f}<br>ER: %{y:.3f}<extra></extra>",
            ))

        # Punto cercato
        fig.add_trace(go.Scatter(
            x=[row["Distance from Replicome"]],
            y=[row["Distance from ER"]],
            mode="markers+text",
            name=row["Gene"],
            marker=dict(color=ORANGE, size=14, symbol="star",
                        line=dict(width=1.5, color="white")),
            text=[row["Gene"]],
            textposition="top center",
            textfont=dict(size=11, color=ORANGE),
            hoverinfo="skip",
        ))

        # Linea diagonale
        max_val = max(df["Distance from Replicome"].max(), df["Distance from ER"].max()) + 0.1
        fig.add_shape(type="line", x0=0, y0=0, x1=max_val, y1=max_val,
                      line=dict(color="#3fb950", dash="dash", width=1.5))

        fig.update_layout(
            xaxis_title="Distance from Replicome",
            yaxis_title="Distance from ER",
            plot_bgcolor="#0d1117",
            paper_bgcolor="#0d1117",
            font=dict(color="#e6edf3"),
            height=420,
            margin=dict(l=60, r=30, t=30, b=60),
            xaxis=dict(gridcolor="#21262d", zerolinecolor="#484f58"),
            yaxis=dict(gridcolor="#21262d", zerolinecolor="#484f58"),
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        )
        st.plotly_chart(fig, use_container_width=True)

else:
    # Mostra tabella generale quando non si cerca niente
    st.markdown('<div class="section-title">📊 Browse all proteins</div>', unsafe_allow_html=True)

    col1, col2 = st.columns(2)
    with col1:
        filter_org = st.selectbox(
            "Filter by closest organelle",
            ["All", "Close to Replicome", "Close to ER"]
        )
    with col2:
        filter_sel = st.selectbox(
            "Filter by selection",
            ["All", "Final selected hits only"]
        )

    display = df.copy()
    if filter_org != "All":
        display = display[display["Closest Organelle"] == filter_org]
    if filter_sel == "Final selected hits only":
        display = display[display["Final Selected"] == True]

    show_cols = ["Gene", "Protein Description", "Distance from Replicome",
                 "Q-value Replicome", "Distance from ER", "Q-value ER",
                 "Closest Organelle", "Delta Distance", "Quantile", "Final Selected"]
    show_cols = [c for c in show_cols if c in display.columns]

    st.dataframe(display[show_cols].reset_index(drop=True),
                 use_container_width=True, hide_index=True)

    st.caption(f"Showing {len(display):,} proteins")
