import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import umap

st.set_page_config(page_title="Traslocome", page_icon="🔬", layout="wide")

# ── Palette organelli ─────────────────────────────────────────────────────────
ORGANELLE_COLORS = {
    "ER":               "#4e9af1",
    "PM":               "#f4a261",
    "Cytoplasm":        "#a8dadc",
    "Mitochondria":     "#e63946",
    "Nucleus":          "#9b5de5",
    "Translation":      "#f1c40f",
    "Lysosome":         "#2ec4b6",
    "protein-complex":  "#1d3557",
    "endosome-vesicle": "#e9c46a",
    "Golgi apparatus":  "#06d6a0",
    "Peroxisome":       "#ff6b6b",
    "Others":           "#adb5bd",
    "Unknown":          "#ced4da",
}

NUMERIC_COLS = [
    "endosome-vesicle", "protein-complex", "PM", "ER",
    "Cytoplasm", "Mitochondria", "Nucleus", "Translation",
    "Lysosome", "Peroxisome", "Golgi apparatus", "LD",
]

# ── Load & process data ───────────────────────────────────────────────────────
@st.cache_data
def load_and_compute_umap():
    mock = pd.read_csv("data/mock.csv")
    inf  = pd.read_csv("data/inf.csv")

    mock = mock.rename(columns={"C: Winner MOCK": "Localization"})
    inf  = inf.rename(columns={"C: Winner INF":  "Localization"})

    # Proteine comuni
    common = list(set(mock["T: T: Genes"]).intersection(set(inf["T: T: Genes"])))
    mock = mock[mock["T: T: Genes"].isin(common)].drop_duplicates("T: T: Genes").set_index("T: T: Genes").loc[common].reset_index()
    inf  = inf[inf["T: T: Genes"].isin(common)].drop_duplicates("T: T: Genes").set_index("T: T: Genes").loc[common].reset_index()

    # Raggruppa categorie rare in "Others"
    rare = ["LD", "Peroxisome", "endosome-vesicle"]
    mock["Localization"] = mock["Localization"].replace(rare, "Others")
    inf["Localization"]  = inf["Localization"].replace(rare, "Others")

    # UMAP su dati combinati per spazio comune
    mock["_cond"] = "MOCK"
    inf["_cond"]  = "INF"
    combined = pd.concat([mock, inf], ignore_index=True)

    num = combined[NUMERIC_COLS].fillna(0).values
    reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
    coords  = reducer.fit_transform(num)

    combined["UMAP1"] = coords[:, 0]
    combined["UMAP2"] = coords[:, 1]

    mock_umap = combined[combined["_cond"] == "MOCK"].copy()
    inf_umap  = combined[combined["_cond"] == "INF"].copy()

    return mock_umap, inf_umap, common

# ── UI ────────────────────────────────────────────────────────────────────────
st.title("🔬 Traslocome")
st.markdown(
    "Spatial reorganization of the host proteome upon coronavirus infection. "
    "Each dot is a protein, colored by its predicted organelle localization. "
    "Explore how proteins **relocalize** from mock to infected conditions."
)
st.divider()

with st.spinner("Computing UMAP (first load may take ~30s)..."):
    mock_umap, inf_umap, common_genes = load_and_compute_umap()

# ── Sidebar ───────────────────────────────────────────────────────────────────
with st.sidebar:
    st.header("⚙️ Options")

    view_mode = st.radio(
        "View mode",
        ["Side by side", "Animated transition"],
        index=0,
    )

    color_by = st.radio(
        "Color by",
        ["MOCK localization", "INF localization"],
        index=0,
    )

    organelles = list(ORGANELLE_COLORS.keys())
    selected_org = st.multiselect(
        "Highlight organelles",
        options=organelles,
        default=organelles,
    )

    st.divider()
    st.markdown("**Legend**")
    for org, col in ORGANELLE_COLORS.items():
        if org in selected_org:
            st.markdown(f"<span style='color:{col}'>■</span> {org}", unsafe_allow_html=True)

# ── Find your protein ─────────────────────────────────────────────────────────
st.subheader("🔍 Find your protein")
search = st.text_input("Type a gene name (e.g. ACTB, VIM, EGFR...)", "").strip().upper()

highlight_gene = None
if search:
    match = mock_umap[mock_umap["T: T: Genes"].str.upper() == search]
    if len(match) == 0:
        st.warning(f"Protein **{search}** not found in the dataset.")
    else:
        highlight_gene = search
        row_mock = mock_umap[mock_umap["T: T: Genes"].str.upper() == search].iloc[0]
        row_inf  = inf_umap[inf_umap["T: T: Genes"].str.upper() == search].iloc[0]

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Gene", row_mock["T: T: Genes"])
        c2.metric("Localization MOCK", row_mock["Localization"])
        c3.metric("Localization INF",  row_inf["Localization"])
        relocated = row_mock["Localization"] != row_inf["Localization"]
        c4.metric("Relocalized?", "✅ YES" if relocated else "➖ NO")

        if "T: T: First Protein Description" in row_mock:
            st.caption(f"📖 {row_mock['T: T: First Protein Description']}")

st.divider()

# ── Helper: build scatter trace ───────────────────────────────────────────────
def make_trace(df, color_col, title, highlight=None):
    traces = []
    for org in selected_org:
        subset = df[df[color_col] == org]
        if subset.empty:
            continue
        is_hl = (highlight is not None) and (subset["T: T: Genes"].str.upper() == highlight.upper()).any()
        traces.append(go.Scatter(
            x=subset["UMAP1"],
            y=subset["UMAP2"],
            mode="markers",
            name=org,
            marker=dict(
                color=ORGANELLE_COLORS.get(org, "#adb5bd"),
                size=4,
                opacity=0.7,
                line=dict(width=0),
            ),
            text=subset["T: T: Genes"],
            customdata=subset[[color_col, "T: T: First Protein Description"]].values,
            hovertemplate=(
                "<b>%{text}</b><br>"
                "Localization: %{customdata[0]}<br>"
                "<i>%{customdata[1]}</i><extra></extra>"
            ),
            showlegend=(title == "MOCK"),  # mostra legenda solo una volta
        ))

    # Highlight proteina cercata
    if highlight:
        hl = df[df["T: T: Genes"].str.upper() == highlight.upper()]
        if not hl.empty:
            row = hl.iloc[0]
            traces.append(go.Scatter(
                x=[row["UMAP1"]],
                y=[row["UMAP2"]],
                mode="markers+text",
                name=f"★ {highlight}",
                marker=dict(color="white", size=14, symbol="star",
                            line=dict(width=1.5, color="black")),
                text=[row["T: T: Genes"]],
                textposition="top center",
                textfont=dict(size=11, color="white"),
                showlegend=True,
                hoverinfo="skip",
            ))
    return traces

# ── Color column ──────────────────────────────────────────────────────────────
color_col_mock = "Localization"
color_col_inf  = "Localization"
color_col = color_col_mock if color_by == "MOCK localization" else color_col_inf

# ── View: Side by side ────────────────────────────────────────────────────────
if view_mode == "Side by side":
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["MOCK (uninfected)", "INFECTED"],
        horizontal_spacing=0.05,
    )

    for trace in make_trace(mock_umap, color_col, "MOCK", highlight_gene):
        fig.add_trace(trace, row=1, col=1)
    for trace in make_trace(inf_umap, color_col, "INF", highlight_gene):
        t = trace
        t.showlegend = False
        fig.add_trace(t, row=1, col=2)

    fig.update_layout(
        plot_bgcolor="#0d1117",
        paper_bgcolor="#0d1117",
        font=dict(color="#e6edf3"),
        height=560,
        margin=dict(l=40, r=40, t=60, b=40),
        legend=dict(
            orientation="v", x=1.01, y=1,
            font=dict(size=10),
            itemsizing="constant",
        ),
        hovermode="closest",
    )
    fig.update_xaxes(showgrid=False, zeroline=False, showticklabels=False)
    fig.update_yaxes(showgrid=False, zeroline=False, showticklabels=False)
    fig.update_annotations(font=dict(size=13, color="#e6edf3"))

    st.plotly_chart(fig, use_container_width=True)

# ── View: Animated transition ─────────────────────────────────────────────────
else:
    st.markdown("Use the **slider** to transition from MOCK to INFECTED.")

    alpha = st.slider("MOCK ← → INFECTED", 0.0, 1.0, 0.0, 0.01)

    # Interpola posizioni
    anim_df = mock_umap.copy()
    anim_df["UMAP1"] = (1 - alpha) * mock_umap["UMAP1"].values + alpha * inf_umap["UMAP1"].values
    anim_df["UMAP2"] = (1 - alpha) * mock_umap["UMAP2"].values + alpha * inf_umap["UMAP2"].values
    # Localizzazione: cambia a INF quando alpha > 0.5
    anim_df["Localization"] = np.where(
        alpha > 0.5,
        inf_umap["Localization"].values,
        mock_umap["Localization"].values,
    )

    fig2 = go.Figure()
    for trace in make_trace(anim_df, "Localization", "ANIM", highlight_gene):
        fig2.add_trace(trace)

    label = f"{'MOCK' if alpha < 0.5 else 'INFECTED'} ({int(alpha*100)}%)"
    fig2.update_layout(
        title=dict(text=label, font=dict(size=14, color="#e6edf3")),
        plot_bgcolor="#0d1117",
        paper_bgcolor="#0d1117",
        font=dict(color="#e6edf3"),
        height=560,
        margin=dict(l=40, r=160, t=60, b=40),
        legend=dict(orientation="v", x=1.01, y=1, font=dict(size=10), itemsizing="constant"),
        hovermode="closest",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
    )

    st.plotly_chart(fig2, use_container_width=True)

# ── Stats ─────────────────────────────────────────────────────────────────────
st.divider()
with st.expander("📊 Localization statistics"):
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**MOCK**")
        st.dataframe(
            mock_umap["Localization"].value_counts().reset_index().rename(
                columns={"index": "Organelle", "Localization": "Count"}),
            use_container_width=True, hide_index=True,
        )
    with col2:
        st.markdown("**INFECTED**")
        st.dataframe(
            inf_umap["Localization"].value_counts().reset_index().rename(
                columns={"index": "Organelle", "Localization": "Count"}),
            use_container_width=True, hide_index=True,
        )