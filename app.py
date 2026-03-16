import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go

# ── Configurazione pagina ─────────────────────────────────────────────────────
st.set_page_config(
    page_title="Volcano Plot",
    page_icon="🌋",
    layout="wide",
)

st.title("🌋 Volcano Plot — Analisi Proteomica")
st.markdown("Carica i tuoi dati o usa i dati di esempio per esplorare i risultati.")

# ── Sidebar: controlli ────────────────────────────────────────────────────────
with st.sidebar:
    st.header("⚙️ Parametri")

    fc_thresh = st.slider(
        "Soglia Fold Change (log₂FC)",
        min_value=0.1, max_value=5.0, value=1.5, step=0.1,
        help="Punti con |log2FC| > soglia sono considerati significativi"
    )

    p_thresh = st.slider(
        "Soglia p-value",
        min_value=0.001, max_value=0.1, value=0.05, step=0.001,
        format="%.3f",
        help="Punti con p-value < soglia sono considerati significativi"
    )

    st.divider()
    st.header("📂 Dati")
    use_example = st.toggle("Usa dati di esempio", value=True)

    uploaded_file = None
    if not use_example:
        uploaded_file = st.file_uploader(
            "Carica CSV",
            type=["csv"],
            help="Il file deve avere colonne: Gene, log2FC, pvalue"
        )
        st.caption("Colonne richieste: `Gene`, `log2FC`, `pvalue`")

# ── Genera o carica dati ──────────────────────────────────────────────────────
@st.cache_data
def generate_example_data(n=200):
    rng = np.random.default_rng(42)
    genes = [f"GENE_{i:03d}" for i in range(n)]
    # Alcuni nomi reali per rendere più realistico
    real_genes = [
        "ACTB","GAPDH","HSP90AA1","VIM","TUBB","LMNA","MYH9","FLNA","HSPA8","PKM",
        "ENO1","ALDOA","PGK1","TPI1","LDHA","MDH2","FASN","VDAC1","CANX","CALR",
        "RAB7A","LAMP1","CTSD","EGFR","ERBB2","TP53","BCL2","CASP3","IL6","TNF",
        "HIF1A","VEGFA","COL1A1","TGFB1","STAT3","AKT1","MTOR","CDK1","CCNB1","PLK1",
    ]
    genes[:len(real_genes)] = real_genes

    fc = rng.normal(0, 1.5, n)
    # Alcuni hit veri
    fc[:20] = rng.choice([-1, 1], 20) * (2.5 + rng.random(20) * 3)

    pval_raw = rng.uniform(0, 1, n)
    pval = np.where(np.abs(fc) > 2, pval_raw * 0.001, pval_raw * 0.5)
    pval = np.clip(pval, 1e-12, 1)

    return pd.DataFrame({"Gene": genes, "log2FC": fc, "pvalue": pval})


def load_data(file):
    df = pd.read_csv(file)
    required = {"Gene", "log2FC", "pvalue"}
    missing = required - set(df.columns)
    if missing:
        st.error(f"❌ Colonne mancanti nel CSV: {missing}")
        st.stop()
    return df


if use_example or uploaded_file is None:
    df = generate_example_data()
    if not use_example:
        st.warning("⚠️ Nessun file caricato — uso i dati di esempio.")
else:
    df = load_data(uploaded_file)

# ── Classificazione punti ─────────────────────────────────────────────────────
neg_log_p_thresh = -np.log10(p_thresh)
df["negLog10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))

def classify(row):
    if row["log2FC"] >= fc_thresh and row["negLog10p"] >= neg_log_p_thresh:
        return "Upregulated"
    elif row["log2FC"] <= -fc_thresh and row["negLog10p"] >= neg_log_p_thresh:
        return "Downregulated"
    else:
        return "Not significant"

df["Direction"] = df.apply(classify, axis=1)

counts = df["Direction"].value_counts()
up_n   = counts.get("Upregulated", 0)
down_n = counts.get("Downregulated", 0)
ns_n   = counts.get("Not significant", 0)

# ── Metriche rapide ───────────────────────────────────────────────────────────
col1, col2, col3, col4 = st.columns(4)
col1.metric("Totale proteine", len(df))
col2.metric("⬆️ Upregulated", up_n)
col3.metric("⬇️ Downregulated", down_n)
col4.metric("— Non significativo", ns_n)

# ── Volcano Plot ──────────────────────────────────────────────────────────────
color_map = {
    "Upregulated":     "#f85149",
    "Downregulated":   "#58a6ff",
    "Not significant": "#8b949e",
}

# Etichette solo per i top hit significativi
sig = df[df["Direction"] != "Not significant"].copy()
top_labels = sig.nlargest(20, "negLog10p")

fig = go.Figure()

# Punti NS prima (sotto)
for direction, color in color_map.items():
    subset = df[df["Direction"] == direction]
    opacity = 0.4 if direction == "Not significant" else 0.85
    size    = 5   if direction == "Not significant" else 7
    fig.add_trace(go.Scatter(
        x=subset["log2FC"],
        y=subset["negLog10p"],
        mode="markers",
        name=direction,
        marker=dict(color=color, size=size, opacity=opacity,
                    line=dict(width=0.5, color="rgba(255,255,255,0.2)")),
        text=subset["Gene"],
        customdata=np.stack([subset["pvalue"], subset["log2FC"]], axis=1),
        hovertemplate=(
            "<b>%{text}</b><br>"
            "log₂FC: %{customdata[1]:.3f}<br>"
            "p-value: %{customdata[0]:.2e}<br>"
            "-log₁₀p: %{y:.2f}<extra></extra>"
        ),
    ))

# Etichette sui top hit
fig.add_trace(go.Scatter(
    x=top_labels["log2FC"],
    y=top_labels["negLog10p"],
    mode="text",
    text=top_labels["Gene"],
    textposition="top center",
    textfont=dict(size=9, color="white"),
    showlegend=False,
    hoverinfo="skip",
))

# Linee soglia
x_range = [df["log2FC"].min() - 0.5, df["log2FC"].max() + 0.5]
y_range = [0, df["negLog10p"].max() + 0.5]

fig.add_shape(type="line", x0=fc_thresh,  x1=fc_thresh,  y0=0, y1=y_range[1],
              line=dict(color="#f85149", dash="dash", width=1.5))
fig.add_shape(type="line", x0=-fc_thresh, x1=-fc_thresh, y0=0, y1=y_range[1],
              line=dict(color="#58a6ff", dash="dash", width=1.5))
fig.add_shape(type="line", x0=x_range[0], x1=x_range[1],
              y0=neg_log_p_thresh, y1=neg_log_p_thresh,
              line=dict(color="#3fb950", dash="dash", width=1.5))

fig.update_layout(
    xaxis_title="log₂ Fold Change",
    yaxis_title="-log₁₀(p-value)",
    plot_bgcolor="#0d1117",
    paper_bgcolor="#0d1117",
    font=dict(color="#e6edf3"),
    legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    height=550,
    margin=dict(l=60, r=30, t=40, b=60),
    xaxis=dict(gridcolor="#21262d", zerolinecolor="#484f58"),
    yaxis=dict(gridcolor="#21262d", zerolinecolor="#484f58"),
)

st.plotly_chart(fig, use_container_width=True)

# ── Tabella geni significativi ────────────────────────────────────────────────
st.subheader("📋 Proteine significative")

sig_display = df[df["Direction"] != "Not significant"].copy()
sig_display = sig_display.sort_values("negLog10p", ascending=False)
sig_display["pvalue"] = sig_display["pvalue"].map("{:.2e}".format)
sig_display["log2FC"] = sig_display["log2FC"].round(3)
sig_display["negLog10p"] = sig_display["negLog10p"].round(2)
sig_display = sig_display.rename(columns={
    "log2FC": "log₂FC", "pvalue": "p-value", "negLog10p": "-log₁₀p", "Direction": "Direzione"
})

st.dataframe(
    sig_display[["Gene", "log₂FC", "p-value", "-log₁₀p", "Direzione"]],
    use_container_width=True,
    hide_index=True,
)

# Download CSV
csv = sig_display.to_csv(index=False).encode("utf-8")
st.download_button(
    label="⬇️ Scarica CSV proteine significative",
    data=csv,
    file_name="significant_proteins.csv",
    mime="text/csv",
)
