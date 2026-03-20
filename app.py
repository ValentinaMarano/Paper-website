import streamlit as st

[data-testid="stSidebarNavItems"] li:first-child a span {
    visibility: hidden;
    position: relative;
}
[data-testid="stSidebarNavItems"] li:first-child a span::after {
    content: "Home";
    visibility: visible;
    position: absolute;
    left: 0;
}

st.set_page_config(
    page_title="Organelle Proteomics Atlas",
    page_icon="🧬",
    layout="wide",
)

# ── Palette ───────────────────────────────────────────────────────────────────
TEAL   = "#0D869B"
ORANGE = "#E66F02"

# ── Custom CSS ────────────────────────────────────────────────────────────────
st.markdown(f"""
<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@300;400;500;600&family=DM+Serif+Display:ital@0;1&display=swap');

html, body, [class*="css"] {{
    font-family: 'DM Sans', sans-serif;
}}

.block-container {{
    padding-top: 2rem;
    padding-bottom: 4rem;
    max-width: 900px;
}}

.hero-title {{
    font-family: 'DM Serif Display', serif;
    font-size: clamp(2rem, 5vw, 3.2rem);
    line-height: 1.15;
    margin-bottom: 0.3rem;
    color: var(--text-color);
}}
.hero-title span.teal   {{ color: {TEAL}; }}
.hero-title span.orange {{ color: {ORANGE}; }}

.hero-sub {{
    font-size: 1.05rem;
    font-weight: 300;
    opacity: 0.75;
    margin-bottom: 1.8rem;
    line-height: 1.6;
}}

.pill {{
    display: inline-block;
    padding: 4px 14px;
    border-radius: 20px;
    font-size: 0.78rem;
    font-weight: 500;
    margin-right: 6px;
    margin-bottom: 6px;
}}
.pill-teal   {{ background: {TEAL}22;   color: {TEAL};   border: 1px solid {TEAL}55; }}
.pill-orange {{ background: {ORANGE}22; color: {ORANGE}; border: 1px solid {ORANGE}55; }}
.pill-grey   {{ background: #88888822;  color: #888;     border: 1px solid #88888855; }}

.stat-grid {{
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: 16px;
    margin: 2rem 0;
}}
.stat-card {{
    border-radius: 12px;
    padding: 20px 24px;
    border-left: 4px solid {TEAL};
    background: {TEAL}0D;
}}
.stat-card.orange {{ border-left-color: {ORANGE}; background: {ORANGE}0D; }}
.stat-number {{
    font-family: 'DM Serif Display', serif;
    font-size: 2.4rem;
    color: {TEAL};
    line-height: 1;
    margin-bottom: 4px;
}}
.stat-card.orange .stat-number {{ color: {ORANGE}; }}
.stat-label {{
    font-size: 0.82rem;
    font-weight: 500;
    opacity: 0.7;
    text-transform: uppercase;
    letter-spacing: 0.06em;
}}

.section-title {{
    font-family: 'DM Serif Display', serif;
    font-size: 1.5rem;
    margin-bottom: 0.8rem;
    margin-top: 2.5rem;
}}

.nav-grid {{
    display: grid;
    grid-template-columns: repeat(2, 1fr);
    gap: 16px;
    margin-top: 1rem;
}}
.nav-card {{
    border-radius: 12px;
    padding: 22px 24px;
    border: 1px solid #88888833;
}}
.nav-card-icon  {{ font-size: 1.8rem; margin-bottom: 8px; }}
.nav-card-title {{
    font-weight: 600;
    font-size: 1rem;
    margin-bottom: 4px;
    color: {TEAL};
}}
.nav-card-desc  {{ font-size: 0.85rem; opacity: 0.65; line-height: 1.5; }}

.abstract-box {{
    border-left: 3px solid {TEAL};
    padding: 16px 22px;
    border-radius: 0 10px 10px 0;
    background: {TEAL}08;
    font-size: 0.92rem;
    line-height: 1.75;
    opacity: 0.9;
    margin-top: 1rem;
}}

.footer {{
    margin-top: 4rem;
    padding-top: 1.5rem;
    border-top: 1px solid #88888833;
    font-size: 0.78rem;
    opacity: 0.5;
    text-align: center;
    line-height: 2;
}}

.coming-soon {{
    display: inline-block;
    background: {ORANGE}22;
    color: {ORANGE};
    border: 1px solid {ORANGE}55;
    border-radius: 6px;
    padding: 3px 10px;
    font-size: 0.75rem;
    font-weight: 600;
    letter-spacing: 0.05em;
    vertical-align: middle;
    margin-left: 8px;
}}
</style>
""", unsafe_allow_html=True)

# ── Hero ──────────────────────────────────────────────────────────────────────
st.markdown("""
<div class="hero-title">
    An organelle-resolved<br>
    <span class="teal">proteomics atlas</span> reveals<br>
    host factors for <span class="orange">coronavirus</span> infection
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hero-sub">
    A spatial proteomics resource mapping the subcellular reorganisation of human cells
    during β-coronavirus HCoV-OC43 infection.
</div>
""", unsafe_allow_html=True)

st.markdown("""
<span class="pill pill-teal">Organelle Proteomics</span>
<span class="pill pill-teal">HCoV-OC43</span>
<span class="pill pill-teal">Protein Correlation Profiling</span>
<span class="pill pill-orange">SARS-CoV-2</span>
<span class="pill pill-orange">siRNA screen</span>
<span class="pill pill-grey">Mass Spectrometry</span>
<span class="pill pill-grey">UMAP</span>
""", unsafe_allow_html=True)

# ── Figure ────────────────────────────────────────────────────────────────────
st.image(
    "data/figure_home.jpeg",
    caption="Experimental strategy: organelle fractionation + LC-MS/MS + Protein Correlation Profiling",
    use_column_width=True,
)

# ── Stats ─────────────────────────────────────────────────────────────────────
st.markdown("""
<div class="stat-grid">
    <div class="stat-card">
        <div class="stat-number">773</div>
        <div class="stat-label">Host proteins relocating upon infection</div>
    </div>
    <div class="stat-card orange">
        <div class="stat-number">&gt;⅓</div>
        <div class="stat-label">Candidates confirmed as dependency or restriction factors</div>
    </div>
    <div class="stat-card">
        <div class="stat-number">224</div>
        <div class="stat-label">Host proteins enriched near the viral Replicome</div>
    </div>
</div>
""", unsafe_allow_html=True)

# ── Abstract ──────────────────────────────────────────────────────────────────
st.markdown('<div class="section-title">About this study</div>', unsafe_allow_html=True)

with st.expander("Read abstract", expanded=False):
    st.markdown("""
    <div class="abstract-box">
    Human Coronaviruses (HCoVs) are enveloped, positive-sense single-stranded RNA viruses capable of
    infecting birds, mammals, and humans. A key aspect of coronavirus replication is the formation of
    specialised compartments called viral replication organelles (ROs), primarily composed of
    double-membrane vesicles (DMVs) that originate from the endoplasmic reticulum (ER).
    <br><br>
    By combining density gradient fractionation, mass spectrometry, and computational classification
    based on Protein Correlation Profiling (PCP), we mapped the subcellular landscape of human cells
    during infection with HCoV-OC43. We identified <strong>773 host proteins</strong> that dynamically
    relocate upon infection — most of which would not have been detected by conventional whole-proteome
    analysis. A customised siRNA screen with SARS-CoV-2 revealed that over a third function as
    dependency or restriction factors.
    <br><br>
    Furthermore, by analysing the distribution of non-structural viral proteins enriched in the ER,
    we derived a novel profile termed the <strong>Replicome</strong> — a specific suborganelle ER domain
    reflecting the viral replication organelle in both protein composition and molecular function.
    </div>
    """, unsafe_allow_html=True)

# ── Explore ───────────────────────────────────────────────────────────────────
st.markdown('<div class="section-title">Explore the data</div>', unsafe_allow_html=True)

st.markdown("""
<div class="nav-grid">
    <div class="nav-card">
        <div class="nav-card-icon">🔬</div>
        <div class="nav-card-title">Traslocome</div>
        <div class="nav-card-desc">
            Explore the 773 host proteins that relocalize upon HCoV-OC43 infection.
            Visualise the UMAP of the subcellular proteome in mock vs infected conditions
            and search for your protein of interest.
        </div>
    </div>
    <div class="nav-card">
        <div class="nav-card-icon">🧫</div>
        <div class="nav-card-title">Replicome</div>
        <div class="nav-card-desc">
            Discover the 224 host proteins spatially enriched near the viral replication organelle.
            Explore their proximity to the Replicome and the ER, and search for your protein of interest.
        </div>
    </div>
</div>
""", unsafe_allow_html=True)

# ── Publication ───────────────────────────────────────────────────────────────
st.markdown('<div class="section-title">Publication</div>', unsafe_allow_html=True)

st.markdown(f"""
<div style="padding: 16px 20px; border-radius: 10px; border: 1px solid #88888833; font-size: 0.9rem; line-height: 1.7;">
    <strong>An organelle-resolved proteomics atlas reveals host factors required for coronavirus infection</strong>
    <span class="coming-soon">⏳ Coming soon</span><br>
    <span style="opacity:0.6;">Valentina Marano · Cortese Lab · TIGEM, Naples</span>
</div>
""", unsafe_allow_html=True)

# ── Footer ────────────────────────────────────────────────────────────────────
st.markdown(f"""
<div class="footer">
    <a href="https://www.tigem.it/research/research-faculty/cortese" target="_blank"
       style="color:{TEAL}; text-decoration:none; font-weight:500;">Cortese Lab</a>
    · TIGEM, Naples · 2025
    <br>
    Built by Valentina Marano
</div>
""", unsafe_allow_html=True)