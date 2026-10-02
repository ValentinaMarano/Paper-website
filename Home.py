import streamlit as st

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
    font-family: 'DM Sans', sans-serif;
    font-size: clamp(2rem, 5vw, 3.2rem);
    line-height: 1.15;
    margin-bottom: 0.3rem;
    color: var(--text-color);
}}
.hero-title span.teal   {{ color: {TEAL}; }}
.hero-title span.orange {{ color: {ORANGE}; }}

.hero-sub {{
    font-size: 1.25rem;
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
    font-family: 'DM Sans', sans-serif;
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
    font-family: 'DM Sans', sans-serif;
    font-size: 1.5rem;
    margin-bottom: 0.8rem;
    margin-top: 2.5rem;
}}

.nav-grid {{
    display: grid;
    grid-template-columns: 1fr;
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
.nav-card-desc  {{ font-size: 0.9rem; opacity: 0.65; line-height: 1.5; }}

.abstract-box {{
    border-left: 3px solid {TEAL};
    padding: 16px 22px;
    border-radius: 0 10px 10px 0;
    background: {TEAL}08;
    font-size: 0.9rem;
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
    Spatial proteomics resolves <br>
    <span class="teal"> virus-host  </span> interaction at <br>
     <span class="orange">sub-organelle </span> resolution <br>
</div>
""", unsafe_allow_html=True)

st.markdown("""
<div class="hero-sub">
    A spatial proteomics resource to map virus-host interactions and identify replication-associated host factors
</div>
""", unsafe_allow_html=True)

st.markdown("""
<span class="pill pill-teal">Organelle Proteomics</span>
<span class="pill pill-teal">Protein Correlation Profiling</span>
<span class="pill pill-orange">SARS-CoV-2</span>
<span class="pill pill-orange">HCoV-OC43</span>
<span class="pill pill-orange">Arbovirus</span>
""", unsafe_allow_html=True)

# ── Figure ────────────────────────────────────────────────────────────────────
fig_path = "data/figure_home_light.png"

st.image(
    fig_path,
    caption="Experimental strategy: organelle fractionation + LC-MS/MS + Protein Correlation Profiling",
    use_container_width=True,
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
        <div class="stat-number">197</div>
        <div class="stat-label">Host proteins enriched near viral replication organelle</div>
    </div>
</div>
""", unsafe_allow_html=True)

# ── Abstract ──────────────────────────────────────────────────────────────────
st.markdown('<div class="section-title">About this study</div>', unsafe_allow_html=True)

with st.expander("Read abstract", expanded=False):
    st.markdown("""
    <div style="text-align: justify;">
    Viruses reshape host cell organization to support their replication. Here, we combined
organelle-resolved spatial proteomics with machine learning and functional screening to map
host protein redistribution during coronavirus infection through two complementary strategies.

In the first arm, we identified <strong>773 proteins</strong> that relocalize upon HCoV-OC43
infection, the <strong>translocome</strong>, largely invisible to conventional proteomics.
Functional screening of 166 candidates revealed <strong>58 conserved host factors</strong>
that promote or restrict infection of SARS-CoV-2 and other RNA viruses.

In the second arm, we used viral replication proteins as endogenous spatial landmarks to define
a replication-associated sub-organellar signature, the <strong>Replicome</strong>, identifying
<strong>197 host proteins</strong> enriched at coronavirus replication sites. Selected lipid
metabolism enzymes were functionally validated as restriction or dependency factors.
    </div>
    """, unsafe_allow_html=True)

# ── Explore ───────────────────────────────────────────────────────────────────
st.markdown('<div class="section-title">Explore the data and <span style="color: #E66F02;">find your protein of interest:</span></div>', unsafe_allow_html=True)
st.markdown("""
<div class="nav-grid">
    <div class="nav-card">
        <div class="nav-card-title">Traslocome</div>
        <div class="nav-card-desc">
            Explore the 773 host proteins that relocalize upon HCoV-OC43 infection.
            Visualise the subcellular proteome in mock vs infected conditions.
        </div>
    </div>
    <div class="nav-card">
        <div class="nav-card-title">Screen</div>
        <div class="nav-card-desc">
            Uncover the functional role of 166 translocating proteins in viral infection.
            Browse siRNA knockdown results across SARS-CoV-2 and five arboviruses,
            with immunofluorescence images for selected candidates.
        </div>
    </div>
    <div class="nav-card">
        <div class="nav-card-title">Replicome</div>
        <div class="nav-card-desc">
            Discover the 197 host proteins spatially enriched near the viral replication organelle.
            Explore their proximity to the Replicome and the ER.
        </div>
    </div>
</div>
""", unsafe_allow_html=True)

# ── Publication ───────────────────────────────────────────────────────────────
st.markdown('<div class="section-title">Publication</div>', unsafe_allow_html=True)

st.markdown(f"""
<div style="padding: 16px 20px; border-radius: 10px; border: 1px solid #88888833; font-size: 1rem; line-height: 1.7;">
    <strong>Spatial proteomics resolves virus-host interaction at sub-organelle resolution</strong>
    <span class="coming-soon"> Submitted</span><br>
    <span style="opacity:0.6;">V. Marano et al. · Cortese Lab · TIGEM</span>
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