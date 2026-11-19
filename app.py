import streamlit as st
from detector import parse_fasta, detect_duplications

st.set_page_config(page_title="Gene Duplication Detector", page_icon="🧬", layout="wide")

st.title("🧬 Gene Duplication Detector")
st.markdown("Pairwise FASTA sequence identity analyzer for detecting potential gene duplications.")

SAMPLE_FASTA = """>1L2Y_1|Chain A|TC5b|synthetic construct
NLYIQWLKDGGPSSGRPPPS
>1L2Y_2|Chain A|TC5b_mutant_1|synthetic construct
NLYIQWLKDGGPSSGRPPPA
>1L2Y_3|Chain A|TC5b_divergent|synthetic construct
MLYIQWLKDAAPSSGRAAAS
>1L2Y_4|Chain A|TC5b_short|synthetic construct
NLYIQWLKDGGPSS"""

if "fasta_text" not in st.session_state:
    st.session_state.fasta_text = SAMPLE_FASTA

col1, col2 = st.columns([1.1, 0.9], gap="large")

with col1:
    st.subheader("1. Input FASTA Sequences")

    btn_col1, btn_col2 = st.columns(2)
    with btn_col1:
        if st.button("Load Sample Data", use_container_width=True):
            st.session_state.fasta_text = SAMPLE_FASTA
            st.rerun()
    with btn_col2:
        if st.button("Clear Input", use_container_width=True):
            st.session_state.fasta_text = ""
            st.rerun()

    uploaded_file = st.file_uploader("Upload FASTA file (.fasta, .fa, .fna, .txt)", type=["fasta", "fa", "fna", "txt"])
    if uploaded_file is not None:
        st.session_state.fasta_text = uploaded_file.read().decode("utf-8")

    fasta_input = st.text_area(
        "FASTA Sequences",
        value=st.session_state.fasta_text,
        height=240,
        placeholder=">gene_1\nATGCGATCGATCGATCG\n>gene_2\nATGCGATCGATCGATCA"
    )
    st.session_state.fasta_text = fasta_input

    limit = st.slider("Identity Threshold (Limit %)", min_value=1, max_value=100, value=80)
    run_clicked = st.button("Run Duplication Detection", type="primary", use_container_width=True)

with col2:
    st.subheader("2. Final Results")

    genes = parse_fasta(st.session_state.fasta_text)
    summary = detect_duplications(genes, limit=limit)

    m1, m2, m3 = st.columns(3)
    m1.metric("Genes", summary["gene_count"])
    m2.metric("Comparisons", summary["total_comparisons"])
    m3.metric("Duplicates", summary["duplicate_pairs"])

    st.divider()

    if summary["gene_count"] < 2:
        st.info("Please enter at least 2 FASTA sequences (starting with `>`) to run pairwise comparisons.")
    else:
        for item in summary["results"]:
            status = item["status"]
            pair_label = f"**{item['gene_1']}** vs **{item['gene_2']}**"
            if status == "duplicate":
                st.success(f"{pair_label} — **{item['identity']}% Identity** ({item['matches']}/{item['length']} matches). {item['message']}")
            elif status == "no_duplication":
                st.info(f"{pair_label} — **{item['identity']}% Identity** ({item['matches']}/{item['length']} matches). {item['message']}")
            else:
                st.warning(f"{pair_label} — {item['message']}")
