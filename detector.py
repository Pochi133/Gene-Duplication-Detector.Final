from io import StringIO
from Bio import SeqIO


def parse_fasta(fasta_text: str) -> dict[str, str]:
    """Parse FASTA formatted text into a dictionary mapping sequence IDs to uppercase sequences."""
    genes = {}
    if not fasta_text or not fasta_text.strip():
        return genes

    handle = StringIO(fasta_text.strip())
    for record in SeqIO.parse(handle, "fasta"):
        genes[record.id] = str(record.seq).upper()
    return genes


def parse_fasta_file(file_path: str) -> dict[str, str]:
    """Parse a FASTA file from disk into a dictionary mapping sequence IDs to uppercase sequences."""
    genes = {}
    for record in SeqIO.parse(file_path, "fasta"):
        genes[record.id] = str(record.seq).upper()
    return genes


def detect_duplications(genes: dict[str, str], limit: float = 80.0) -> dict:
    """Compare all gene pairs for exact position-by-position sequence identity."""
    gene_names = list(genes.keys())
    duplicates = 0
    total = 0
    results = []

    for i in range(len(gene_names)):
        for j in range(i + 1, len(gene_names)):
            name1 = gene_names[i]
            name2 = gene_names[j]
            seq1 = genes[name1]
            seq2 = genes[name2]
            total += 1

            if len(seq1) != len(seq2):
                results.append({
                    "gene_1": name1,
                    "gene_2": name2,
                    "status": "length_mismatch",
                    "identity": None,
                    "matches": 0,
                    "length": None,
                    "message": f"Sequences have different lengths ({len(seq1)} vs {len(seq2)}). Alignment is needed to compare them."
                })
            elif len(seq1) == 0:
                results.append({
                    "gene_1": name1,
                    "gene_2": name2,
                    "status": "empty_sequence",
                    "identity": None,
                    "matches": 0,
                    "length": 0,
                    "message": "Cannot compare empty sequences."
                })
            else:
                matches = sum(1 for k in range(len(seq1)) if seq1[k] == seq2[k])
                identity = round((matches / len(seq1)) * 100, 2)
                is_dup = identity >= limit
                if is_dup:
                    duplicates += 1

                results.append({
                    "gene_1": name1,
                    "gene_2": name2,
                    "status": "duplicate" if is_dup else "no_duplication",
                    "identity": identity,
                    "matches": matches,
                    "length": len(seq1),
                    "message": "Possible gene duplication detected!" if is_dup else "No potential duplication detected."
                })

    return {
        "gene_count": len(gene_names),
        "total_comparisons": total,
        "duplicate_pairs": duplicates,
        "results": results
    }
