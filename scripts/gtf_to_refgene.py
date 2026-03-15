#!/usr/bin/env python3
"""Convert an NCBI GTF annotation file to CoRAL's refGene format.

Designed for assemblies whose GTF uses RefSeq accession IDs (NC_xxxxxx.x)
as chromosome names, such as T2T-CHM13v2.0.  A sequence report TSV (from
the NCBI Datasets page or FTP) provides the authoritative accession → UCSC
chromosome name mapping.

Output format matches the UCSC refGene table: one row per transcript, with
real exon coordinates in columns 9/10.  CoRAL's parse_genes deduplicates by
gene name so only the first transcript per gene is used for plotting, but
all transcripts are included for compatibility.

Usage
-----
python scripts/gtf_to_refgene.py \\
    --gtf     ~/Downloads/GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf \\
    --seqrep  ~/Downloads/sequence_report.tsv \\
    --output  coral/supplemental_data/refGene_t2t.txt

The sequence report TSV is downloadable from NCBI Datasets with columns
including "RefSeq seq accession" and "UCSC style name".  The one used to
generate refGene_t2t.txt was:
  GCF_009914755.1  sequence_report.tsv  (T2T-CHM13v2.0, retrieved 2026-04)
"""

import argparse
import pathlib
import re
import sys
from collections import defaultdict

PSEUDOGENE_BIOTYPES = {
    "pseudogene",
    "transcribed_pseudogene",
    "ncRNA_pseudogene",
}

_ATTR_RE = re.compile(r'(\w+)\s+"([^"]*)"')


def parse_attributes(attr_str: str) -> dict[str, str]:
    return dict(_ATTR_RE.findall(attr_str))


def parse_sequence_report(path: pathlib.Path) -> dict[str, str]:
    """Return {refseq_accession: ucsc_name} from an NCBI sequence report TSV."""
    mapping: dict[str, str] = {}
    with path.open() as fh:
        header = next(fh).rstrip("\n").split("\t")
        try:
            acc_col = header.index("RefSeq seq accession")
            ucsc_col = header.index("UCSC style name")
        except ValueError as e:
            sys.exit(
                f"Sequence report is missing expected column: {e}\n"
                f"Columns found: {header}"
            )
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            acc = fields[acc_col].strip()
            ucsc = fields[ucsc_col].strip()
            if acc and ucsc:
                mapping[acc] = ucsc
    return mapping


def hpv16_lines(hg_refgene: pathlib.Path) -> list[str]:
    """Extract the HPV16 header lines prepended to the hg refGene files."""
    lines = []
    with hg_refgene.open() as fh:
        for line in fh:
            if "hpv16ref_1" in line:
                lines.append(line)
            elif lines:
                break  # hpv16 block is over
    return lines


def convert(
    gtf_path: pathlib.Path,
    seqrep_path: pathlib.Path,
    output_path: pathlib.Path,
    hg_refgene: pathlib.Path | None,
) -> None:
    acc_to_chr = parse_sequence_report(seqrep_path)
    print(
        f"Loaded {len(acc_to_chr)} accession → chromosome mappings "
        f"from {seqrep_path.name}",
        file=sys.stderr,
    )

    # transcript_id -> metadata dict
    transcripts: dict[str, dict] = {}
    # transcript_id -> list of (exon_start, exon_end), 0-based half-open
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    # track insertion order so output preserves GTF order
    transcript_order: list[str] = []

    skipped_acc: set[str] = set()
    n_pseudogene = 0

    with gtf_path.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            feature = fields[2]
            if feature not in ("transcript", "exon"):
                continue

            acc = fields[0]
            chr_name = acc_to_chr.get(acc)
            if chr_name is None:
                skipped_acc.add(acc)
                continue

            attrs = parse_attributes(fields[8])
            transcript_id = attrs.get("transcript_id", "")
            if not transcript_id:
                continue

            # GTF is 1-based closed; refGene is 0-based half-open
            start = int(fields[3]) - 1
            end = int(fields[4])

            if feature == "transcript":
                biotype = attrs.get("gene_biotype", "")
                if biotype in PSEUDOGENE_BIOTYPES:
                    n_pseudogene += 1
                    continue
                transcripts[transcript_id] = {
                    "chrom": chr_name,
                    "strand": fields[6],
                    "tx_start": start,
                    "tx_end": end,
                    "gene_id": attrs.get("gene_id", transcript_id),
                }
                transcript_order.append(transcript_id)

            elif feature == "exon":
                # only collect exons for transcripts we're keeping
                if transcript_id in transcripts:
                    exons[transcript_id].append((start, end))

    if skipped_acc:
        print(
            f"Warning: {len(skipped_acc)} accession(s) not found in sequence "
            f"report and were skipped: {sorted(skipped_acc)}",
            file=sys.stderr,
        )
    # Deduplicate: keep first transcript seen per gene (matching parse_genes behaviour)
    seen_genes: set[str] = set()
    deduped_order: list[str] = []
    for tid in transcript_order:
        gene_id = transcripts[tid]["gene_id"]
        if gene_id not in seen_genes:
            seen_genes.add(gene_id)
            deduped_order.append(tid)

    print(f"Skipped {n_pseudogene} pseudogene transcripts.", file=sys.stderr)
    print(
        f"Writing {len(deduped_order)} gene entries (first transcript per gene) "
        f"to {output_path}",
        file=sys.stderr,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as out:
        if hg_refgene is not None:
            hpv = hpv16_lines(hg_refgene)
            if hpv:
                out.writelines(hpv)
                print(
                    f"Prepended {len(hpv)} HPV16 lines from {hg_refgene.name}",
                    file=sys.stderr,
                )
            else:
                print(
                    f"Warning: no HPV16 lines found in {hg_refgene}",
                    file=sys.stderr,
                )

        for tid in deduped_order:
            info = transcripts[tid]
            exon_list = sorted(exons.get(tid, [(info["tx_start"], info["tx_end"])]))
            exon_count = len(exon_list)
            exon_starts = ",".join(str(s) for s, _ in exon_list) + ","
            exon_ends = ",".join(str(e) for _, e in exon_list) + ","
            exon_frames = ",".join(["-1"] * exon_count) + ","

            row = (
                f"0\t{tid}\t{info['chrom']}\t{info['strand']}\t"
                f"{info['tx_start']}\t{info['tx_end']}\t"
                f"{info['tx_start']}\t{info['tx_end']}\t"
                f"{exon_count}\t{exon_starts}\t{exon_ends}\t"
                f"0\t{info['gene_id']}\tunk\tunk\t{exon_frames}"
            )
            out.write(row + "\n")


def main() -> None:
    here = pathlib.Path(__file__).parent
    default_hg = here.parent / "coral" / "supplemental_data" / "refGene_hg38.txt"

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", required=True, type=pathlib.Path,
                        help="NCBI GTF annotation file")
    parser.add_argument("--seqrep", required=True, type=pathlib.Path,
                        help="NCBI sequence report TSV with accession→UCSC name mapping")
    parser.add_argument("--output", "-o", type=pathlib.Path,
                        default=here.parent / "coral" / "supplemental_data" / "refGene_t2t.txt",
                        help="Output refGene file "
                             "(default: coral/supplemental_data/refGene_t2t.txt)")
    parser.add_argument("--hg-refgene", type=pathlib.Path, default=default_hg,
                        help="Existing hg refGene file to copy HPV16 header lines from "
                             "(default: coral/supplemental_data/refGene_hg38.txt)")
    args = parser.parse_args()

    for p, label in [(args.gtf, "--gtf"), (args.seqrep, "--seqrep")]:
        if not p.exists():
            sys.exit(f"Error: {label} file not found: {p}")

    hg_refgene = args.hg_refgene if args.hg_refgene.exists() else None
    if args.hg_refgene and not args.hg_refgene.exists():
        print(
            f"Warning: --hg-refgene file not found ({args.hg_refgene}), "
            "skipping HPV16 lines.",
            file=sys.stderr,
        )

    convert(args.gtf, args.seqrep, args.output, hg_refgene)


if __name__ == "__main__":
    main()
