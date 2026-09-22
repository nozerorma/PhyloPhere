#!/usr/bin/env python3
"""
generate_ensembl_mapping.py  —  Auto-generate the gene Ensembl mapping file
from a gene list, via a BioMart query against Ensembl human genes.

Replaces the need for a user-supplied --gene_ensembl_file in the common
case (gene names match HGNC symbols). Only genes matched by
external_gene_name are returned; genes with no BioMart hit are reported to
--unresolved and are simply absent from the output (only 'gene' and
'length' are hard-required downstream, per filter_caas_genes.py).

Model: extract_bg.py (Malignancy_Primates/Scripts/AdHoc-Scripts), using direct
Ensembl BioMart XML queries via requests.

Output (--output) schema, per validation/fixtures/tier1/pepc/build.py:
    gene  chr  start  end  strand  length  human_protein_id
Coordinate length (end - start + 1) is used as 'length', since this file
supplies per-gene genomic length, not alignment-column length.

Caching: the BioMart response is cached under --cache-dir (default:
.cache/biomart/ next to --output), keyed by the sorted gene-ID set, so
repeated runs against the same gene list don't re-hit Ensembl.
"""

import argparse
import hashlib
from io import StringIO
import os
import sys

import pandas as pd
import requests

_COLUMNS = ["gene", "chr", "start", "end", "strand", "length", "human_protein_id"]


def load_gene_list(path: str) -> list:
    genes = set()
    with open(path) as fh:
        for line in fh:
            gene = line.strip()
            if gene:
                genes.add(gene)
    return sorted(genes)


def cache_key(genes: list) -> str:
    return hashlib.sha256("\n".join(genes).encode("utf-8")).hexdigest()


def query_biomart(genes: list, dataset: str = "hsapiens_gene_ensembl", chunk_size: int = 500) -> pd.DataFrame:
    url = "http://www.ensembl.org/biomart/martservice"
    frames = []
    for i in range(0, len(genes), chunk_size):
        chunk = genes[i : i + chunk_size]
        gene_list_str = ",".join(chunk)
        query_xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" datasetConfigVersion="0.6">
    <Dataset name="{dataset}" interface="default">
        <Filter name="external_gene_name" value="{gene_list_str}"/>
        <Attribute name="external_gene_name"/>
        <Attribute name="chromosome_name"/>
        <Attribute name="start_position"/>
        <Attribute name="end_position"/>
        <Attribute name="strand"/>
        <Attribute name="ensembl_peptide_id"/>
    </Dataset>
</Query>"""
        resp = requests.get(url, params={"query": query_xml}, timeout=120)
        resp.raise_for_status()
        if "Query ERROR" in resp.text:
            raise RuntimeError(resp.text.strip())
        chunk_df = pd.read_csv(StringIO(resp.text), sep="\t")
        chunk_df = chunk_df.rename(columns={
            "Gene name": "gene",
            "Chromosome/scaffold name": "chr",
            "Gene start (bp)": "start",
            "Gene end (bp)": "end",
            "Strand": "strand",
            "Protein stable ID": "human_protein_id",
        })
        frames.append(chunk_df)

    if not frames:
        return pd.DataFrame(columns=["gene", "chr", "start", "end", "strand", "human_protein_id"])
    return pd.concat(frames, ignore_index=True)


def derive_dataset(ref_species: str) -> str:
    """e.g. 'Homo_sapiens' -> 'hsapiens_gene_ensembl'"""
    parts = ref_species.strip().lower().split("_")
    if len(parts) >= 2:
        return f"{parts[0][0]}{parts[1]}_gene_ensembl"
    return "hsapiens_gene_ensembl"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-list", required=True,
                         help="File with one gene name per line (e.g. alignment file basenames)")
    parser.add_argument("--output", required=True, help="Output TSV")
    parser.add_argument("--unresolved", required=True,
                         help="Output file listing genes with no BioMart hit")
    parser.add_argument("--dataset", default="",
                         help="BioMart dataset name (e.g. hsapiens_gene_ensembl)")
    parser.add_argument("--ref-species", default="Homo_sapiens",
                         help="Reference species name for dataset derivation (default: Homo_sapiens)")
    parser.add_argument("--cache-dir", default=None,
                         help="Directory to cache raw BioMart responses (default: alongside --output)")
    args = parser.parse_args()

    dataset = args.dataset or derive_dataset(args.ref_species)

    genes = load_gene_list(args.gene_list)
    if not genes:
        print(f"Error: no genes found in {args.gene_list}", file=sys.stderr)
        sys.exit(1)

    cache_dir = args.cache_dir or os.path.join(os.path.dirname(os.path.abspath(args.output)), ".cache", "biomart")
    os.makedirs(cache_dir, exist_ok=True)
    cache_file = os.path.join(cache_dir, f"{cache_key(genes)}_{dataset}.tsv")

    if os.path.exists(cache_file):
        print(f"Using cached BioMart response: {cache_file}", file=sys.stderr)
        df = pd.read_csv(cache_file, sep="\t")
    else:
        print(f"Querying Ensembl BioMart ({dataset}) for {len(genes)} genes...", file=sys.stderr)
        try:
            df = query_biomart(genes, dataset=dataset)
        except Exception as exc:
            print(f"Error: Ensembl BioMart query failed ({exc}). The service may be "
                  "temporarily unavailable — retry later, or supply --gene_ensembl_file "
                  "directly.", file=sys.stderr)
            sys.exit(1)
        df.to_csv(cache_file, sep="\t", index=False)

    df = df.drop_duplicates(subset="gene", keep="first")
    df["length"] = (df["end"] - df["start"] + 1).astype(int)
    df["strand"] = df["strand"].map({1: "+", -1: "-"}).fillna(df["strand"])
    df = df[_COLUMNS]
    df.to_csv(args.output, sep="\t", index=False)

    resolved_genes = set(df["gene"])
    unresolved = [g for g in genes if g not in resolved_genes]
    with open(args.unresolved, "w") as fh:
        fh.write("\n".join(unresolved) + ("\n" if unresolved else ""))

    print(f"Resolved {len(resolved_genes)}/{len(genes)} genes via BioMart.", file=sys.stderr)
    if unresolved:
        print(f"WARNING: {len(unresolved)} genes unresolved, see {args.unresolved}", file=sys.stderr)


if __name__ == "__main__":
    main()
