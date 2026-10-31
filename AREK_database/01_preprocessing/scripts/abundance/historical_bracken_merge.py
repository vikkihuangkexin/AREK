#!/usr/bin/env python3

from pathlib import Path
import argparse
import re
import pandas as pd

RANK_NAME = {
    "S": "species",
    "G": "genus",
    "F": "family",
    "P": "phylum",
}

RANK_PREFIX = {
    "D": "d__",
    "K": "k__",
    "P": "p__",
    "C": "c__",
    "O": "o__",
    "F": "f__",
    "G": "g__",
    "S": "s__",
}

MAIN_RANKS = ["D", "K", "P", "C", "O", "F", "G", "S"]


def clean_name(value):
    value = str(value).strip()
    return re.sub(r"\s+", "_", value)


def parse_kraken_report(path):
    lineage_map = {}
    stack = []

    with open(path, "r", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            parts = line.split("\t")
            if len(parts) < 6:
                continue

            rank_code = parts[3].strip()
            taxid = parts[4].strip()
            raw_name = parts[5]
            name = raw_name.strip()

            indent = len(raw_name) - len(raw_name.lstrip(" "))
            depth = indent // 2
            rank_simple = rank_code if rank_code in RANK_PREFIX else ""  # fix: G1/F1 are intermediate ranks, not G/F

            while stack and stack[-1][0] >= depth:
                stack.pop()

            stack.append((depth, rank_simple, taxid, name))

            seen = {}
            for _, rank, _, tax_name in stack:
                if rank in RANK_PREFIX:
                    seen[rank] = clean_name(tax_name)

            lineage = [
                RANK_PREFIX[rank] + seen[rank]
                for rank in MAIN_RANKS
                if rank in seen
            ]

            if lineage:
                lineage_map[taxid] = ";".join(lineage)

    return lineage_map


def read_bracken(path):
    try:
        frame = pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return pd.DataFrame()

    required = [
        "name",
        "taxonomy_id",
        "taxonomy_lvl",
        "kraken_assigned_reads",
        "added_reads",
        "new_est_reads",
        "fraction_total_reads",
    ]

    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{path} missing columns: {missing}")

    frame["taxonomy_id"] = frame["taxonomy_id"].astype(str)
    frame["Feature ID"] = "taxid_" + frame["taxonomy_id"]
    frame["new_est_reads"] = pd.to_numeric(
        frame["new_est_reads"], errors="coerce"
    ).fillna(0)
    frame["fraction_total_reads"] = pd.to_numeric(
        frame["fraction_total_reads"], errors="coerce"
    ).fillna(0)
    return frame


def merge_tables(tables):
    merged = tables[0]
    for table in tables[1:]:
        merged = merged.merge(table, on="Feature ID", how="outer")
    return merged.fillna(0)


def write_level(
    level,
    sample_ids,
    bracken_dir,
    kraken_report_dir,
    export_dir,
    tax_dir,
):
    rank_name = RANK_NAME[level]
    level_dir = bracken_dir / level

    count_tables = []
    relative_tables = []
    taxonomy_records = {}
    long_tables = []
    lineage_map = {}
    used_samples = []

    for sample in sample_ids:
        report = kraken_report_dir / f"{sample}.kreport"
        if report.exists() and report.stat().st_size > 0:
            lineage_map.update(parse_kraken_report(report))

    for sample in sample_ids:
        bracken_file = level_dir / f"{sample}.bracken.tsv"

        if not bracken_file.exists() or bracken_file.stat().st_size == 0:
            print(f"[WARNING] Missing Bracken file: {bracken_file}")
            continue

        frame = read_bracken(bracken_file)
        if frame.empty:
            print(f"[WARNING] Empty Bracken table: {bracken_file}")
            continue

        used_samples.append(sample)
        frame["SampleID"] = sample

        count_tables.append(
            frame[["Feature ID", "new_est_reads"]].rename(
                columns={"new_est_reads": sample}
            )
        )
        relative_tables.append(
            frame[["Feature ID", "fraction_total_reads"]].rename(
                columns={"fraction_total_reads": sample}
            )
        )

        for _, row in frame.iterrows():
            feature_id = str(row["Feature ID"])
            taxid = str(row["taxonomy_id"])
            name = str(row["name"])
            taxon = lineage_map.get(taxid)

            if not taxon:
                prefix = {
                    "S": "s__",
                    "G": "g__",
                    "F": "f__",
                    "P": "p__",
                }[level]
                taxon = prefix + clean_name(name)

            taxonomy_records[feature_id] = {
                "Feature ID": feature_id,
                "Taxon": taxon,
                "Confidence": "",
                "TaxID": taxid,
                "Name": name,
                "Rank": level,
            }

        frame["Taxon"] = frame["taxonomy_id"].map(lineage_map)
        long_tables.append(frame)

    if not count_tables:
        raise RuntimeError(f"No non-empty Bracken tables found for level {level}")

    count_merged = merge_tables(count_tables)
    relative_merged = merge_tables(relative_tables)
    taxonomy = pd.DataFrame(
        list(taxonomy_records.values())
    ).sort_values(["Taxon", "Feature ID"])
    long_all = pd.concat(long_tables, ignore_index=True)

    feature_order = taxonomy["Feature ID"].tolist()

    count_merged = (
        count_merged.set_index("Feature ID")
        .reindex(feature_order)
        .fillna(0)
        .reset_index()
    )
    relative_merged = (
        relative_merged.set_index("Feature ID")
        .reindex(feature_order)
        .fillna(0)
        .reset_index()
    )

    if level == "S":
        count_out = export_dir / "feature-table.tsv"
        relative_out = export_dir / "feature-table-relative.tsv"
        taxonomy_out = tax_dir / "taxonomy.tsv"
        qiime_taxonomy_out = tax_dir / "taxonomy_qiime_import.tsv"
        long_out = export_dir / "long-table.tsv"
        sample_count_out = export_dir / "sample-by-feature-table.tsv"
        sample_relative_out = export_dir / "sample-by-feature-table-relative.tsv"
    else:
        count_out = export_dir / f"{rank_name}-feature-table.tsv"
        relative_out = export_dir / f"{rank_name}-feature-table-relative.tsv"
        taxonomy_out = tax_dir / f"{rank_name}-taxonomy.tsv"
        qiime_taxonomy_out = tax_dir / f"{rank_name}-taxonomy_qiime_import.tsv"
        long_out = export_dir / f"{rank_name}-long-table.tsv"
        sample_count_out = (
            export_dir / f"{rank_name}-sample-by-feature-table.tsv"
        )
        sample_relative_out = (
            export_dir / f"{rank_name}-sample-by-feature-table-relative.tsv"
        )

    count_merged.rename(
        columns={"Feature ID": "#OTU ID"}
    ).to_csv(count_out, sep="\t", index=False)

    relative_merged.rename(
        columns={"Feature ID": "#OTU ID"}
    ).to_csv(relative_out, sep="\t", index=False)

    taxonomy[
        ["Feature ID", "Taxon", "Confidence"]
    ].to_csv(taxonomy_out, sep="\t", index=False)

    taxonomy[
        ["Feature ID", "Taxon"]
    ].to_csv(
        qiime_taxonomy_out,
        sep="\t",
        index=False,
        header=False,
    )

    long_all.to_csv(long_out, sep="\t", index=False)

    count_sample_by_feature = count_merged.set_index("Feature ID").T
    count_sample_by_feature.index.name = "SampleID"
    count_sample_by_feature.to_csv(sample_count_out, sep="\t")

    relative_sample_by_feature = relative_merged.set_index("Feature ID").T
    relative_sample_by_feature.index.name = "SampleID"
    relative_sample_by_feature.to_csv(sample_relative_out, sep="\t")

    print(
        f"[INFO] {rank_name}: "
        f"{len(feature_order)} features across {len(used_samples)} samples"
    )
    print(f"[INFO] Wrote: {count_out}")
    print(f"[INFO] Wrote: {relative_out}")
    print(f"[INFO] Wrote: {taxonomy_out}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--bracken_dir", required=True)
    parser.add_argument("--kraken_report_dir", required=True)
    parser.add_argument("--export_dir", required=True)
    parser.add_argument("--tax_dir", required=True)
    parser.add_argument(
        "--levels",
        nargs="+",
        default=["S", "G", "F", "P"],
    )
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t")
    sample_ids = manifest["sample_id"].astype(str).tolist()

    bracken_dir = Path(args.bracken_dir)
    kraken_report_dir = Path(args.kraken_report_dir)
    export_dir = Path(args.export_dir)
    tax_dir = Path(args.tax_dir)

    export_dir.mkdir(parents=True, exist_ok=True)
    tax_dir.mkdir(parents=True, exist_ok=True)

    for level in args.levels:
        write_level(
            level,
            sample_ids,
            bracken_dir,
            kraken_report_dir,
            export_dir,
            tax_dir,
        )


if __name__ == "__main__":
    main()
