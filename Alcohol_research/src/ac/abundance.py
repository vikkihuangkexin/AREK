"""Supported-rank aggregation, with explicit scale and no invented lineages."""
from .inputs import *

def build(args, out):
    if not args.profile:
        if args.rank != "genus":
            raise ValueError("Frozen canonical panels support genus-level features only; supply genuine rank evidence.")
        t, _, family = Inputs(args.input_root).panel(args.panel)
        save(t, out/"genus_relative_abundance.tsv.gz", True)
        save(family, out/"taxonomy.tsv")
        return dict(stage="integrated_genus_profile", semantics="relative_original_TSS",
            counts="not_inferred", ranks=["genus"], specimens=len(t), features=t.shape[1])
    table = read(args.profile, index_col=0)
    taxonomy = read(args.taxonomy)
    require_columns(taxonomy, ["feature_id", args.rank])
    if taxonomy.feature_id.duplicated().any():
        raise ValueError("Taxonomy feature_id must be unique.")
    keys = taxonomy.set_index("feature_id")[args.rank].reindex(table.columns)
    if keys.isna().all():
        raise ValueError("No source-supported lineage at requested rank.")
    if (table < 0).any().any():
        raise ValueError("Negative abundance is invalid.")
    grouped = table.T.groupby(keys, dropna=True, sort=True).sum(min_count=1).T
    if args.semantics in ("counts", "estimated_reads"):
        save(grouped, out/(args.rank+"_"+args.semantics+".tsv.gz"), True)
        denominator = table.sum(axis=1)
        relative = grouped.div(denominator.where(denominator.gt(0)), axis=0)
    elif args.semantics == "percent":
        relative = grouped/100
    else:
        relative = grouped
    save(relative, out/(args.rank+"_relative_abundance.tsv.gz"), True)
    save(taxonomy, out/"source_taxonomy.tsv")
    save(pd.DataFrame(dict(feature_id=table.columns, supported=keys.notna().to_numpy())), out/"rank_support.tsv")
    return dict(rank=args.rank, semantics=args.semantics, specimens=len(table),
        features=grouped.shape[1], denominator="all input feature abundance, not reclosed retained subset",
        unsupported_features=int(keys.isna().sum()), ASV_or_species_inference=False)
