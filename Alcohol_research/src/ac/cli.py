"""Independent module interfaces, shared output protection and dependency selection."""
import argparse
import importlib
import sys
from .inputs import *
from .io import execute

NAMES = {
    "01": "audit_inputs", "02": "collect_harmonize_metadata", "03": "preprocess_fastq",
    "04": "profile_16s", "05": "profile_shotgun", "06": "build_taxonomic_abundance",
    "07": "harmonize_profiles", "08": "define_comparisons", "09": "alpha_diversity",
    "10": "beta_diversity", "11": "alcohol_vs_health_M1", "12": "disease_comparisons_M1",
    "13": "clinical_subset_M1", "14": "clinical_interactions", "15": "source_robustness",
    "16": "functional_annotation", "17": "functional_comparisons",
    "18": "taxon_function_contributions", "19": "patient_patterns", "20": "candidate_overlap",
    "21": "make_figures", "22": "build_presentation"}

def parser(module):
    ap = argparse.ArgumentParser(description=NAMES[module].replace("_", " ")+
        ". Preserve the completed frozen analysis; see docs/MODULE_GUIDE.md.")
    ap.add_argument("--input-root", required=True, help="Extracted authoritative M1 package root, or external original verification-input root.")
    ap.add_argument("--clinical-root", help="Extracted authoritative clinical increment root.")
    ap.add_argument("--output-root", required=True, help="Separate writable output root, never an authoritative input root.")
    ap.add_argument("--resume", action="store_true", help="Reuse only matching dependency/code hashes and intact outputs.")
    ap.add_argument("--overwrite", action="store_true", help="Replace this module's output node only.")
    ap.add_argument("--dry-run", action="store_true", help="List checked inputs and target; do not execute.")
    if module in {"06", "07", "09", "10", "11", "12", "13", "14", "15", "17", "19", "20", "21"}:
        ap.add_argument("--panel", default="gut__common_all" if module in {"10", "19"} else "gut__full_source_reference",
            help="Exact frozen panel identifier. Clinical --panel all selects every relevant panel.")
    if module in {"09", "10", "11", "12", "13", "15", "17", "18", "19"}:
        ap.add_argument("--comparison", choices=COMPARISONS,
            default="04_ALD_vs_Health" if module == "19" else COMPARISONS[0])
        ap.add_argument("--scope", choices=["M1_all", "M1_clinical_complete", "supplemented"],
            default="supplemented" if module == "13" else "M1_all")
    if module == "01":
        ap.add_argument("--archives", nargs="*", help="Optional completed ZIP files for CRC and identity checks.")
    if module in {"07", "09"}:
        ap.add_argument("--source-root", help="Original v3 root with provenance/source_inputs and canonical crosswalk.")
    if module in {"03", "04", "05"}:
        ap.add_argument("--source-config", help="JSON with recovered adapter, project, original inputs and isolated output roots.")
    if module == "06":
        ap.add_argument("--profile", help="Optional specimen-by-feature TSV; omit for frozen genus TSS.")
        ap.add_argument("--taxonomy", help="Required with --profile: TSV with feature_id and source-supported rank columns.")
        ap.add_argument("--rank", choices=["phylum", "class", "order", "family", "genus", "species"], default="genus")
        ap.add_argument("--semantics", choices=["counts", "estimated_reads", "fraction", "percent"], default="fraction")
    if module == "14":
        ap.add_argument("--variable", default="all", help="One exact frozen variable key, or all.")
        ap.add_argument("--body", choices=["gut", "oral", "all"], default="all")
        ap.add_argument("--scope", choices=["frozen", "exclude_explicit_disease_references", "all"], default="all")
        ap.add_argument("--outcome", choices=["taxa", "alpha", "function", "all"], default="all")
        ap.add_argument("--coverage", choices=["base", "global_coverage", "trait_coverage", "all"], default="all")
    if module == "15":
        ap.add_argument("--outcome", choices=["taxa", "function"], default="taxa")
    if module == "16":
        ap.add_argument("--source-root", required=True, help="Original genus/function v2 root with inputs_v1 and functional/faprotax.")
        ap.add_argument("--database-root", help="Original FAPROTAX-mapper_v2.0 directory.")
        ap.add_argument("--rerun-mapper", action="store_true", help="Execute the official mapper before recomputing source contributions.")
    if module == "17":
        ap.add_argument("--coverage", choices=["base", "global_coverage", "trait_coverage"], default="base")
        ap.add_argument("--variable", help="Optional clinical variable; omit for group M1.")
        ap.add_argument("--body", choices=["gut", "oral"], default="gut")
        ap.add_argument("--clinical-scope", choices=["frozen", "exclude_explicit_disease_references"], default="frozen")
    return ap

def dependencies(module, args):
    data = Inputs(args.input_root)
    base = [data.i/"metadata_4603.tsv.gz", data.i/"comparison_membership_4603.tsv.gz"]
    if module in {"09", "11", "12", "13", "17"}:
        from .groups import model_inputs
        files = model_inputs(args, "alpha" if module == "09" else "function" if module == "17" else "taxa")
        if module == "17" and args.variable:
            c = Path(args.clinical_root)
            files += [c/"reproducibility/inputs/metadata.tsv.gz", c/"harmonised_clinical_values.tsv.gz"]
            files += [p for p in (data.i/"functions"/args.body).glob("*") if p.is_file()]
        if module == "09" and args.source_root:
            files += sorted((Path(args.source_root)/"provenance/source_inputs").glob("*/TSS.tsv.gz"))
        return files
    if module in {"01", "06", "07", "08", "10", "15", "19", "20", "21"}:
        folders = list((data.i/"panels").glob("*")) if module == "01" else [data.i/"panels"/getattr(args, "panel", "gut__full_source_reference")]
        for folder in folders:
            base += sorted(folder.glob("*"))
    if module == "01":
        base += [Path(x) for x in (args.archives or [])]
    if module == "06" and args.profile:
        base += [Path(args.profile), Path(args.taxonomy)]
    if module in {"03", "04", "05"} and args.source_config:
        base += [Path(args.source_config)]
    if module == "14":
        c = Path(args.clinical_root)
        base += [c/"reproducibility/inputs/metadata.tsv.gz", c/"harmonised_clinical_values.tsv.gz"]
        base += [x for folder in ["panels", "diversity", "functions"]
            for x in (data.i/folder).rglob("*") if x.is_file()]
    if module == "02":
        c = Path(args.clinical_root)
        base += [c/"evidence/individual_value_evidence.tsv.gz", c/"harmonised_clinical_values.tsv.gz",
            c/"reproducibility/inputs/metadata.tsv.gz", c/"evidence/unresolved_conflicts.tsv"]
    if module == "08" and args.clinical_root:
        c = Path(args.clinical_root)
        base += [c/"harmonised_clinical_values.tsv.gz", c/"reproducibility/inputs/metadata.tsv.gz"]
    if module == "07" and args.source_root:
        source = Path(args.source_root)
        base += [source/"provenance/taxon_crosswalk.tsv.gz"]
        base += [p for pattern in ("*/TSS.tsv.gz", "*/CLR_unadjusted.tsv.gz")
            for p in (source/"provenance/source_inputs").glob(pattern)]
    if module == "16":
        from .functions import annotation_inputs
        base += annotation_inputs(args)
    if module == "18":
        base += [data.i/"candidate_taxon_trait_mapping.tsv.gz",
            data.i/"candidate_taxon_trait_contribution_all_supported_samples.tsv.gz"]
        base += sorted((data.i/"candidate_function_original_TSS").glob("*"))
    if module in {"10", "19"}:
        folder = data.root/args.comparison/args.panel.split("__")[0]/args.scope/(
            "diversity" if module == "10" else "patterns")/args.panel
        base += [p for p in folder.rglob("*") if p.is_file() and "__pycache__" not in p.parts]
        if args.scope == "supplemented":
            base += [data.root/"evidence_supplemented/metadata_evidence_supplemented.tsv.gz"]
            original = data.root/"evidence_supplemented/analyses"/args.comparison/args.panel.split("__")[0]/"M1_clinical_complete_evidence_supplemented"
            base += [p for p in (original/("diversity" if module == "10" else "patterns")/args.panel).rglob("*") if p.is_file()]
    if module in {"15", "20", "21"}:
        base += sorted(data.root.glob("0*/*/M1*/taxa/*/all_results.tsv.gz"))
        base += sorted(data.root.glob("0*/*/M1*/functions/*.tsv.gz"))
        if module == "20":
            base += sorted((data.root/"evidence_supplemented/analyses").glob("0*/*/M1*/taxa/*/all_results.tsv.gz"))
    if module in {"20", "21"} and args.clinical_root:
        c = Path(args.clinical_root)
        base += [c/"interaction_models/all_interaction_tests.tsv.gz", c/"interaction_models/all_model_coefficients.tsv.gz"]
    if module == "21" and (data.root/"beta_completed_summary.tsv").exists():
        base += [data.root/"beta_completed_summary.tsv"]
        base += sorted(data.root.glob("01_Alcohol_vs_Health/gut/M1_all/diversity/*/alpha_M1.tsv"))
        base += sorted(data.root.glob("04_ALD_vs_Health/gut/M1_all/patterns/gut__common_all/PAM_completion/case_PCA_scores.tsv.gz"))
        base += sorted(data.root.glob("04_ALD_vs_Health/gut/M1_all/patterns/gut__common_all/PAM_completion/original_K2.tsv"))
    if module == "22":
        c = Path(args.clinical_root)
        base += [c/"Alcohol_Microbiome_Main_Results.deck.html"]
        base += [x for x in (c/"presentation/source").rglob("*") if x.is_file()]
        base += list((c/"interaction_figures").glob("*.data.tsv.gz"))
    return [p for p in base if not p.is_dir()]

def node_name(module, args):
    fields = [module]
    for key in ("comparison", "scope", "panel", "variable", "body", "outcome", "coverage", "rank"):
        if module == "17" and key == "body" and not args.variable:
            continue
        value = getattr(args, key, None)
        if value:
            fields.append(value)
    return "__".join(fields)

def invoke(module, args, out):
    if module in {"01", "02", "07", "08"}:
        from . import early
        return {"01": early.audit, "02": early.harmonize_metadata,
            "07": early.profiles, "08": early.comparisons}[module](args, out)
    if module in {"03", "04", "05"}:
        from .raw import run
        return run(args, out, module)
    if module == "06":
        from .abundance import build
        return build(args, out)
    if module in {"09", "11", "12", "13", "17"}:
        from .groups import run_m1, alpha
        if module == "17" and args.variable:
            from .clinical import run
            from types import SimpleNamespace
            run(SimpleNamespace(input_root=args.input_root, clinical_root=args.clinical_root,
                pilot=False, variable=args.variable, body=args.body, scope=args.clinical_scope,
                panel=args.body, outcome="function", coverage=args.coverage,
                output_prefix=str(out.resolve())))
            return json.loads((out/"run_receipt.json").read_text())
        return alpha(args, out) if module == "09" else run_m1(args, out, "function" if module == "17" else "taxa")
    if module in {"10", "15", "19"}:
        os.environ["AC_INPUT_ROOT"] = str(Path(args.input_root).resolve())
        os.environ["AC_OUTPUT_ROOT"] = str(out.resolve())
        if args.scope == "supplemented":
            if module == "15":
                raise ValueError("No supplemented source-robustness node was completed in the authoritative source module.")
            from . import supplemented
            return supplemented.beta(args, out) if module == "10" else supplemented.patterns(args, out)
        if module == "10":
            from .beta import run_node
            return dict(results=run_node(args.comparison, args.scope, args.panel))
        if module == "19":
            from .patterns import run_node
            return run_node(args.comparison, args.scope, args.panel, {})
        from .sources import nodes, one
        matches = [x for x in nodes() if (x[0], x[2], x[4], x[3]) ==
            (args.comparison, args.scope, args.panel if args.outcome == "taxa" else args.panel.split("__")[0],
             "taxa" if args.outcome == "taxa" else "functions")]
        if not matches:
            raise ValueError("No frozen source robustness node matched.")
        result = one(*matches[0])
        return dict(features=len(result), robust_primary=int((result.source_robust&result.ordinary_primary_hit).sum()))
    if module == "14":
        from .clinical import run
        from types import SimpleNamespace
        values = {k: (None if v == "all" else v) for k, v in vars(args).items()}
        values.update(pilot=False, output_prefix=str(out.resolve()))
        run(SimpleNamespace(**values))
        return json.loads((out/"run_receipt.json").read_text())
    if module in {"16", "18"}:
        from .functions import annotation, contributions
        return annotation(args, out) if module == "16" else contributions(args, out)
    if module == "20":
        from .overlap import candidate_overlap
        return candidate_overlap(args, out)
    if module == "21":
        from .figures import make
        return make(args, out)
    if module == "22":
        from .presentation import build
        return build(args, out)
    raise ValueError("Unknown module.")

def entry(module):
    args = parser(module).parse_args()
    if module in {"02", "14", "22"} and not args.clinical_root:
        parser(module).error("--clinical-root is required for this module.")
    if module == "17" and args.variable and not args.clinical_root:
        parser(module).error("--clinical-root is required for clinical functional models.")
    if module == "11" and args.comparison != COMPARISONS[0]:
        parser(module).error("Module 11 is only the full Alcohol vs Health analysis; use module 12.")
    if module == "11" and args.scope != "M1_all":
        parser(module).error("Use module 13 for clinical-complete subsets.")
    if module == "12" and args.comparison == COMPARISONS[0]:
        parser(module).error("Select one disease comparison explicitly.")
    if module == "06" and args.profile and not args.taxonomy:
        parser(module).error("--taxonomy is required with --profile.")
    try:
        execute(args, node_name(module, args), dependencies(module, args),
            lambda out: invoke(module, args, out))
    except Exception as error:
        print(type(error).__name__+": "+str(error), file=sys.stderr)
        raise SystemExit(2)
