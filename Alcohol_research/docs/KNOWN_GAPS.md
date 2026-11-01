# Validation boundaries and known gaps

## Blocked, not validated

1. **03–05 original raw-read processing:** real recovered adapter implementations are included and syntax/preflight tested, but original reads, matching classifier/index assets and required profiling programs were not available together. A passing preflight-rejection test is not a successful read-processing run.
2. **Official FAPROTAX mapper invocation:** command and database interface are implemented; castor/optparse were absent and original optional package versions are not recovered. The accepted official mapping was reused/hash-checked, while source contributions were freshly reconstructed and checked.
3. Original raw program versions and database hashes are not established. Substituting a current classifier/database cannot be labelled exact reproduction.

## Implemented but not exhaustively rerun

- Beta: actual 9,999-permutation reproduction covers primary gut shared-7. A supplemented gut 16S node was checked and correctly stopped at 47 case / 2 reference working units with NA statistics and zero permutations. Other panel/comparison/scope combinations are not certified by these checks.
- PAM: primary ALD and Alcohol gut shared-7 were actually replayed. The supplemented ALD node correctly stopped at 18 case units with zero bootstrap attempts. Other nodes remain NOT_RUN in this release.
- Module 06: frozen genus export was exercised; optional supported-rank aggregation is implemented, but a broad rank-by-source raw validation was not performed.
- Module 14 full 406-node replay was performed. The independent single-variable branch and module 17 sex/function branch were also exercised; those partial branches withhold incomplete joint-family q-values.
- Source robustness validation covers primary gut taxa and functions, not every disease/panel/scope.
- The 22 interfaces expose frozen nodes, not arbitrary new cohorts or an open-ended analysis framework. A missing original node is an error or unavailable comparison, not permission to invent one.

## Reused inputs versus fresh calculations

- Frozen metadata, comparison labels, feature families and result authorities are not reselected.
- Module 02 reconstructs clinical values from accepted evidence but retains technical metadata.
- Module 07 without `V3` reuses canonical mapping/TSS and freshly checks CLR. With `V3`, source integration is recomputed.
- Module 08 reuses frozen labels and recomputes membership counts.
- Module 09 without `V3` reuses full-source alpha measures, then refits M1. With `V3`, measures are also recomputed and compared.
- Modules 20/21 consume frozen completed statistics; overlaps/plots themselves are fresh.
- Module 22 retains immutable narrative/runtime/assets and rebuilds HTML. No standalone PPTX or PDF was generated or claimed.
- The generic successful execution envelope does not override these component-level distinctions.

## Scope and packaging limits

The full original archives, raw reads, database trees, large source objects, clinical support files and immutable Chinese template are external. The code ZIP intentionally does not duplicate them. Acquisition URLs are provenance, not proof of local reads. Additional recovered historical routes are indexed rather than all ported. MICOM's original small pilot is outside this package's active 22-module scope and has not been expanded.

Python dependency pins record the validation environment. R and shell records are explicit about missing pins; they are not a fully locked historical raw environment. Native Windows can run the Python profile workflow with external inputs; full R/Bash adapters require appropriate environments. Short ZIP paths cannot guarantee Windows compatibility of arbitrary user-selected output roots.

## Resolved validation issues

During development, checks exposed three implementation mismatches: canonical-pooled rather than source-composition alpha, incorrect fitted feature scope for disease/original versus supplemented subset tables, and duplicate overlap-export arguments. These were corrected against the originals and the affected checks rerun. Frozen scientific results were not changed. Final receipts distinguish completed successful checks from unavailable branches; earlier development failures are not represented as successful executions.

## GitHub preparation scope

The original full-data demo, acquisition/source-object inventories, and detailed per-node execution logs remain in the supplied Biomni archive. Compact original summaries are under validation/biomni_summary/. Current checks are recorded separately in validation/repository_checks.json. The repository contains no supervised Alcohol/AUD/ALD multiclass model.
