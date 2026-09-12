# Published-result reproducibility audit

This directory keeps the published targets separate from checked-in result
snapshots and future reruns.

## Files

- `published_reference.json` transcribes the quantitative values printed in the
  main article and Supporting Information. It also declares transcription
  tolerances equal to one half-unit in the final displayed digit.
- `result_manifest.json` maps every main-text and supporting figure/table to its
  notebook, archived data, and required comparison method.
- `phase_a_coverage.json` gives the final evidence-backed disposition of all 22
  published figures and tables and states the Phase A PR gate limitations.
- `audit_published_tables.py` compares Tables S4-S6 with the two checked-in CSV
  families.
- `audit_source_tables.py` compares the model inputs in Tables S1-S3 with the
  public Python, notebook, and CSV sources.
- `runs/` contains immutable structured records for executed reproduction
  attempts, including failed or interrupted attempts that affect provenance.

The transcription tolerances are not solver regression tolerances. Phase A will
add tighter, metric-specific solver tolerances after repeat-run variability has
been measured with the reconstructed environment.

## Run the current audit

From the repository root:

```bash
python reproducibility/audit_published_tables.py
python reproducibility/audit_published_tables.py --snapshot migrated
python reproducibility/audit_source_tables.py
python -m unittest discover -s tests -v
```

After creating a historical candidate environment, rerun the M5/Bakken tax
sequence without modifying archived checkpoints or CSVs:

```bash
python reproducibility/run_m5_bakken_tax_series.py \
  --ipopt /absolute/path/to/hsl-enabled/ipopt \
  --linear-solver ma27 \
  --tee \
  --output /path/outside/the/repository/m5-bakken-tax-series.json
```

Tax rates are expressed in USD/kg CO2e, matching the optimization notebook;
the corresponding reported table values are 1,000 times larger in USD/tonne.
To isolate a difficult sequential transition, `--initial-optimal-tax 0.19
--tax-rates 0.41` loads the archived USD 190/tonne optimum before solving the
USD 410/tonne case.

Regenerate all 23 archived optimal composite-curve PDFs into a separate
directory, together with a machine-readable record of their numeric
coordinates and source-workbook hash:

```bash
python reproducibility/regenerate_composite_curves.py \
  --output-dir /path/outside/the/repository/composite-curves \
  --record /path/outside/the/repository/composite-curves.json
```

This covers Main Figure 4 and Supporting Figures S3, S4, and S6-S8. The command
reads `results/solution_data.xlsx` and the checked-in postprocessed result CSVs,
uses the original composite-curve algorithm, and never overwrites the archived
publication plots.

Regenerate the other eight quantitative figures from the immutable optimal
checkpoints and checked-in CSV inputs:

```bash
python reproducibility/regenerate_archived_figures.py \
  --output-dir /path/outside/the/repository/archived-figures \
  --record /path/outside/the/repository/archived-figures.json
```

This covers Main Figures 3 and 5-8 and Supporting Figures S1, S2, and S5.
The driver builds each ROK topology once, loads the archived optimal states,
and reuses the M5 topology for all regional checkpoints. The external JSON
record preserves every component series and SHA-256 source hash. Published
emissions-label differences remain classified under
`EMISSIONS-NORMALIZATION-001`; the EF-7 and EF-10 Figure 7 labels retain the
already documented publication-snapshot drift.

Rerun one M5 regional case with the initialization and temperature
perturbation used by its published notebook:

```bash
python reproducibility/run_m5_region.py \
  --region EF-1 \
  --ipopt /absolute/path/to/hsl-enabled/ipopt \
  --linear-solver ma27 \
  --tee \
  --output /path/outside/the/repository/m5-ef-1.json
```

Use `run_m5_region_series.py --regions EF-1 EF-2 ... EF-12` for an incremental
multi-zone run. The coordinator reloads the archived EF-Basin optimum before
each zone, except EF-9, whose published notebook reloads EF-8. It writes the
JSON record after every completed zone so an interrupted run retains completed
evidence. The single-case runner also accepts `--initial-optimum` as a
compatibility control when the published initialization path is solver
sensitive; this is a diagnostic restart and not a substitute for reproducing
the notebook sequence. Completed regional records include the fresh stream
table, heat-exchanger table, composite-curve coordinates, liquid-product
component flows and LHV contributions, and upstream emissions needed to audit
Figures 6-8 and S6-S8 as well as Table S6.

Audit either regional pre-solve state without invoking IPOPT:

```bash
python reproducibility/diagnose_m5_region_state.py \
  --region EF-11 \
  --output /path/outside/the/repository/m5-ef-11-state-diagnostics.json
```

Add `--initial-optimum` to inspect the archived target-region checkpoint. The
report uses the diagnostics available in the historical IDAES commit and lists
the largest constraint residuals, badly scaled free variables, variables near
bounds, and variables or constraints without scaling factors. This is intended
to distinguish a poor regional-substitution state from a feasible checkpoint
that is no longer stationary under the reconstructed solver stack.

Compare a completed regional run record with the archived stream and
heat-integration sheets plus the numeric labels transcribed from Figures 6 and
7:

```bash
python reproducibility/compare_regional_figure_data.py \
  --run-record /path/outside/the/repository/m5-regions.json \
  --output /path/outside/the/repository/m5-regions-comparison.json \
  --quiet
```

The comparison reports numerical table differences separately from unit-label
text differences because the historical IDAES table API serializes unit objects
differently from the labels already stored in the workbook.

The curated regional evidence is in
`runs/2026-09-12-candidate-a1-m5-regional-ma27.json`. Nine of the thirteen
EF-Basin/EF-1-through-EF-12 cases closely reproduce the archived scalar state on
the exact notebook path. EF-8 reaches a lower alternate local optimum, EF-10
has a smaller solver-path discrepancy, and EF-2 and EF-11 do not yield accepted
exact-path solutions with either MA27 or MA57. Checkpoint restarts are recorded
as controls and are not counted as reproductions. The record also preserves the
maximum detailed stream-table and composite-curve differences so matching a
single rounded figure label cannot mask a different stationary point.

## Phase B baseline

The IDAES 2.12 modernization branch starts from the Phase A head. Its first
baseline is recorded as `B-IDAES-2.12-BASELINE-001`: all repository modules
import, M5/Bakken builds, the historical unit/constrained/costed checkpoint
chain loads, all 24 tests pass under both the historical and modern
environments, and MA27/MA57 smoke solves terminate optimally. The modern build
contains 156 fewer component data objects than Candidate A1, so later numerical
comparisons must not assume internal object identity even when named model
states agree.

The default `postprocessed` snapshot uses the root-level CSV files created by
commit `957e363`, which recalculated TAC and MSP using cooling water above 303 K
and refrigerated water from 288-303 K. The `migrated` snapshot uses the CSV files
under `results/` that arrived in the initial public migration commit `2295e03`.
Neither snapshot is silently treated as the paper.

## Known provenance findings

The published process/downstream emissions in Tables S4-S6 and Figures 3 and 6
do not match the CSV values. For example, Table S4 reports 9.17 g CO2e/MJ for M2,
while both checked-in CSV snapshots report 8.876.

The public implementation constrains both upstream and downstream emissions
using this denominator in `src/emissions_calculations.py`:

```text
fuel_energy + gas_energy + h2_energy
```

The paper describes emissions per unit of liquid fuel, and the printed values
appear to use an earlier normalization convention. The public repository starts
with a bulk migration from a private research repository, so its Git history
does not contain the equation version that generated the paper numbers. The
audit labels this unresolved issue `EMISSIONS-NORMALIZATION-001` and continues
to fail it under `--strict`.

The audit also labels 22 postprocessed decision/cost values
`PUBLISHED-SNAPSHOT-DRIFT-001`. They exceed literal print-rounding tolerance,
which establishes that the SI tables and the checked-in higher-precision CSVs
were not generated from precisely the same saved optimization snapshot. The
differences are reported rather than erased by widening tolerances.

The Tables S1-S3 audit records two additional provenance findings. Table S1
prints a 295-500 K bound for H104, but the optimization update never applies
those bounds; a reconstructed A1 model retains the property-package range
273.15-1500 K. Table S2's cooling-water and refrigerated-water prices are split
between the optimization source and the publication postprocessing added in
commit `957e363`. The printed methane-recovery fraction and methane GWP do not
appear as semantic inputs anywhere in the public computational source or its
history back to the initial private-repository migration. They may have been
used upstream to prepare the regional emissions factors, but that derivation
was not published in this repository.

Finally, the migrated snapshot has 29 additional MSP/TAC differences labeled
`LEGACY-ECONOMICS-POSTPROCESSING-001`. Commit `957e363` explains these: it
recalculated economics using separate cooling-water and refrigerated-water
grades, producing the root-level `postprocessed` CSV family.

## Evidence sources

- Ghosh et al. (2024), DOI `10.1021/acssuschemeng.4c00933`.
- Main article Figures 1-8, especially the numeric labels in Figures 3, 6, and 7.
- Supporting Information Figures S1-S8 and Tables S1-S6.
- `results/solution_data.xlsx`, containing 23 stream sheets and 46 initialized/
  optimal heat-integration sheets plus a README sheet.
- Saved notebook metadata and IPOPT output.

The source PDFs are not redistributed in this repository.
