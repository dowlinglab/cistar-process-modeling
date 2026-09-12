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

To test whether two IDAES/Pyomo environments construct the same nonlinear NLP
away from an archived point, generate a detailed derivative fingerprint in
each environment and compare the resulting files by component name:

```bash
python reproducibility/fingerprint_m5_bakken_derivatives.py \
  --perturb H105_temperature=1.0 \
  --output /path/outside/the/repository/fingerprint.json \
  --detail-output /path/outside/the/repository/fingerprint.npz

python reproducibility/compare_derivative_fingerprints.py \
  /path/to/reference.npz /path/to/candidate.npz \
  --output /path/outside/the/repository/comparison.json
```

The comparison aligns all variables and constraints by name and canonicalizes
the Hessian triangle before reporting sparsity and scale-aware differences.
Record `B-M5-BAKKEN-PERTURBED-DERIVATIVE-COMPARISON-001` finds matching named
NLPs and derivative sparsity under Candidate A1 and IDAES 2.12. At the archived
tax-zero M5/Bakken point with H105 displaced by 1 K, no Jacobian or Lagrangian
Hessian entry differs by more than `1e-10` on the recorded scale-aware metric.
Together with the checkpoint stationarity control, this makes a changed model
equation or higher derivative an unsupported explanation for the modern solve
failure. Writer/solver path sensitivity amplified by the severely conditioned
state Jacobian remains the leading diagnosis.

A full 100-iteration `nl_v1`/MA27 control is recorded as
`B-M5-BAKKEN-NLV1-FULL-001`. It terminates at `maxIterations` after moving
materially away from the archived solution, so selecting the legacy-compatible
writer alone is not a modernization fix. The writer also reports 440 inherited
scaling-suffix keys for components absent from the exported NL and 1,750 keys
whose component types cannot be exported; this points to explicit active-NLP
scaling, rather than writer selection, as the next intervention.

That intervention is available as the diagnostic-only
`--active-nlp-autoscale` flag and recorded in
`B-M5-BAKKEN-ACTIVE-NLP-SCALING-SCREEN-001`. It clears 7,638 inherited entries
and rebuilds factors for exactly 5,555 active variables and 5,618 active
constraints, eliminating all skipped-key writer warnings. The numerical result
is substantially worse: the first MA57 step raises the unscaled constraint
violation to `9.15e4`, and dual infeasibility remains `2.49e6` at iteration ten.
Thus skipped suffix keys are not the cause, and this active-row policy is not a
recommended scaling fix. The difference from the more stable whole-model
AutoScaler control narrows the next comparison to how inequalities and
inherited constraint factors are treated.

The subsequent apples-to-apples and factor-level controls are recorded in
`B-M5-BAKKEN-SCALING-FACTOR-COMPARISON-001` and supersede that preliminary
narrowing. Matching IPOPT `user-scaling` and preserving all 71 inherited
inequality factors still fail to reproduce the whole-model AutoScaler path.
All 5,555 active variable factors match, but 5,302 of 5,547 equality factors
differ between the policies, including 1,974 by more than a factor of ten. A
three-constraint reproducer identifies an IDAES 2.12 row-index mismatch: the
block-level AutoScaler enumerates equality constraints while indexing row norms
from the full Jacobian, so an interspersed inequality shifts factors assigned
to later equalities. The more stable whole-model trajectory therefore relies on
accidentally misaligned equality factors and is diagnostic evidence, not a
defensible modernization policy.

The optimization runner also supports `--free-design-variables` for
unit/decision-level gates and `--file-determinism sort-symbols` for a
cross-version ordering control. Record
`B-M5-BAKKEN-DOF-ORDERING-CONTROLS-001` shows that the default-order modern
model fails even with only H103/R102 temperature free, while symbol ordering
converts that one-DOF case to an accepted optimum. The complete eight-DOF
modern case still fails with symbol ordering, whereas Candidate A1 converges
and closely reproduces the archived result under the same requested rule.
Ordering is therefore an important sensitivity but the high-level writer
setting is not a complete explanation; actual emitted NL symbol maps and
writer grouping must be compared across Pyomo versions.

That direct comparison is recorded in
`B-M5-BAKKEN-NL-SYMBOL-MAP-COMPARISON-001`. The named component sets are
identical and all 5,618 constraints retain the same row position, but only 780
of 5,555 variables retain the same column position; the common prefix is just
11 variables and the largest displacement is 2,460 columns in an R102 outlet
state. Pyomo 6.10.1 also emits 548 reusable and 7,309 single-use constraint
common expressions where Pyomo 6.4.2 emits none, reducing the NL file from
661 MB to 167 MB with unchanged dimensions and Jacobian nonzero count. This
isolates a real writer representation change capable of altering IPOPT's path,
while the matching named equations and derivative fingerprints continue to
rule against a changed mathematical model as the explanation.

The causal follow-up is recorded in
`B-M5-BAKKEN-HISTORICAL-COLUMN-ORDER-001`. The tax-series runner accepts a
digest-verified `--column-order-from-symbol-map` diagnostic. When the modern
writer is given Candidate A1's exported variable sequence, its emitted map
matches all 5,555 historical column positions (and all rows already matched)
while retaining the modern common-expression representation. The full modern
eight-DOF tax-zero case then terminates optimally in 46 iterations and matches
the migrated result to `1.1e-10` USD/MJ in MSP, `5.6e-4` K in R102 temperature,
and 23 USD/year in TAC. This establishes column permutation as causal for this
failure, not merely correlated with it. The pinned sequence is a defensible
cross-version reproduction control; representation-robust scaling remains the
preferred long-term modernization outcome.

The first tax-continuation screen is recorded in
`B-M5-BAKKEN-HISTORICAL-COLUMN-TAX-CONTINUATIONS-001`. Continuing directly
from the fresh modern tax-zero solution reaches an accepted but incorrect local
optimum at USD 0.01/tonne, despite the tax-zero headline metrics matching the
archive. Reloading the exact archived tax-zero checkpoint instead closely
recovers USD 0.01/tonne. The next transition, to USD 1/tonne, fails from its
exact preceding archived checkpoint under MA27, MA57, and MA57's internal
automatic scaling. MA57 without internal scaling ends closest to feasibility
(`2.32e-2`) but remains in restoration. The `nl_v1` writer reproduces the
failing MA27 trace, so writer selection is not a remedy for this transition.
This is a partial modern reproduction: checkpoint state and linear-system
regularization remain consequential even after column order is pinned.

Complete named-state and pseudo-zero-phase follow-up is recorded in
`B-M5-BAKKEN-STATE-DRIFT-PSEUDO-ZERO-REGULARIZATION-001`. Restricting the
checkpoint comparison to the 5,555 active NL variables isolates 66 dominant
differences in four structurally absent inlet phases. Their component flows are
near `1e-8` mol/s, so the associated composition-normalization equations are
too weak to prevent large composition drift at the global feasibility
tolerance. An opt-in guarded regularization removes those 66 equations and
variables, but a direct USD 1/tonne solve then converges to the same incorrect
low-temperature local optimum. Even a USD 0.02/tonne step crosses a nearly
feasible, nearly stationary iterate before restoration failure. Explicitly
bounded acceptable termination retains the first two small continuation
points but fails by USD 0.10/tonne; it is diagnostic evidence, not the Phase B
reproduction policy. The runner now stops a continuation after a failed solve
by default so later cases cannot inherit an invalid state.

The decisive writer control is recorded in
`B-M5-BAKKEN-DEFINED-VARIABLE-INLINING-001`. Pyomo 6.10's NL writer exposes an
`export_defined_variables` option that defaults to true; Pyomo 6.4.2 emitted no
common expressions for this model. With the historical variable-order request
already supplied, `--inline-defined-variables` sets that option false. The
modern NL keeps the same 5,555 variables, 5,618 constraints, component sets,
constraint ordering, and Jacobian sparsity, but its reported
Lagrangian-Hessian nonzeros fall from 26,211 to 17,776. Inlining changes
nonlinear-variable classification from 2,305 to 2,317 columns, so the emitted
variable order has a 2,305-column historical prefix and 2,573 exact positions
rather than being completely identical. Starting from the archived tax-zero
checkpoint, the complete ascending
USD 0.01, 1, 17, 45, 190, and 410/tonne sequence then reaches strict optima in
12-14 iterations per point. The largest absolute discrepancy from the migrated
CSV is `1.20e-9` USD/MJ in MSP, `4.83e-4` K in R102 temperature, and 39.01
USD/year in TAC. Even USD 410/tonne converges directly from the newly solved
USD 190/tonne point, without an archived restart.

A complete control then removed the historical column-order request. Inlining
alone reproduces all six points as strict optima, with maximum absolute
differences of `3.87e-10` USD/MJ in MSP, `4.89e-4` K in R102 temperature, and
41.29 USD/year in TAC. The historical symbol map is therefore unnecessary in
the recommended modern workflow; it remains useful only as a diagnostic
ordering control.

The corresponding fresh tax-zero optimization, started from the ordinary
costed initialization checkpoint rather than an archived optimum, also reaches
a strict optimum with inlining alone in 20 iterations. Its MSP, R102
temperature, and TAC differences are `8.82e-10` USD/MJ, `7.77e-4` K, and
40.15 USD/year. Thus the same compatibility option covers both the base solve
and the full positive-tax continuation.

The same inlining-only policy reproduces the M2-M4 Bakken model-comparison
cases in record `B-ROK-MODEL-COMPARISON-DEFINED-VARIABLE-INLINING-001`. All
three terminate at strict MA27 optima in 17-19 iterations. Across the three
cases, the largest absolute differences from the migrated CSV are
`1.28e-10` USD/MJ in MSP, `4.00e-4` K in R102 temperature, and 1.22 USD/year
in TAC. No historical column map or model-specific scaling change is used.

The same record preserves the negative controls that led to this result. The
archived and fresh tax-zero R102 scaled row and column norms differ by at most
about six parts per million, so checkpoint drift does not materially reshape
the local reactor Jacobian. R102-only heat-variable rescaling and row
normalization either fail badly or converge to an incorrect local optimum.
Those controls remain available for diagnostics; they are not part of the
modern reproduction policy. For this model, defined-variable inlining is the
validated compatibility policy.

The regional extension is recorded in
`B-M5-REGIONAL-DEFINED-VARIABLE-INLINING-001`. Eleven of twelve M5 regional
cases reach IPOPT optimal termination using inlining alone. This includes EF-2
and EF-11, which did not produce accepted exact-path solutions under the
reconstructed historical environment, and EF-8, which had previously entered
an alternate optimum. EF-10 retains the largest small path discrepancy: an
unscaled constraint violation of `5.72e-6`, an R102-temperature difference of
`5.69e-2` K, and a TAC difference of 6,182 USD/year.

EF-9 is the sole ordering-sensitive exception. Its inlining-only MA27 path was
stopped after 19 minutes in restoration, and MA57 reached its 100-iteration
limit with unscaled constraint and dual infeasibilities of `0.821` and `33.9`.
Inlining plus the digest-verified Candidate A1 column-order request reaches an
MA27 optimum and matches the migrated result within `2.95e-9` USD/MJ in MSP,
`3.43e-3` K in R102 temperature, and 506 USD/year in TAC. The historical map
is therefore retained as a narrow EF-9 compatibility exception, not a global
modernization dependency.

Modern Pyomo also required one lifecycle correction before regional switching:
generated index sets are not always registered as named block components.
Region cleanup now resolves each component by name and safely skips absent
generated components. Remaining figure-generation and all-paper audit gates
must still pass before Phase B is complete.

The fresh regional figure-input comparison is recorded in
`B-M5-REGIONAL-FRESH-FIGURE-COMPARISON-001`. All twelve accepted modern states
yield stream and heat-exchanger tables, component product series, emissions,
and composite-curve coordinates. Ten fresh Figure 7 totals match their printed
integer labels within half a unit. EF-7 retains the existing published-snapshot
difference, whereas EF-10's 0.78 MW label miss is classified separately as a
fresh-solution difference. EF-10 also has the largest detailed drift: a 0.082%
material stream-table difference and 3.76 GJ/h maximum hot-curve heat shift.
Figure 6's printed emissions retain `EMISSIONS-NORMALIZATION-001` for every
region. The comparison does not turn a reproduced archived plot into proof that
the fresh solution and printed paper are exactly identical.

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
