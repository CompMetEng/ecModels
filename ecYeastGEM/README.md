# ecYeastGEM 9.1.1

The enzyme-constrained version of [yeast-GEM](https://github.com/SysBioChalmers/yeast-GEM)
**9.1.1**, built with [geckopy](https://github.com/SysBioChalmers/geckopy)
and tuned against measured growth rates and exchange fluxes with CMA-ES.
It is derived from the `v9.1.1` release of yeast-GEM and describes that
release only: the reaction and gene identifiers of later yeast-GEM
versions are not guaranteed to match it.

The model is `models/ecYeastGEM.yml`. This folder also holds everything
needed to rebuild it from the raw inputs: one script (`build.py`), one
parameter file (`model_adapter.toml`), the pinned geckopy version
(`requirements.txt`), and the database snapshots and experimental data in
`data/`. geckopy's `tutorials/full_ecModel` walks
through the same protocol step by step; its model is a teaching example,
not for reuse.

## Using the model

```bash
git clone --branch gecko4 --depth 1 https://github.com/SysBioChalmers/ecModels.git
pip install -r ecModels/ecYeastGEM/requirements.txt
```

```python
from geckopy import ModelAdapter, load_ec_model

adapter = ModelAdapter.from_folder("ecModels/ecYeastGEM")
model = load_ec_model("ecYeastGEM.yml", adapter=adapter)  # read from ecModels/ecYeastGEM/models/
model.reactions.get_by_id(adapter.params.c_source).lower_bound = -1000  # unconstrained glucose
model.objective = adapter.params.bio_rxn
print(model.optimize().objective_value)                                  # 0.4111 /h
```

## Status

| | |
|---|---|
| Base model | yeast-GEM 9.1.1 (release `v9.1.1`): 4105 reactions, 2748 metabolites, 1143 genes |
| ecModel | 8010 reactions, 1144 enzymes (1143 gene products and the `standard` pseudo-enzyme), 4850 enzyme-constrained reaction entries |
| kcats | 4850, all positive; median 11 1/s |
| Protein budget | `prot_pool_exchange` upper bound 125 mg/gDW |
| Maximum growth, glucose unconstrained | 0.1897 /h before tuning, **0.4111 /h** after (measured 0.41) |
| Fit to the tuning data | RMSE 6.11 before tuning, **0.974** after (a second seed: 0.892) |
| Kcats tuned | 108 searched, 107 moved, 31 by more than 2-fold |

The RMSE is measured on the same 33 flux conditions and 8 maximum growth
rates the kcats were fitted to. It is an in-sample fit, not a held-out
validation.

## Layout

| path | what it is |
|------|------------|
| `model_adapter.toml` | organism parameters, kcat trust tiers, CMA-ES settings |
| `build.py` | the pipeline, one subcommand per stage |
| `models/yeast-GEM.yml` | yeast-GEM 9.1.1, fetched by `build.py` and not tracked |
| `anaerobic.py` | the anaerobic switch the tuner applies to oxygen-free conditions |
| `requirements.txt` | the geckopy version the model is built with |
| `data/` | database snapshots, curated kcats and the experimental data tuned against |
| `models/` | `ecYeastGEM.yml` is the tuned model; the yeast-GEM base and the stage checkpoints are built by `build.py` and not tracked |
| `output/` | what the build reports about itself |

## Rebuilding

```bash
pip install -r requirements.txt
python build.py all --n-proc 32
```

Run from this folder (`ecYeastGEM`). `requirements.txt` pins geckopy to
the version that reproduces the model.

Each step skips itself when its output is already on disk, so a run can
be interrupted and resumed, or driven one step at a time. The tuned model
and the reports are tracked; `screen` and `tune` are skipped only when
their run files (`output/*.pkl`) exist, so a fresh clone rebuilds
everything and overwrites them:

| step | does | output | time at 32 cores |
|------|------|--------|-----------------|
| `stage1` | fetches yeast-GEM 9.1.1, expands it into an ecModel, applies Complex Portal stoichiometry | `models/yeast-GEM.yml`, `models/ecYeastGEM_stage1.yml` | under a minute |
| `kcats` | fills EC codes, matches BRENDA, resolves metabolite SMILES | `output/kcat_list_brenda.*`, `models/ecYeastGEM_stage2_smiles.yml` | about a minute |
| `okp-submit` | submits enzyme-substrate pairs to OpenKineticsPredictor | `data/OKP.csv`, `data/OKP_job.txt` | depends on the predictor's queue |
| `okp-fetch` | merges every kcat source, finishes the model | `models/ecYeastGEM_stage2.yml`, `output/kcat_list_*.tsv` | about a minute |
| `screen` | ranks kcats by how much the measured data can speak to them | `output/kcat_leverage_screen.*` | about 55 minutes |
| `tune` | fits the highest-leverage kcats with CMA-ES | `models/ecYeastGEM.yml`, `output/corrections.tsv` and the other reports | about 40 minutes |

`tune` takes `--seed` (default 0) and `--target-impact-share` (default
0.95, see [Choosing how many kcats to tune](#choosing-how-many-kcats-to-tune)).
Only `screen` and `tune` are parallel; `--n-proc` sets their worker
count. A `gurobi.env` next to `build.py` containing `Threads 1` keeps
each worker's LP solve single-threaded so the workers do not
oversubscribe the machine.

`okp-submit` needs an OpenKineticsPredictor API key (free, generated at
<https://predictor.openkinetics.org/api-docs>) in `OKP_API_KEY` or
`data/okpApiKey.txt`; both stay out of git. The predictor's output is
stored in `data/OKP_output.csv`, and `okp-fetch` uses that file when it
is present, so a rebuild needs neither the key nor the network.
`data/PhylDist.mat` is downloaded from RAVEN on the first `kcats` run and
is likewise untracked.

## Inputs

| input | what | where |
|-------|------|-------|
| yeast-GEM 9.1.1 | the conventional model, downloaded from the `v9.1.1` release of yeast-GEM and checked against a sha256 | `models/yeast-GEM.yml` |
| UniProt | reviewed *S. cerevisiae* proteome UP000002311: sequences, molecular weights, EC codes | `data/uniprot.tsv` |
| KEGG | organism `sce`: gene links and EC codes | `data/kegg.tsv` |
| Complex Portal | complexes of taxon 559292, for subunit stoichiometry | `data/ComplexPortal.json` |
| BRENDA | release 2026.1 (generated 2026-05-28): about 39 000 kcat, 28 000 molecular-weight and 18 000 specific-activity entries, bundled with geckopy | `geckopy/data/brenda/` |
| KEGG phylogenetic distances | to pick the closest organism when BRENDA has no exact match | `data/PhylDist.mat` |
| OpenKineticsPredictor | EITLEM kcat predictions | `data/OKP.csv`, `data/OKP_output.csv` |
| Curated kcats | 92 entries with a stated source, covering 219 model kcats | `data/customKcats.tsv` |
| Metabolite lists | currency and ignored metabolites excluded from predictor input; pseudo-reactions excluded from the standard kcat | `data/DLKcat*.tsv`, `data/pseudoRxns.tsv` |
| Tuning data | see [Tuning data](#tuning-data) | `data/evotune*.tsv` |

## How the model is built

### Stage 1: ecModel structure

`make_ec_model` expands yeast-GEM into 8009 reactions. Reversible
reactions get a `_REV` copy, and a reaction catalysed by several isozymes
or complexes becomes one `_EXP_<n>` entry per alternative. Each of the
1143 genes becomes an enzyme with a molecular weight from UniProt and a
usage reaction that draws on a shared protein pool. `apply_complex_data`
then sets subunit stoichiometry from Complex Portal wherever a complex
matches the model's subunits. The result has 4679 enzyme-constrained
reaction entries.

### Stage 2: assigning kcats

**EC codes.** `fill_eccodes_from_database` takes each enzyme's EC codes
from UniProt and KEGG.

**BRENDA.** `load_brenda_data` reads the bundled snapshot and drops any
row whose kcat is more than four orders of magnitude below the median of
its own EC class (the class widens to a shorter EC prefix when fewer
than five reference values exist). Such rows are almost always errors in
a specific-activity × molecular-weight join, and a kcat that small makes
an enzyme's protein cost enormous. `fuzzy_kcat_matching` then matches each
reaction to BRENDA by progressively relaxing the EC number (wildcards
from the right), the substrate and the organism, using the KEGG
phylogenetic distance to pick the closest organism when the exact one is
absent. Each match carries an *origin* (1-6) and a wildcard level:

| origin | matched on |
|-------:|------------|
| 1 | organism, substrate, kcat |
| 2 | any organism, substrate, kcat |
| 3 | organism, any substrate, kcat |
| 4 | any organism, any substrate, kcat |
| 5 | organism, specific activity × molecular weight |
| 6 | any organism, specific activity × molecular weight |

**Predictions.** `find_met_smiles` resolves metabolite SMILES (cached in
`data/smilesDB.tsv`). `submit_open_kinetics_predictor` sends the model's
enzyme-substrate pairs (2744 of them) to OpenKineticsPredictor with the
EITLEM method; sequences beyond the predictor's 1024-residue limit are
truncated server-side. Parsing the answer gives 4204 prediction rows.

**Merging.** `merge_kcats` keeps, for each reaction, the highest-ranked
source available:

1. `database_exact`: a measured value the predictor service returned from
   BRENDA or Sabio-RK,
2. `database_top`: a fuzzy BRENDA match with no EC wildcard and origin
   ≤ 4 (`top_origin_limit=4`),
3. `uniprot`: a UniProt-curated value,
4. `EITLEM-Kinetics`: the prediction,
5. `database_bottom`: any other BRENDA match, up to origin 6 and three
   wildcards.

Origins 5 and 6 are specific-activity-derived, so they rank below a
prediction; they are used only where no prediction exists.

**Curation and gap filling**, in this order:

1. `apply_custom_kcats` overrides reactions listed in
   `data/customKcats.tsv` (three matching modes, by enzyme set and/or
   reaction list),
2. `fill_kcats_from_isozymes` gives a reaction with no kcat the mean of
   its sibling isozymes' kcats,
3. `assign_standard_kcat` adds a `standard` pseudo-enzyme (median UniProt
   molecular weight) and a standard kcat (the median of the assigned
   kcats, per subsystem where there are enough). It attaches them to
   reactions that have no gene association, and gives the standard kcat
   to enzyme-associated reactions that still have none, so every
   constrained reaction carries a protein cost.

**Constraints and budget.** `apply_kcat_constraints` writes each kcat into
the stoichiometry: reaction *r* consumes `subunits · MW / (kcat · 3600)`
of its enzyme's usage metabolite, with MW in Da and kcat in 1/s.
`set_prot_pool_size` bounds the shared protein pool at
`p_tot · f · sigma · 1000` = 0.5 · 0.5 · 0.5 · 1000 = **125 mg/gDW**,
where `p_tot` is total protein per dry weight, `f` the fraction that is
model enzymes and `sigma` their average saturation.

The finished stage 2 model grows at 0.1897 /h on unlimited glucose, with
the protein pool fully used.

### Where the 4850 kcats come from

| source (`ec.source`) | kcats | prior width in tuning |
|--------|------:|----------------------:|
| BRENDA, fuzzy-matched | 2909 | 0.2 |
| EITLEM prediction | 1464 | 0.25 |
| `data/customKcats.tsv` | 219 | 0.1 |
| standard kcat | 255 | 0.3 |
| filled from an isozyme | 3 | 0.3 |

The BRENDA-matched kcats by origin: 214, 687, 693, 635, 542 and 138 for
origins 1 to 6. Twenty-four kcats are below 10⁻³ 1/s (the lowest is
1.6×10⁻⁴) and are worth checking against the literature. Curate any of
them at the source in `data/customKcats.tsv`.

## Tuning

### Tuning data

| file | content |
|------|---------|
| `evotuneFluxData.tsv` | 33 conditions (30 on glucose, 2 on fructose, 1 on sucrose): growth rate and exchange fluxes |
| `evotuneMaxGrowth.tsv` | 8 maximum growth rates, one carbon source each, uptake unconstrained (glucose, fructose, sucrose, galactose, maltose, mannose, acetate, ethanol) |
| `evotuneZeroExch.tsv` | 219 reactions assumed to carry no flux in any condition |

Oxygen-free conditions are simulated through `anaerobic.py`.

### Screen

`screen` moves each group of tied isozyme copies up and down two-fold and
records the largest change in RMSE, for each of the 4312 groups. The
leverage is weighted by the source's prior width, so a trusted kcat needs
proportionally more measured effect to rank alongside an untrusted one.
`output/kcat_leverage_screen.tsv` is useful for curation on its own: it
lists what the measured data can and cannot speak to.

### Choosing how many kcats to tune

`select_tunable_mask` keeps the fewest highest-ranked groups whose
combined leverage reaches `target_impact_share` of the total. The share
sets how many kcats are searched, and the fit depends steeply on it
(mean of two seeds; measured with `sweep_tunable_mask_size`):

| `target_impact_share` | kcats | RMSE |
|------:|------:|-----:|
| 0.2 | 2 | 5.95 |
| 0.5 | 8 | 5.73 |
| 0.6 | 10 | 4.75 |
| 0.7 | 17 | 2.88 |
| 0.8 | 32 | 2.25 |
| 0.9 | 55 | 1.46 |
| 0.93 | 76 | 1.43 |
| **0.95** | **108** | **0.93** |
| 0.97 | 207 | 0.92 |
| 0.99 | 815 | 0.82 |

The fit improves steeply up to about 100 kcats and slowly after. 0.95
gives 108 kcats and is the default; going to 815 kcats improves the RMSE
by 0.1 and moves seven times as many kcats. The adapter's
`max_generations = 212` corresponds to a search over about 110 kcats.

### CMA-ES

`cmaes_kcat_tuning` searches the selected kcats in log space, starting
from their assigned values. The objective is the fit,
`(RMSE_flux + w · RMSE_maxgrowth) / (w + 1)` with `w = 2` so the eight
maximum-growth conditions count double against the 33 flux conditions,
plus a prior term `0.03 · mean((log(k/k0) / sigma0_log)²)` that charges
for moving a kcat away from its assigned value `k0`. `sigma0_log` is the
source's prior width (custom 0.1, database 0.2, EITLEM 0.25, other 0.3),
so a trusted source is charged more for the same move. The RMSE quoted
throughout is the fit term alone. The search runs 212 generations without
reaching the 0.2 RMSE threshold, and `models/ecYeastGEM.yml` holds the
best kcat vector it found.

### Result

The RMSE falls from 6.11 to 0.974 (seed 0). Of the 108 searched kcats,
107 changed, with a median 1.42-fold move; 31 moved by more than 2-fold
and 8 by more than 5-fold. Median fold change over the searched kcats, by source: database
1.25, EITLEM 1.78, custom 1.50, unlabelled 1.34. The order does not
follow the prior widths: the custom kcats, which have the narrowest
prior, moved further than the database ones. The custom kcats that moved
by more than 2-fold are on ATP synthase, phosphofructokinase,
fructose-bisphosphate aldolase and glyceraldehyde-3-phosphate
dehydrogenase.

Leverage is concentrated: the top 5 corrections carry 46% of it, the top
10 67%, the top 25 81% and the top 50 95%. Four of the top five are
isozyme copies of `r_1166`.

## Reading the tuning result

`output/` holds the artifacts worth reviewing, rather than the kcat vector
itself:

- `corrections.tsv`: one row per changed kcat, ranked by leverage, with
  its prior value, tuned value, fold change and source.
- `kcat_leverage_screen.tsv`: what the data could speak to at all.
- `cmaes_trace.tsv`: best-so-far RMSE and objective per generation.
- `tuning_summary.txt`: how far the fit moved, by source.
- `seed_reproducibility.txt`: a second seed's run against the first.
- `kcat_list_brenda.tsv`, `kcat_list_okp.tsv`, `kcat_list_merged.tsv`:
  the candidate kcats before and after merging.

**Which corrections to trust.** A second seed reaches a similar fit
(0.892) but moves some kcats differently, so a kcat the seeds disagree on
is undetermined by the data. Of the 31 corrections beyond 2-fold, 22
(71%) go the same direction at the second seed and 29 (94%) land within
5-fold of each other. Of the eight beyond 5-fold, three agree closely
across seeds: `r_0079` (phosphoribosylformylglycinamidine synthase, 0.05
→ 0.37 and 0.39), `r_0109` (acetyl-CoA carboxylase, 1.23 → 8.1 and 7.1)
and `r_1134_EXP_6` (1.23 → 6.4 and 7.0). `r_0961` (pyruvate
dehydrogenase, 37.9 → 761 and 529) agrees on direction and rough size.
`r_2491` moves in opposite directions and the adenylate kinase pair
(`r_0149`, `r_0150`) disagrees on size, so those are undetermined.
`output/seed_reproducibility.txt` has the full table.

A correction that both seeds make, to a database value, is a candidate
for curation at the source: put it in `data/customKcats.tsv` with its
reference and rebuild. A kcat the seeds disagree on carries the fit
without saying anything about that enzyme.
