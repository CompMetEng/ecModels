"""Build and tune ecYeastGEM from yeast-GEM 9.1.1.

Stages match the GECKO Nature Protocols paper and geckopy's
``tutorials/full_ecModel/protocol.py``; the differences are that kcat
predictions come from OpenKineticsPredictor (EITLEM) rather than a
shipped DLKcat file, and that Stage 3 fits kcats to measured growth
rates and exchange fluxes with CMA-ES instead of sensitivity tuning.

Each step checkpoints to ``models/`` or ``output/`` and is skipped when
its output already exists, so the long ones can be run separately::

    python build.py stage1      # yeast-GEM 9.1.1 -> ecModel structure
    python build.py kcats       # EC codes, BRENDA, SMILES
    python build.py okp-submit  # submit the EITLEM prediction job
    python build.py okp-fetch   # collect it, merge, finish Stage 2
    python build.py screen      # kcat leverage screen
    python build.py tune        # CMA-ES, writes models/ecYeastGEM.yml
    python build.py all         # every step above, in order

``okp-submit`` needs an OpenKineticsPredictor API key, from
``OKP_API_KEY`` or ``data/okpApiKey.txt`` (both untracked).
"""
from __future__ import annotations

import argparse
import logging
import hashlib
import multiprocessing
import pickle
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from geckopy import (
    ModelAdapter,
    apply_complex_data,
    apply_custom_kcats,
    apply_kcat_constraints,
    apply_kcat_list,
    assign_standard_kcat,
    fill_eccodes_from_database,
    fill_kcats_from_isozymes,
    find_met_smiles,
    fuzzy_kcat_matching,
    load_brenda_data,
    load_conventional_gem,
    load_ec_model,
    load_phyl_dist,
    load_uniprot_tsv,
    make_ec_model,
    merge_kcats,
    save_ec_model,
    set_prot_pool_size,
)
from geckopy.databases import download_phyl_dist, load_kegg_tsv
from geckopy.gather_kcats.open_kinetics_predictor import (
    fetch_open_kinetics_predictor,
    submit_open_kinetics_predictor,
)
from geckopy.kcat_tuning import annotate_from_model, corrections_tsv
from geckopy.kcat_tuning.evotune import (
    cmaes_kcat_tuning,
    parsimony,
    screen_kcat_leverage,
    select_tunable_mask,
)

ROOT = Path(__file__).parent
DATA = ROOT / "data"
MODELS = ROOT / "models"
OUTPUT = ROOT / "output"

#: Reaction kcats are taken from the first tier that has one: an exact
#: experimental measurement OKP found, then a high-confidence fuzzy
#: BRENDA match, then a UniProt-curated value, then the EITLEM
#: prediction, and only then a weak fuzzy BRENDA match.
#:
#: The predictor's entry is its *output* label, not the method id in
#: `[okp] method`: OKP answers an EITLEM job with "Prediction from
#: EITLEM-Kinetics", which `parse_okp_output` tags "OKP-EITLEM-Kinetics".
#: A name that does not match is dropped from the merge, with a warning.
SOURCE_PRIORITY = [
    "database_exact", "database_top", "uniprot", "EITLEM-Kinetics",
    "database_bottom",
]

#: The yeast-GEM release this ecModel is built from, and the checksum of
#: its ``model/yeast-GEM.yml``.
BASE_TAG = "v9.1.1"
BASE_SHA256 = "2ddfaaa73ead41ab243127d1b34fa56db2c86052f6eabcd12bb4eea93947dc17"
BASE_URL = (
    "https://raw.githubusercontent.com/SysBioChalmers/yeast-GEM/"
    f"{BASE_TAG}/model/yeast-GEM.yml"
)

logger = logging.getLogger("build")

sys.path.insert(0, str(ROOT))
from anaerobic import make_anaerobic  # noqa: E402


def ensure_base_model() -> Path:
    """Put the ``BASE_TAG`` release of yeast-GEM at ``models/yeast-GEM.yml``.

    Downloaded from the release on GitHub, and checked against
    ``BASE_SHA256`` whether it was just downloaded or was already there,
    so the ecModel cannot be built from a different version of the model
    than the one it describes.
    """
    path = MODELS / "yeast-GEM.yml"
    if not path.is_file():
        logger.info("Downloading %s", BASE_URL)
        with urllib.request.urlopen(BASE_URL, timeout=60) as response:
            data = response.read()
        MODELS.mkdir(exist_ok=True)
        path.write_bytes(data)
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    if digest != BASE_SHA256:
        raise SystemExit(
            f"{path} is not yeast-GEM {BASE_TAG} (sha256 {digest}); "
            "delete it to fetch the release again."
        )
    return path


def get_adapter() -> ModelAdapter:
    return ModelAdapter.from_folder(ROOT)


def _uniprot(adapter: ModelAdapter):
    return load_uniprot_tsv(DATA / "uniprot.tsv")


# --------------------------------------------------------------- stage 1


def stage1(adapter: ModelAdapter) -> None:
    """Expand yeast-GEM into an ecModel structure (STEP 8-14)."""
    out = MODELS / "ecYeastGEM_stage1.yml"
    if out.is_file():
        logger.info("stage1: %s exists, skipping.", out.name)
        return

    ensure_base_model()
    model = load_conventional_gem(adapter)
    logger.info(
        "yeast-GEM: %d reactions, %d metabolites, %d genes",
        len(model.reactions), len(model.metabolites), len(model.genes),
    )

    ec_model = make_ec_model(model, adapter, uniprot_db=_uniprot(adapter))
    logger.info(
        "ecModel: %d reactions, %d enzymes",
        len(ec_model.reactions), len(ec_model.ec.enzymes),
    )

    apply_complex_data(ec_model, path=DATA / "ComplexPortal.json")
    save_ec_model(ec_model, out.name, adapter=adapter)


# ---------------------------------------------------------------- kcats


def kcats(adapter: ModelAdapter) -> None:
    """EC numbers, fuzzy BRENDA matching and metabolite SMILES (STEP 16-22).

    Writes the BRENDA kcat list and a SMILES-annotated model, both of
    which the OKP steps below build on.
    """
    brenda_out = OUTPUT / "kcat_list_brenda.pkl"
    model_out = MODELS / "ecYeastGEM_stage2_smiles.yml"
    if brenda_out.is_file() and model_out.is_file():
        logger.info("kcats: outputs exist, skipping.")
        return

    ec_model = load_ec_model("ecYeastGEM_stage1.yml", adapter=adapter)

    uniprot_db = _uniprot(adapter)
    kegg_db = load_kegg_tsv(DATA / "kegg.tsv")
    fill_eccodes_from_database(ec_model, uniprot_db, kegg_db=kegg_db)

    phyl_dist_path = adapter.get_phyl_dist_path()
    if not phyl_dist_path.is_file():
        download_phyl_dist(phyl_dist_path)
    brenda = load_brenda_data(adapter.get_brenda_db_folder())
    kcat_list_brenda = fuzzy_kcat_matching(
        ec_model, brenda, load_phyl_dist(phyl_dist_path),
    )
    logger.info("Fuzzy BRENDA matches: %d rows", len(kcat_list_brenda))
    # Pickle is the checkpoint the next step reads: the substrate and gene
    # columns hold lists, which do not survive a TSV round trip. The TSV
    # beside it is for reading, not for reloading.
    kcat_list_brenda.to_pickle(brenda_out)
    kcat_list_brenda.to_csv(
        OUTPUT / "kcat_list_brenda.tsv", sep="\t", index=False,
    )

    find_met_smiles(ec_model, cache_path=DATA / "smilesDB.tsv")
    save_ec_model(ec_model, model_out.name, adapter=adapter)


# ------------------------------------------------------------------ OKP


def okp_submit(adapter: ModelAdapter) -> None:
    """Submit the EITLEM kcat-prediction job (STEP 23-24).

    Skipped when a job is already recorded or the predictions are already
    stored in ``data/OKP_output.csv``, which is all ``okp-fetch`` needs.
    """
    job_file = DATA / "OKP_job.txt"
    stored = DATA / "OKP_output.csv"
    for done in (job_file, stored):
        if done.is_file():
            logger.info("okp-submit: %s exists, skipping.", done.name)
            return
    ec_model = load_ec_model("ecYeastGEM_stage2_smiles.yml", adapter=adapter)
    job_id = submit_open_kinetics_predictor(ec_model, adapter=adapter)
    logger.info("Submitted OKP job %s", job_id)


def okp_fetch(adapter: ModelAdapter) -> None:
    """Collect the OKP result, merge every kcat source, finish Stage 2.

    Covers STEP 25-32: merge the prediction with the BRENDA matches,
    apply the curated kcats, propagate across isozymes, assign the
    standard kcat, and write the constraints into the stoichiometry.
    """
    out = MODELS / "ecYeastGEM_stage2.yml"
    if out.is_file():
        logger.info("okp-fetch: %s exists, skipping.", out.name)
        return

    ec_model = load_ec_model("ecYeastGEM_stage2_smiles.yml", adapter=adapter)
    done, kcat_list_okp = fetch_open_kinetics_predictor(
        ec_model,
        use_stored=(DATA / "OKP_output.csv").is_file(),
        adapter=adapter,
        wait=True,
        timeout=6 * 3600,
    )
    if not done or kcat_list_okp is None:
        raise SystemExit("OKP job has not finished; rerun okp-fetch later.")
    logger.info("OKP (EITLEM) predictions: %d rows", len(kcat_list_okp))
    kcat_list_okp.to_csv(OUTPUT / "kcat_list_okp.tsv", sep="\t", index=False)

    kcat_list_brenda = pd.read_pickle(OUTPUT / "kcat_list_brenda.pkl")
    merged = merge_kcats(
        kcat_list_brenda, kcat_list_okp, source_priority=SOURCE_PRIORITY,
    )
    merged.to_csv(OUTPUT / "kcat_list_merged.tsv", sep="\t", index=False)

    apply_kcat_list(ec_model, merged)
    apply_custom_kcats(ec_model, path=DATA / "customKcats.tsv")
    fill_kcats_from_isozymes(ec_model)
    assign_standard_kcat(ec_model, _uniprot(adapter))
    apply_kcat_constraints(ec_model)

    params = adapter.params
    set_prot_pool_size(
        ec_model, p_tot=params.p_tot, f=params.f, sigma=params.sigma,
    )
    save_ec_model(ec_model, out.name, adapter=adapter)

    ec_model.reactions.get_by_id(params.c_source).lower_bound = -1000
    ec_model.objective = params.bio_rxn
    growth = ec_model.optimize().fluxes[params.bio_rxn]
    logger.info(
        "Untuned growth rate on unlimited glucose: %.4f /hour (reference %.2f)",
        growth, params.gr_exp,
    )


# ------------------------------------------------------- Stage 3: tuning


def screen(adapter: ModelAdapter, n_proc: int) -> None:
    """Rank kcats by how much the measured data can speak to them.

    Skipped when the screen is already on disk as
    ``output/kcat_leverage_screen.pkl``, the copy ``tune`` reads. The
    tracked TSV beside it is the report and does not count.
    """
    out = OUTPUT / "kcat_leverage_screen.pkl"
    if out.is_file():
        logger.info("screen: %s exists, skipping.", out.name)
        return
    ec_model = load_ec_model("ecYeastGEM_stage2.yml", adapter=adapter)
    table = screen_kcat_leverage(
        ec_model, adapter=adapter, make_anaerobic=make_anaerobic,
        n_proc=n_proc,
    )
    table.drop(columns="_positions").to_csv(
        OUTPUT / "kcat_leverage_screen.tsv", sep="\t", index=False,
    )
    table.to_pickle(out)
    logger.info("Screened %d kcat groups", len(table))


def tune(
    adapter: ModelAdapter, n_proc: int, seed: int, target_impact_share: float,
) -> None:
    """Fit the high-leverage kcats to the measured data with CMA-ES.

    ``target_impact_share`` is the fraction of the total achievable
    improvement the searched kcats must carry (see ``select_tunable_mask``).
    """
    out = MODELS / "ecYeastGEM.yml"
    result_file = OUTPUT / "evotune_result.pkl"
    if result_file.is_file():
        logger.info("tune: %s exists, skipping.", result_file.name)
        return

    ec_model = load_ec_model("ecYeastGEM_stage2.yml", adapter=adapter)
    table = pd.read_pickle(OUTPUT / "kcat_leverage_screen.pkl")
    mask = select_tunable_mask(
        ec_model, table, target_impact_share=target_impact_share,
    )
    logger.info("Tuning %d of %d kcats", int(mask.sum()), len(mask))

    result = cmaes_kcat_tuning(
        ec_model, adapter=adapter, screen=table, tunable_mask=mask,
        make_anaerobic=make_anaerobic, n_proc=n_proc, seed=seed,
    )
    save_ec_model(ec_model, out.name, adapter=adapter)
    # Before anything that can fail: hours of search live in `result`,
    # and nothing else on disk holds the kcat vector it settled on.
    with open(OUTPUT / "evotune_result.pkl", "wb") as fh:
        pickle.dump(result, fh)
    _report(result, ec_model, table, seed)


def _report(result, ec_model, table, seed: int) -> None:
    """Write the artifacts a reviewer reads instead of the kcat vector."""
    leverage_by_rxn = dict(zip(table["rxn_id"], table["leverage"]))
    leverage = np.array([leverage_by_rxn.get(r, 0.0) for r in result.rxns])

    # annotate_from_model also returns `sources`, which corrections()
    # takes from result.groups itself and rejects as a keyword.
    annotations = annotate_from_model(ec_model, result.rxns)
    rows = result.corrections(
        leverage=leverage,
        names=annotations["names"], ec_codes=annotations["ec_codes"],
    )
    (OUTPUT / "corrections.tsv").write_text(corrections_tsv(rows))

    pd.DataFrame({
        "generation": np.arange(1, len(result.rmse_trace) + 1),
        "rmse": result.rmse_trace,
        "objective": result.objective_trace,
    }).to_csv(OUTPUT / "cmaes_trace.tsv", sep="\t", index=False)

    n_changed = parsimony.n_changed(result.new_kcat, result.old_kcat)
    summary = [
        f"seed: {seed}",
        f"generations: {result.n_generations}",
        f"converged: {result.converged}",
        f"rmse: {result.rmse_trace[0]:.4f} -> {min(result.rmse_trace):.4f}",
        f"kcats searched: {len(result.rxns)}",
        f"kcats changed: {n_changed}",
        f"corrections listed: {len(rows)}",
    ]
    for group, moves in parsimony.source_movement(
        result.new_kcat, result.old_kcat, result.groups,
        order=sorted(set(result.groups)),
    ).items():
        summary.append(f"movement[{group}]: {moves}")
    (OUTPUT / "tuning_summary.txt").write_text("\n".join(summary) + "\n")
    logger.info("\n%s", "\n".join(summary))


# ------------------------------------------------------------------ main

STEPS = {
    "stage1": lambda a, args: stage1(a),
    "kcats": lambda a, args: kcats(a),
    "okp-submit": lambda a, args: okp_submit(a),
    "okp-fetch": lambda a, args: okp_fetch(a),
    "screen": lambda a, args: screen(a, args.n_proc),
    "tune": lambda a, args: tune(
        a, args.n_proc, args.seed, args.target_impact_share,
    ),
}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("step", choices=[*STEPS, "all"])
    parser.add_argument("--n-proc", type=int, default=8)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--target-impact-share", type=float, default=0.95,
        help="fraction of the achievable improvement the tuned kcats carry",
    )
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    # Python 3.14 defaults to the forkserver start method, whose helper
    # process dies under Slurm and takes cobra's ProcessPool with it
    # (ConnectionResetError mid-run). The screen and the tuner are the
    # only parallel steps and both are plain forks of a read-only model.
    if "fork" in multiprocessing.get_all_start_methods():
        multiprocessing.set_start_method("fork", force=True)
    OUTPUT.mkdir(exist_ok=True)

    adapter = get_adapter()
    for name in (STEPS if args.step == "all" else [args.step]):
        logger.info("=== %s ===", name)
        STEPS[name](adapter, args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
