"""Clark et al. 2021 -- 1,850 measured synthetic human gut communities.

The reference benchmark for muODE: 25 sequenced human gut strains, assembled into
631 distinct communities (richness 1-23), grown for 48 h in the chemically defined
medium **DM38**, with **measured** endpoint acetate, butyrate, lactate and
succinate (mM) and **measured** species relative abundances (16S).

    Clark RL, Connors BM, Stevenson DM, Hromada SE, Hamilton JJ, Amador-Noguez D,
    Venturelli OS.  "Design of synthetic human gut microbiome assembly and
    butyrate production."  Nat Commun 12, 3254 (2021).
    https://doi.org/10.1038/s41467-021-22938-y   (CC-BY 4.0)

Why this dataset and not a stool cohort
---------------------------------------
Everything is *known*: which strains are present (not inferred from reads), what
they were fed (DM38, published composition, not a guessed diet), and how long they
grew.  A prediction error is therefore a model error, not an uncertainty in the
input.  A stool cohort confounds the two.  It also has a published precedent --
community-scale metabolic models were benchmarked on this same data in Nat
Microbiol 7, 1-12 (2024), doi:10.1038/s41564-024-01728-4 -- so muODE's numbers are
comparable to a method already in the literature rather than to nothing.

Three tiers of increasing difficulty, deliberately separable
------------------------------------------------------------
1. **Monoculture phenotype** (`monoculture_phenotypes`) -- does each reconstructed
   GEM secrete the fermentation products its organism actually secretes?  No
   community dynamics involved, so a failure here is a *reconstruction* failure.
   This is the cleanest signal in the dataset and the first thing to get right.
2. **Pairwise** -- does the engine predict the interaction between two strains?
3. **High-richness assembly** -- does it predict composition and butyrate for
   communities of up to 23 species?

Two traps this module handles for you
-------------------------------------
* **Lactate is a medium component** (28.3 mM sodium lactate in DM38), not only a
  product.  Scoring predicted against measured *endpoint* lactate would mostly be
  scoring the medium.  :func:`observations` therefore also exposes each metabolite
  as a change from its DM38 baseline (`net_*`), which is what a model predicts.
* **Contaminated wells are flagged in the source data** and are excluded by
  default (`Contamination? != "No"`), as are wells with too few 16S reads.
"""

from __future__ import annotations

import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from muode.diet import Diet

#: The full measurement table from the authors' analysis repository.
DATA_URL = (
    "https://raw.githubusercontent.com/RyanLincolnClark/"
    "DesignSyntheticGutMicrobiomeAssemblyFunction/master/"
    "commonfiles/2020_02_28_MasterDF.csv"
)

#: The four metabolites quantified by HPLC, as muODE (BiGG) exchange metabolites.
#: Used when the community is reconstructed with CarveMe (BiGG namespace).
METABOLITES: Dict[str, str] = {
    "Acetate": "ac_e",
    "Butyrate": "but_e",
    "Lactate": "lac__L_e",
    "Succinate": "succ_e",
}

#: The same four, as ModelSEED exchange metabolites, for a gapseq-reconstructed
#: community.  gapseq is the namespace that can actually secrete butyrate (CarveMe's
#: BiGG universe cannot, so a CarveMe run scores AC/CC/ER/RI as non-producers).  A
#: prediction and its measured baseline must be compared in ONE namespace, so the medium
#: (dm38 vs dm38_modelseed) and this map switch together -- see ``metabolite_map`` /
#: ``medium`` and ``examples/benchmarks/clark2021/derive_dm38_modelseed.py``.
METABOLITES_MODELSEED: Dict[str, str] = {
    "Acetate": "cpd00029_e0",
    "Butyrate": "cpd00211_e0",
    "Lactate": "cpd00159_e0",
    "Succinate": "cpd00036_e0",
}


def metabolite_map(namespace: str = "bigg") -> Dict[str, str]:
    """The measured-column -> exchange-id map for the reconstruction's namespace."""
    if namespace == "bigg":
        return METABOLITES
    if namespace == "modelseed":
        return METABOLITES_MODELSEED
    raise ValueError(f"namespace must be 'bigg' or 'modelseed', not {namespace!r}")


@dataclass(frozen=True)
class Strain:
    """One strain of the panel.

    ``strain`` is load-bearing, not decoration.  Metabolic phenotype is a *strain*
    property: NCBI's species reference for *F. prausnitzii* is M21/2, but Clark et
    al. measured A2-165, and reconstructing the wrong one would validate a model
    of an organism nobody grew.  `fetch_genomes.py` matches on this field and
    shouts when it cannot.

    ``reference`` is the genome label the authors used themselves (their
    ``metadata_2019_06_17.py``), kept verbatim: it is what the 16S reads were
    mapped against, so it is what the measured abundances actually refer to.
    """

    abbreviation: str
    species: str
    strain: str
    reference: str
    #: The name to search NCBI under, when the organism has since been
    #: reclassified and the name in the 2021 paper no longer resolves.  Empty
    #: means "same as ``species``".
    ncbi_taxon: str = ""

    @property
    def taxon(self) -> str:
        return self.ncbi_taxon or self.species


#: Every strain code that appears in the measurement table, from the authors'
#: `namedict` and Supplementary Data 1.  Note **26**, not 25: ``HB``
#: (*Holdemanella biformis*) was used in the pairwise and COMM2/COMM3 experiments
#: but is not part of the authors' 25-species design space (:data:`DESIGN_SPECIES`).
#: It is included because dropping it would silently discard every community that
#: contains it.
STRAINS: Dict[str, Strain] = {s.abbreviation: s for s in [
    Strain("PC", "Prevotella copri", "DSM 18205", "Prevotella_copri_DSM_18205"),
    Strain("PJ", "Parabacteroides johnsonii", "DSM 18315", "Parabacteroides_johnsonii_DSM_18315_NZ_ABYH01000014"),
    Strain("BV", "Phocaeicola vulgatus", "ATCC 8482", "Bacteroides_vulgatus_ATCC_8482_NC_009614"),
    Strain("BF", "Bacteroides fragilis", "NCTC 9343", "Bacteroides_fragilis_NCTC_9343"),
    Strain("BO", "Bacteroides ovatus", "ATCC 8483", "Bacteroides_ovatus_ATCC_8483"),
    Strain("BT", "Bacteroides thetaiotaomicron", "VPI-5482", "Bacteroides_thetaiotaomicron_VPI-5482_NC_004663"),
    Strain("BC", "Bacteroides caccae", "ATCC 43185", "Bacteroides_caccae_ATCC_43185"),
    Strain("BY", "Bacteroides cellulosilyticus", "DSM 14838", "Bacteroides_cellulosilyticus_DSM_14838"),
    Strain("BU", "Bacteroides uniformis", "ATCC 8492", "Bacteroides_uniformis_ATCC_8492"),
    Strain("DP", "Desulfovibrio piger", "ATCC 29098", "Desulfovibrio_piger_ATCC_29098"),
    Strain("BL", "Bifidobacterium longum subsp. infantis", "ATCC 15697", "Bifidobacterium_longum_subsp_infantis"),
    Strain("BA", "Bifidobacterium adolescentis", "ATCC 15703", "Bifidobacterium_adolescentis_ATCC_15703_NC_008618"),
    Strain("BP", "Bifidobacterium pseudocatenulatum", "DSM 20438", "Bifidobacterium_pseudocatenulatum_DSM20438"),
    Strain("CA", "Collinsella aerofaciens", "ATCC 25986", "Collinsella_aerofaciens_ATCC_25986"),
    Strain("EL", "Eggerthella lenta", "DSM 2243", "Eggerthella_lenta_DSM_2243_NC_013204"),
    # Reclassified since the paper: F. prausnitzii A2-165 is now *Faecalibacterium
    # duncaniae* A2-165.  Searching NCBI under the 2021 name cannot find it, and
    # silently returns the species reference (M21/2) -- a different organism.
    Strain("FP", "Faecalibacterium prausnitzii", "A2-165", "Faecalibacterium_prausnitzii_A2_165",
           ncbi_taxon="Faecalibacterium duncaniae"),
    Strain("CH", "Clostridium hiranonis", "DSM 13275", "Clostridium_hiranonis_DSM_13275"),
    Strain("AC", "Anaerostipes caccae", "DSM 14662", "Anaerostipes_caccae_DSM_14662"),
    Strain("BH", "Blautia hydrogenotrophica", "DSM 10507", "Blautia_hydrogenotrophica_DSM_10507"),
    Strain("CG", "Clostridium asparagiforme", "DSM 15981", "Clostridium_asparagiforme_DSM_15981"),
    Strain("ER", "Eubacterium rectale", "ATCC 33656", "Eubacterium_rectale_ATCC_33656_NC_012781"),
    Strain("RI", "Roseburia intestinalis", "L1-82", "Roseburia_intestinalis_L1_82"),
    Strain("CC", "Coprococcus comes", "ATCC 27758", "Coprococcus_comes_ATCC_27758"),
    Strain("DL", "Dorea longicatena", "DSM 13814", "Dorea_longicatena_DSM_13814"),
    Strain("DF", "Dorea formicigenerans", "ATCC 27755", "Dorea_formicigenerans_ATCC_27755"),
    Strain("HB", "Holdemanella biformis", "DSM 3989", "Holdemanella_biformis_DSM_3989"),
]}

#: The authors' 25-species design space (their ``allspecies``): the strains the
#: community-design experiments were drawn from.  Excludes ``HB``.
DESIGN_SPECIES: Tuple[str, ...] = (
    "ER", "FP", "AC", "CC", "RI", "EL", "CH", "DP", "BH", "CA", "PC", "PJ", "DL",
    "CG", "BF", "BO", "BT", "BU", "BV", "BC", "BY", "DF", "BL", "BP", "BA",
)

#: The authors' designated butyrate producers (their ``bpbspecies``) and succinate
#: producers (``spbspecies``).  Kept only as a *cross-check* on the phenotypes
#: :func:`monoculture_phenotypes` derives from the measurements -- the data, not
#: this tuple, is the ground truth a model is scored against.
AUTHORS_BUTYRATE_PRODUCERS: Tuple[str, ...] = ("ER", "FP", "AC", "CC", "RI")
AUTHORS_SUCCINATE_PRODUCERS: Tuple[str, ...] = ("PJ", "BT", "BF", "BC", "BO", "BV", "PC")


def dm38() -> Diet:
    """The DM38 defined medium the communities were actually grown in (BiGG namespace).

    Derived from Supplementary Data 4 -- see
    `examples/benchmarks/clark2021/derive_dm38.py`.
    """
    from muode.diet import _DIET_DIR

    return Diet.from_csv(_DIET_DIR / "dm38.csv", name="DM38")


def dm38_modelseed() -> Diet:
    """DM38 in the ModelSEED namespace, for feeding a gapseq-reconstructed community.

    Not committed: it is derived from the gapseq GEMs' own annotations (a workstation
    step), so it depends on which strains were reconstructed.  Build it with
    ``examples/benchmarks/clark2021/derive_dm38_modelseed.py --gems <dir>``.
    """
    from muode.diet import _DIET_DIR

    path = _DIET_DIR / "dm38_modelseed.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} does not exist.  It is a workstation artifact: reconstruct the "
            "clark strains with gapseq, then run "
            "examples/benchmarks/clark2021/derive_dm38_modelseed.py --gems <dir>."
        )
    return Diet.from_csv(path, name="DM38_modelseed")


def medium(namespace: str = "bigg") -> Diet:
    """The DM38 medium in the reconstruction's namespace (see :func:`metabolite_map`)."""
    if namespace == "bigg":
        return dm38()
    if namespace == "modelseed":
        return dm38_modelseed()
    raise ValueError(f"namespace must be 'bigg' or 'modelseed', not {namespace!r}")


def fetch(dest: str | Path) -> Path:
    """Download the measurement table (~1 MB) to ``dest``.

    Not vendored into the repository: it is the authors' data, and fetching it
    keeps the provenance a URL rather than a copy.
    """
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(DATA_URL) as fh:  # noqa: S310  (fixed, https)
        dest.write_bytes(fh.read())
    return dest


def load(path: str | Path, drop_contaminated: bool = True):
    """Load the measurement table as a DataFrame.

    Parameters
    ----------
    drop_contaminated:
        Exclude wells the authors flagged as contaminated or as having too few
        16S reads.  On by default -- these are known-bad measurements, and
        scoring a model against them measures the pipetting, not the model.
    """
    import pandas as pd

    df = pd.read_csv(path, index_col=0)
    if drop_contaminated and "Contamination?" in df.columns:
        df = df[df["Contamination?"] == "No"].copy()

    # `Treatment` is a '-'-joined set of strain codes: "AC", "DP-BL-FP", ...
    # The authors' named pools ("COMM6") and their leave-one-out dropouts
    # ("COMM6*AC") are excluded: the label alone does not say who is in them, so
    # we cannot build the community to simulate.  ~200 wells; they are the only
    # measurements this loader ignores, and it ignores them because they are not
    # actionable, not because they are inconvenient.
    members = df["Treatment"].astype(str).str.split("-")
    known = members.apply(lambda ms: all(m in STRAINS for m in ms))
    return df[known].copy()


@dataclass(frozen=True)
class Observation:
    """What was measured for one community (mean over replicate wells)."""

    community: str                      #: e.g. "AC-BT-RI"
    species: Tuple[str, ...]            #: strain codes present
    n_replicates: int
    metabolites: Dict[str, float]       #: BiGG id -> measured endpoint, mM
    net_metabolites: Dict[str, float]   #: BiGG id -> measured minus DM38 baseline
    abundances: Dict[str, float]        #: strain code -> measured relative abundance
    od: float

    @property
    def richness(self) -> int:
        return len(self.species)


def observations(df, diet: Optional[Diet] = None,
                 metabolites: Optional[Dict[str, str]] = None) -> List[Observation]:
    """Collapse replicate wells into one :class:`Observation` per community.

    ``net_metabolites`` subtracts the DM38 starting concentration, so a positive
    value is net production and a negative one net consumption.  This is the
    quantity a model predicts; the raw endpoint of lactate is dominated by the
    28.3 mM already in the medium.

    ``metabolites`` (and ``diet``) select the namespace: pass
    ``metabolite_map("modelseed")`` + ``dm38_modelseed()`` to key the measured values by
    the same ModelSEED exchange ids a gapseq community's prediction uses.  They MUST
    agree -- a prediction keyed cpd00211_e0 cannot be compared to a truth keyed but_e.
    """
    diet = diet or dm38()
    metabolites = metabolites or METABOLITES
    out: List[Observation] = []

    for community, grp in df.groupby("Treatment", sort=True):
        species = tuple(str(community).split("-"))
        mets, net = {}, {}
        for column, exch in metabolites.items():
            if column not in grp:
                continue
            measured = float(grp[column].mean())
            mets[exch] = measured
            net[exch] = measured - diet.initial_concentration(exch)

        abundances = {}
        for code in species:
            col = f"{code} Fraction"
            if col in grp:
                value = grp[col].mean()
                if value == value:                      # not NaN
                    abundances[code] = float(value)

        out.append(Observation(
            community=str(community),
            species=species,
            n_replicates=int(len(grp)),
            metabolites=mets,
            net_metabolites=net,
            abundances=abundances,
            od=float(grp["OD"].mean()) if "OD" in grp else float("nan"),
        ))
    return out


def monoculture_growth(df) -> Dict[str, float]:
    """Median monoculture OD per strain -- who actually grew in DM38.

    Needed because a strain that did not grow cannot secrete anything, so its
    measured phenotype says nothing about its metabolism.  The clearest case is
    *F. prausnitzii* (FP): a canonical butyrate producer that reaches OD ~0.02 in
    DM38 alone and therefore measures as a non-producer.  Scoring a model against
    that without excluding it punishes the model for being right.
    """
    mono = df[~df["Treatment"].astype(str).str.contains("-")]
    return {
        str(code): float(grp["OD"].median())
        for code, grp in mono.groupby("Treatment", sort=True)
        if str(code) in STRAINS
    }


def non_growers(df, min_od: float = 0.1) -> Tuple[str, ...]:
    """Strains that did not appreciably grow in DM38 monoculture (default OD < 0.1).

    Pass to :func:`score_phenotypes` as ``exclude`` -- and *say in the paper* that
    you did, with the list.  Quietly dropping them would be tuning the benchmark.
    """
    return tuple(sorted(c for c, od in monoculture_growth(df).items() if od < min_od))


def monoculture_phenotypes(df, threshold: float = 5.0, diet: Optional[Diet] = None,
                           metabolites: Optional[Dict[str, str]] = None
                           ) -> Dict[str, Dict[str, bool]]:
    """Ground-truth secretion phenotype per strain, **derived from the data**.

    For each strain grown in monoculture, whether it produced each metabolite
    above ``threshold`` mM *net* of the DM38 baseline.  Nothing here is asserted
    from memory or from the literature -- it is read off the measurements -- so
    the reference a model is scored against cannot be a transcription error.

    This is the sharpest test in the dataset.  It is also where a reconstruction
    fails loudly and interpretably: a GEM for *Bacteroides thetaiotaomicron* that
    secretes butyrate is wrong, and no amount of community dynamics will fix it.

    Caveat worth carrying into any figure: *F. prausnitzii* (FP) barely grew in
    DM38 (OD ~0.06), so its measured butyrate is ~0 even though it is a canonical
    butyrate producer in vivo.  That is a property of the *medium*, not of the
    organism -- and a model that grows FP well will "fail" this test for a
    defensible reason.  Report it; do not tune it away.
    """
    diet = diet or dm38()
    metabolites = metabolites or METABOLITES
    mono = df[~df["Treatment"].astype(str).str.contains("-")]

    phenotypes: Dict[str, Dict[str, bool]] = {}
    for code, grp in mono.groupby("Treatment", sort=True):
        code = str(code)
        if code not in STRAINS:
            continue
        phenotypes[code] = {
            exch: bool(float(grp[column].median()) - diet.initial_concentration(exch) > threshold)
            for column, exch in metabolites.items() if column in grp
        }
    return phenotypes


# ---------------------------------------------------------------------------
# Scoring predictions against the measurements
# ---------------------------------------------------------------------------


def score_metabolites(
    predicted: Dict[str, Dict[str, float]],
    observed: Sequence[Observation],
    metabolite: str = "but_e",
    net: bool = True,
) -> dict:
    """Score predicted vs measured concentrations of one metabolite, across communities.

    Parameters
    ----------
    predicted:
        ``{community: {bigg_id: mM}}`` -- muODE's endpoint (or net) prediction.
    net:
        Compare net production (measured minus DM38 baseline) rather than the raw
        endpoint.  Leave this on unless you know why you are turning it off: for
        lactate the raw endpoint is mostly medium.

    Returns Pearson r, Spearman rho, MAE and RMSE across communities -- the
    metrics this dataset is scored on in the literature, so the numbers are
    directly comparable.
    """
    import numpy as np

    pairs = []
    for obs in observed:
        if obs.community not in predicted:
            continue
        truth = (obs.net_metabolites if net else obs.metabolites).get(metabolite)
        guess = predicted[obs.community].get(metabolite)
        if truth is None or guess is None:
            continue
        pairs.append((float(guess), float(truth)))

    if len(pairs) < 2:
        return {"metabolite": metabolite, "n": len(pairs), "error": "too few paired communities"}

    pred = np.array([p for p, _ in pairs])
    obs_ = np.array([o for _, o in pairs])
    residual = pred - obs_

    result = {
        "metabolite": metabolite,
        "net": net,
        "n": len(pairs),
        "mae": float(np.mean(np.abs(residual))),
        "rmse": float(np.sqrt(np.mean(residual ** 2))),
        "bias": float(np.mean(residual)),          # +ve: the model over-predicts
        "observed_mean": float(np.mean(obs_)),
        "observed_std": float(np.std(obs_)),
    }
    try:
        import warnings

        from scipy.stats import pearsonr, spearmanr

        with warnings.catch_warnings():
            # a model that predicts the same value everywhere is degenerate, not
            # exceptional: correlation is undefined, so report None rather than warn
            warnings.simplefilter("ignore")
            r = float(pearsonr(pred, obs_).statistic)
            rho = float(spearmanr(pred, obs_).correlation)
        result["pearson_r"] = None if r != r else r
        result["spearman_rho"] = None if rho != rho else rho
    except Exception:                               # scipy is optional
        result["pearson_r"] = result["spearman_rho"] = None
    return result


def score_phenotypes(
    predicted: Dict[str, Dict[str, bool]],
    observed: Dict[str, Dict[str, bool]],
    metabolite: str = "but_e",
    exclude: Sequence[str] = (),
) -> dict:
    """Confusion matrix for a secretion phenotype across strains (tier 1).

    ``predicted``/``observed`` are ``{strain_code: {bigg_id: secretes?}}``.
    ``exclude`` drops strains from scoring -- use it for the DM38 non-growers
    (:func:`non_growers`), whose measured phenotype reflects the medium rather
    than their metabolism.  The excluded set is echoed back in the result so it
    always travels with the score.
    """
    codes = sorted((set(predicted) & set(observed)) - set(exclude))
    tp = fp = tn = fn = 0
    wrong: List[str] = []
    for code in codes:
        want = observed[code].get(metabolite)
        got = predicted[code].get(metabolite)
        if want is None or got is None:
            continue
        if got and want:
            tp += 1
        elif got and not want:
            fp += 1
            wrong.append(f"{code} (predicted to secrete, does not)")
        elif not got and want:
            fn += 1
            wrong.append(f"{code} (secretes, not predicted)")
        else:
            tn += 1

    total = tp + fp + tn + fn
    precision = tp / (tp + fp) if (tp + fp) else None
    recall = tp / (tp + fn) if (tp + fn) else None
    f1 = (2 * precision * recall / (precision + recall)) if (precision and recall) else 0.0
    return {
        "metabolite": metabolite,
        "n_strains": total,
        "excluded": sorted(exclude),
        "accuracy": (tp + tn) / total if total else None,
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "true_positive": tp, "false_positive": fp,
        "true_negative": tn, "false_negative": fn,
        "misclassified": wrong,
    }
