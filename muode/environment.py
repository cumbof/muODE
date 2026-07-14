"""What kind of place are we simulating?

muODE used to answer that question with two bare constants in ``qc.py``::

    MAX_PLAUSIBLE_GROWTH = 1.5              # justified in its docstring by "gut anaerobes"
    FORBIDDEN_RESCUES = frozenset({"o2_e"}) # i.e. "this environment is anaerobic"

Both are **statements about a human colon**, and neither said so anywhere a caller
could see.  Run muODE on a soil or marine or skin sample and they are silently wrong:
a legitimate fast aerobe gets flagged as a modelling artefact, and a medium that
genuinely *should* contain oxygen never gets told that oxygen is what it is missing.
A framework that takes arbitrary genomes must take arbitrary environments too, and an
environmental assumption belongs in a named object you can print, not in a constant
nobody reads.

Read this before adding a preset
--------------------------------
``max_plausible_growth`` is a **ceiling on the absurd**, not a measured maximum.  Its
job is to catch a *modelling* failure -- a GEM handed 100 carbon sources at once and
"growing" at 6.5/h -- and it does that job by sitting comfortably ABOVE anything the
environment could really do, so that tripping it means the medium is broken rather
than the organism unusual.  Set it too tight and it censors biology; that is a worse
failure than setting it loose.

There is one number here that is not a judgement call:

    Nothing has ever been observed to grow faster than *Vibrio natriegens*, whose
    doubling time was measured at **9.8 minutes** in BHI at 37 C (Eagon RG, *J
    Bacteriol* 83:736-737, 1962) and which remains the record.  ln(2) / (9.8/60 h)
    = 4.24 / h.

So :data:`ABSOLUTE_MAX_GROWTH` is a hard physical ceiling for *any* environment, and
:class:`Environment` refuses to be constructed above it.  Everything below that is an
argument you have to make, which is why ``why`` is a required field and an
:class:`Environment` will not construct without one.

**Only environments we can actually justify are shipped.**  There is no ``rumen`` or
``marine`` preset here, not because they are uninteresting but because we have not
done the reading, and shipping a preset with a number we made up would be exactly the
failure this module exists to prevent -- it would just move the invented constant
somewhere more official-looking.  Write your own (the class is three fields), state
the ``why``, and cite it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, FrozenSet

#: ln(2) / (9.8 min / 60) -- the fastest doubling ever measured, in any organism,
#: under optimal laboratory conditions.  Eagon RG, "Pseudomonas natriegens, a marine
#: bacterium with a generation time of less than 10 minutes."  J Bacteriol 83:736-737
#: (1962).  No environment may claim a plausible growth rate above this.
ABSOLUTE_MAX_GROWTH = 4.24


@dataclass(frozen=True)
class Environment:
    """The non-chemical assumptions a simulation makes about its setting.

    The *chemistry* of an environment is the :class:`~muode.diet.Diet` -- what is
    there, how much, and how fast a cell may take it.  This is everything else: how
    fast anything could possibly grow here, and what a diagnosis is allowed to propose
    when a model fails to grow.

    ``forbidden_rescues`` is the subtle one.  :func:`muode.qc.diagnose_no_growth`
    completes a medium by asking "what single nutrient would restore growth?", and
    oxygen restores growth in almost any GEM.  In an anaerobic environment, proposing
    it does not *complete* the environment -- it *overturns* it, and it would drown
    every real finding.  In an aerobic environment there is nothing to forbid: oxygen
    is supposed to be there, and if the diet forgot it then saying so is exactly the
    right diagnosis.  So the set is a property of the setting, not of the code.
    """

    name: str
    max_plausible_growth: float
    why: str
    oxygen: bool = False
    forbidden_rescues: FrozenSet[str] = field(default_factory=frozenset)

    def __post_init__(self) -> None:
        if not self.why.strip():
            raise ValueError(
                f"{self.name}: an Environment must say WHY its growth ceiling is what it "
                "is.  This field exists because two unexplained constants in qc.py "
                "quietly assumed a human colon for every sample muODE ever ran."
            )
        if self.max_plausible_growth <= 0:
            raise ValueError(f"{self.name}: max_plausible_growth must be positive")
        if self.max_plausible_growth > ABSOLUTE_MAX_GROWTH:
            raise ValueError(
                f"{self.name}: max_plausible_growth={self.max_plausible_growth} exceeds "
                f"ABSOLUTE_MAX_GROWTH={ABSOLUTE_MAX_GROWTH}/h, the fastest doubling ever "
                "measured in any organism (V. natriegens, 9.8 min, Eagon 1962).  A ceiling "
                "above that cannot catch anything, because nothing can exceed it."
            )
        # An anaerobic environment must refuse oxygen as a rescue; an aerobic one must
        # not, or a diet that genuinely forgot oxygen could never be told so.
        if not self.oxygen and "o2_e" not in self.forbidden_rescues:
            object.__setattr__(
                self, "forbidden_rescues", frozenset(self.forbidden_rescues | {"o2_e"})
            )


#: The default, and the one muODE was silently assuming all along.
#:
#: 1.5/h is a ceiling on the absurd, NOT a measurement, and the honest reason is that
#: the measurement does not exist: precise in vivo growth rates for gut commensals have
#: never been established (Joseph et al., Nat Commun 14:5457, 2023 -- "gut bacterial
#: doubling times likely range from minutes-to-hours, although precise in vivo estimates
#: are not available").  What IS known is a lower bound from first principles: the colon
#: is a flow-through vessel, so anything that persists must grow faster than it washes
#: out.  In vitro, gut anaerobes typically run 0.1-0.5/h.
#:
#: 1.5 therefore sits ~3x above the fast end of what is seen in culture and ~3x below
#: the physical ceiling.  A GEM reporting more than this on a defined medium is being
#: handed more nutrient than a cell can consume -- which is precisely how the 89-MAG run
#: reached 6.5/h, a 6-minute doubling, faster than any organism ever measured.
HUMAN_GUT = Environment(
    name="human_gut",
    max_plausible_growth=1.5,
    oxygen=False,
    why=(
        "Anaerobic. Gut anaerobes run 0.1-0.5/h in culture; precise in vivo rates have "
        "never been measured (Joseph 2023, Nat Commun 14:5457). 1.5/h is a ceiling on the "
        "absurd -- well above any gut commensal, well below the 4.24/h physical limit -- "
        "so tripping it means the medium is wrong, not the organism unusual."
    ),
)

#: For an anaerobic setting we have NOT characterised.  Deliberately useless as a
#: filter: it only catches the physically impossible.  That is the honest default when
#: you do not know your environment -- it will not censor real biology, and it will not
#: pretend to knowledge we do not have.
GENERIC_ANAEROBIC = Environment(
    name="generic_anaerobic",
    max_plausible_growth=ABSOLUTE_MAX_GROWTH,
    oxygen=False,
    why=(
        "An anaerobic environment we have not characterised. The ceiling is the physical "
        "limit (V. natriegens, 9.8 min doubling, Eagon 1962), so this catches ONLY the "
        "impossible. If you know your setting, say so with a real Environment and a "
        "citation; do not tighten this number by guessing."
    ),
)

#: As above, but oxic: oxygen is a legitimate rescue, because it is supposed to be there.
GENERIC_AEROBIC = Environment(
    name="generic_aerobic",
    max_plausible_growth=ABSOLUTE_MAX_GROWTH,
    oxygen=True,
    why=(
        "An aerobic environment we have not characterised. The ceiling is the physical "
        "limit (Eagon 1962). Oxygen is NOT forbidden as a rescue here: an oxic medium "
        "that omits o2_e is a medium with a hole in it, and the diagnosis should say so."
    ),
)

_PRESETS: Dict[str, Environment] = {
    e.name: e for e in (HUMAN_GUT, GENERIC_ANAEROBIC, GENERIC_AEROBIC)
}


def load_environment(name: str) -> Environment:
    """Look up a preset by name.

    Refuses to guess.  If your environment is not here, that is a prompt to write it
    down (with a citation), not a prompt to reach for the nearest-looking preset.
    """
    try:
        return _PRESETS[name]
    except KeyError:
        raise KeyError(
            f"unknown environment {name!r}. Known: {sorted(_PRESETS)}. muODE deliberately "
            "ships only environments it can justify -- define your own Environment with a "
            "`why` and a citation rather than borrowing one that nearly fits."
        ) from None


def known_environments() -> Dict[str, Environment]:
    return dict(_PRESETS)
