"""Perturbation engine -- antibiotics, knockouts and probiotic interventions.

A perturbation modifies the *baseline* bounds of one or more organisms before a
simulation starts (and the change persists for the whole run).  Three primitives
cover the README's use cases:

* **Reaction knockout** -- clamp specific reactions to zero flux.
* **Pathway attenuation** -- reduce the flux capacity of every reaction in a
  subsystem/pathway by a given *efficacy* (e.g. a broad-spectrum antibiotic that
  inhibits folate biosynthesis at 95% efficacy).  ``Vmax`` is multiplied by
  ``1 - efficacy``; an efficacy of 1.0 is a full knockout.
* **Species removal** -- wipe an entire organism (e.g. a targeted bacteriocin),
  letting the simulation reveal *secondary* extinctions when the removed species
  was feeding others via cross-feeding.

The interesting biology -- cascading starvation through a disrupted
cross-feeding network -- is then an emergent property of the dynamic-FBA run,
not something the perturbation hard-codes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, List, Sequence, Union

if TYPE_CHECKING:  # avoid an import cycle at runtime
    from muode.community import Community


def _selected(organism_id: str, selector: Union[str, Sequence[str]]) -> bool:
    if selector == "*" or selector is None:
        return True
    if isinstance(selector, str):
        return organism_id == selector
    return organism_id in selector


@dataclass
class PerturbationTarget:
    """A single perturbation action."""

    #: explicit reaction ids to hit (backend-agnostic)
    reactions: List[str] = field(default_factory=list)
    #: subsystem / pathway name to hit (resolved on cobra GEMs only)
    subsystem: str | None = None
    #: which organisms this applies to: "*", a single id, or a list of ids
    organisms: Union[str, Sequence[str]] = "*"
    #: fraction of capacity removed; 1.0 == full knockout
    efficacy: float = 1.0
    #: if True, remove the whole organism (all reactions to zero)
    remove_organism: bool = False

    def _factor(self) -> float:
        return max(0.0, 1.0 - float(self.efficacy))

    def apply_to(self, organism) -> None:
        if not _selected(organism.id, self.organisms):
            return

        if self.remove_organism:
            # clamp every reaction we can see to zero
            rxn_ids = list(getattr(organism, "reaction_ids", lambda: [])())
            if not rxn_ids and hasattr(organism, "model"):
                rxn_ids = [r.id for r in organism.model.reactions]
            organism.knock_out(rxn_ids)
            return

        targets = list(self.reactions)
        if self.subsystem and hasattr(organism, "reactions_in_subsystem"):
            targets += list(organism.reactions_in_subsystem(self.subsystem))

        if not targets:
            return
        if self.efficacy >= 1.0:
            organism.knock_out(targets)
        else:
            organism.scale_bounds(targets, self._factor())


@dataclass
class Perturbation:
    """A named collection of perturbation targets applied to a community."""

    name: str = "perturbation"
    targets: List[PerturbationTarget] = field(default_factory=list)

    def apply(self, community: "Community") -> None:
        for organism in community.organisms:
            for target in self.targets:
                target.apply_to(organism)

    def describe(self) -> str:
        parts = []
        for t in self.targets:
            if t.remove_organism:
                parts.append(f"remove {t.organisms}")
            elif t.subsystem:
                parts.append(f"{t.subsystem}@{t.efficacy:.0%} on {t.organisms}")
            else:
                parts.append(f"{len(t.reactions)} rxns@{t.efficacy:.0%} on {t.organisms}")
        return f"{self.name}: " + "; ".join(parts)

    # -- convenience constructors mirroring the CLI -------------------------
    @classmethod
    def antibiotic(
        cls,
        target_pathway: str,
        efficacy: float = 0.95,
        organisms: Union[str, Sequence[str]] = "*",
        name: str | None = None,
    ) -> "Perturbation":
        """Attenuate a pathway across (some of) the community."""
        return cls(
            name=name or f"antibiotic:{target_pathway}",
            targets=[
                PerturbationTarget(
                    subsystem=target_pathway, efficacy=efficacy, organisms=organisms
                )
            ],
        )

    @classmethod
    def knockout(
        cls,
        reactions: Sequence[str],
        organisms: Union[str, Sequence[str]] = "*",
        name: str | None = None,
    ) -> "Perturbation":
        """Fully knock out an explicit set of reactions."""
        return cls(
            name=name or "knockout",
            targets=[
                PerturbationTarget(
                    reactions=list(reactions), efficacy=1.0, organisms=organisms
                )
            ],
        )

    @classmethod
    def remove_species(
        cls, organisms: Sequence[str], name: str | None = None
    ) -> "Perturbation":
        """Remove one or more species entirely from the community."""
        return cls(
            name=name or "remove_species",
            targets=[PerturbationTarget(organisms=list(organisms), remove_organism=True)],
        )
