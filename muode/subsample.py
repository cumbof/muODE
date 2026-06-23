"""Abundance-aware community subsampling -- scaling to huge MAG sets.

A metagenome can yield thousands of MAGs, but most of the community biomass (and
most of the metabolic activity that matters for cross-feeding) sits in a handful
of abundant taxa, with a long tail of rare bins.  Simulating every model is both
expensive and noisy.  These helpers pick the subset worth simulating by
abundance, by one of three policies (which compose):

* ``top_n``        -- keep the ``n`` most abundant MAGs;
* ``min_abundance``-- keep MAGs at or above a relative-abundance floor;
* ``coverage``     -- keep the fewest top MAGs whose cumulative relative
  abundance reaches ``coverage`` (e.g. 0.95 of the community).

The dropped tail is reported so nothing is silently discarded.
"""

from __future__ import annotations

from typing import Dict, List, Mapping, Optional, Tuple


def select_by_abundance(
    abundances: Mapping[str, float],
    top_n: Optional[int] = None,
    min_abundance: Optional[float] = None,
    coverage: Optional[float] = None,
) -> List[str]:
    """Return the ids to keep, most-abundant first.

    Thresholds are evaluated on *relative* abundance (the input need not be
    normalised).  When no policy is given, every id is kept.
    """
    pos = {i: float(v) for i, v in abundances.items() if float(v) > 0}
    total = sum(pos.values())
    if total <= 0:
        return list(abundances)
    rel = {i: v / total for i, v in pos.items()}
    ordered = sorted(rel, key=lambda i: rel[i], reverse=True)

    if coverage is not None:
        kept, cum = [], 0.0
        for i in ordered:
            kept.append(i)
            cum += rel[i]
            if cum >= coverage:
                break
        ordered = kept
    if min_abundance is not None:
        ordered = [i for i in ordered if rel[i] >= min_abundance]
    if top_n is not None:
        ordered = ordered[:top_n]
    return ordered


def subsample_abundances(
    abundances: Mapping[str, float],
    top_n: Optional[int] = None,
    min_abundance: Optional[float] = None,
    coverage: Optional[float] = None,
) -> Tuple[Dict[str, float], List[str]]:
    """Split abundances into ``(kept, dropped_ids)`` by :func:`select_by_abundance`."""
    keep = set(select_by_abundance(abundances, top_n, min_abundance, coverage))
    kept = {i: float(v) for i, v in abundances.items() if i in keep}
    dropped = [i for i in abundances if i not in keep]
    return kept, dropped
