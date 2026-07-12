"""Real, published benchmarks -- measured communities muODE must reproduce.

A simulation that is only checked against a toy community it was designed to
reproduce is unfalsifiable.  These modules wire muODE to *measurements*: defined
communities of sequenced strains, grown in a published defined medium, with
species abundances and metabolite concentrations quantified experimentally.

Currently:

* :mod:`muode.benchmarks.clark2021` -- 1,850 synthetic human gut communities
  (25 strains, richness 1-23) grown in the chemically defined medium DM38, with
  measured acetate, butyrate, lactate and succinate.
"""

from muode.benchmarks import clark2021

__all__ = ["clark2021"]
