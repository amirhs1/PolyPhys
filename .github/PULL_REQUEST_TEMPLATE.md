<!--
Title: type(scope): imperative summary -- see docs/NAMING-CONVENTION.md
Put `Closes #N` in the Summary below, never in the title.
-->

## Summary

-

## Checks

- [ ] `flake8 polyphys`
- [ ] `mypy polyphys/analyze polyphys/manage`
- [ ] `pytest polyphys README.md --cov=polyphys --cov-report=term-missing --doctest-modules --doctest-glob="README.md"`
- [ ] `python -m build` if packaging metadata, package data, or package layout changed
- [ ] `python -m sphinx -W -b html docs/source docs/_build` if `docs/source/` changed

## Scientific Correctness

- [ ] No measurement/statistics behavior changed
- [ ] Units, numerical fixtures, and domain assumptions are preserved or explained
- [ ] Parser lineage and organizer vocabulary are unchanged, or matching parser tests/docs were updated

## Notes

-
