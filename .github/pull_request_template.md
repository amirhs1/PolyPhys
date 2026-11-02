<!--
This is the pull request's full report. Keep every section, in this order,
then the trailer block; a section that does not apply says `None`.
-->

<!-- Title: as CONTRIBUTING.md, "Names", sets it. -->

## Summary

What changed and why. The reason comes from a person or an outside report, not
from an AI tool; an agent gives it only as supplied.

## Related issues

`Closes #n` for each issue this pull request completes, or `None`.

## Problem

What was wrong or missing, with evidence.

## What changed

Files as `path:line`, plus reasoning the diff does not show.

## Checks run

Each check actually run, as `command → result`, and any human review or
independent reference used to judge correctness. Then a `Not verified:` line
for anything not checked, and why.

- [ ] `python -m flake8 polyphys`
- [ ] `python -m mypy polyphys/analyze polyphys/manage polyphys/simulate`
- [ ] `python -m pytest polyphys README.md --cov=polyphys --cov-report=term-missing --cov-report=xml --doctest-modules --doctest-glob='README.md'`
- [ ] `python -m build` if packaging metadata, package data, or package layout changed
- [ ] `python -m sphinx -W -b html docs/source docs/_build` if `docs/source/` changed

## Scientific correctness

- [ ] No measurement/statistics behavior changed
- [ ] Units, numerical fixtures, and domain assumptions are preserved or explained
- [ ] Parser lineage and organizer vocabulary are unchanged, or matching parser tests/docs were updated

## Decisions and risks

Choices made, the alternatives rejected, and what could break. Is the change a
one-way door, hard to undo, or a two-way door, and what is its blast radius?

## Notes for review

What needs line-by-line review, and any uncertainty, limitation, or follow-up.

<!--
AI assistance: end the body with one trailer block, after a blank line. Write
one `Assisted-by: <tool>, <model id or not recorded> (<role>)` line for each
AI tool used, then `Checks-run: <check actually run> — <observed result>` for
each check it ran. Do not guess a model; write `not recorded`. With no AI
tool, leave the block out. Do not include prompts, secrets, or personal data.
Responsibility for the contribution remains with the contributor; see
`AI-POLICY.md`.
-->
