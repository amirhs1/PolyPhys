<!--
This is the pull request's full report. Keep every section, in this order,
then the trailer block; a section that does not apply says `None`.
-->

<!-- Title: as CONTRIBUTING.md, "Names", sets it. -->

## Summary

The outcome, in one or two sentences: what changed and why. The reason comes
from a person or an outside report, not from an AI tool; an agent gives it
only as supplied.

## Related issues

`Closes #n` for each issue this pull request completes, `Refs #n` for one it
covers in part, or `None`.

## Problem

What was wrong or missing, with evidence: an issue, a failing check, or an
observed output.

## What changed

How the change was made: files as `path:line`, plus reasoning the diff does
not show. Leave the outcome to the summary.

## Checks run

Each check actually run, as `command → result`, and any human review or
independent reference used to judge correctness. Then a `Not verified:` line
for anything not checked, and why.

- [ ] `python -m flake8 polyphys`
- [ ] `python -m mypy polyphys/analyze polyphys/manage polyphys/probe polyphys/simulate`
- [ ] `python -m pytest polyphys README.md --cov=polyphys --cov-report=term-missing --cov-report=xml --doctest-modules --doctest-glob='README.md'`
- [ ] `python -m build` if packaging metadata, package data, or package layout changed
- [ ] `python -m sphinx -W -b html docs/source docs/_build` if `docs/source/` changed
- [ ] `CHANGELOG.md` has an entry under `[Unreleased]`, or this pull request needs none (`CONTRIBUTING.md`, "Changelog")

## Scientific correctness

- [ ] Measurement and statistics behavior is unchanged, or the change is explained above with how it was validated
- [ ] Units, numerical fixtures, and domain assumptions are preserved or explained
- [ ] Parser lineage and organizer vocabulary are unchanged, or matching parser tests/docs were updated

## Decisions and risks

Trade-offs: the choices made and the alternatives rejected. Then answer:

- How would the change be undone, and what would undoing it cost?
- Which users or callers does it affect, and how?
- Does it change any numerical result, and which?
- Does it change any stored artifact, file name, or file format?

## Notes for review

Where the reviewer should look: what needs line-by-line review, and any
uncertainty, limitation, or follow-up.

<!--
AI assistance: say in this pull request which AI tools you used and for what;
do not guess a model, write `not recorded`. If you wrote this body yourself,
that statement is enough. If an AI agent wrote it, end it with one trailer
block, after a blank line: one
`Assisted-by: <tool>, <model id or not recorded> (<role>)` line for each AI
tool used, then `Checks-run: <check actually run> — <observed result>` for
each check it ran. With no AI tool, leave the block out. Do not include
prompts, secrets, or personal data.
Responsibility for the contribution remains with the contributor; see
`AI-POLICY.md`.
-->
