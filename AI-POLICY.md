# AI Policy - PolyPhys Project

Last reviewed: 2026-10-03

AI tools are welcome here. They don't change who is responsible: whoever
submits a change must understand it, have checked it, and be able to explain
it. This applies to maintainers too. How AI was used in this project, and how
each part was checked, is described in the [README](README.md).

No AI tool is an author or co-author of PolyPhys: none is listed in
[`CITATION.cff`](CITATION.cff) or [`AUTHORS.rst`](AUTHORS.rst), or named in a
`Co-authored-by:` trailer.

## Scope

This policy governs the `amirhs1/PolyPhys` repository and the releases cut from
it. This policy does not govern the writing of manuscripts or the conduct of
the research that PolyPhys supports. Those follow the policies of the relevant
journal and institution.

## Verification and scientific integrity

PolyPhys produces numbers that end up in published work, so the main risk of
generated code here is a plausible-looking value, formula, or reference that no
one verified. These rules close that gap.

- **Expected values are derived, not asserted.** Test values and numerical
  fixtures come from outside the AI's own output: an analytic derivation, the
  literature, measured data, or an independent implementation. If none exists,
  test a property (symmetry, conservation, invariance) and say so.
- **Tests are not weakened to pass.** Do not delete or loosen a test, or
  replace an expectation with the observed output, to make it pass.
- **Doctest output comes from running the code.** Examples in docstrings and in
  the README are executed by the test suite. Their output is copied from a real
  run, never predicted.
- **Citations are verified.** Physical models, statistical methods, algorithms,
  and nontrivial formulas carry a citation to a paper, textbook, standard, or
  official library document that has been confirmed to exist _and_ to support
  the claim being made. Language models fabricate plausible references; an
  unverified citation is treated as a defect, not a formatting detail. An
  AI-suggested reference is a lead, not evidence, until it is checked against
  its source.
- **Units, shapes, and assumptions are preserved.** Physical units, array-shape
  contracts, numerical meaning, and established scientific assumptions are not
  changed by a refactor. A changed expected value must be explained, with its
  scientific basis, in the pull request that changes it.
- **Verification is observed, not assumed.** No check, test, build, or
  benchmark is reported as passing unless its successful result was seen, and
  every reported number traces back to the code or source that produced it.
- **Domain choices are made by a person.** Observables, estimators, and fitting
  ranges are chosen by a person, not an AI tool.

## Disclosure

- Say in the pull request which AI tools you used and for what. A maintainer's
  pull request repeats its commits' `Assisted-by:` lines. If you don't know
  which model was used, write `not recorded`; don't guess.
- Maintainers record substantial AI help in commits with an `Assisted-by:`
  trailer, and every commit an AI agent writes carries one. Add `Checks-run:`
  only for a check actually run, with its observed result. Add
  `Ground-truth-source:` only when a commit adds or changes a reference value,
  naming its independent source. Outside contributors may use these trailers
  too, but their pull-request statement is enough.

  ```text
  Assisted-by: <tool>, <model identifier or not recorded> (<role or extent>)
  Checks-run: <check actually run> — <observed result>
  Ground-truth-source: <independent source of a reference value>
  ```

  Omit trailers that do not apply. A property test without a reference value
  does not need `Ground-truth-source:`.

## Communication

Write issues, pull request descriptions, and replies in your own words. AI may
fix grammar or translate. An agent may draft a pull request description at a
maintainer's request; the maintainer reviews it before the pull request is
marked ready. The reason a change exists — in a commit, pull request, or
changelog — comes from a person, not from the AI.

## Licensing and data

- You must have the right to submit what you submit. AI output that reproduces
  someone else's code is their code: attribute it under its licence or replace
  it.
- Do not give AI tools credentials, private or restricted data, unpublished
  simulation data, draft manuscripts, or other material you are not allowed to
  share.
- Report suspected vulnerabilities privately, as [`SECURITY.md`](SECURITY.md)
  describes, never in a public issue or pull request.

## Agents

AI agents act only with a person's approval: no autonomous issues, pull
requests, or comments. Instructions for agents working in this repository are
in `AGENTS.md`, with tool-specific additions in files like
[`CLAUDE.md`](CLAUDE.md). Agents work within the permissions, sandboxes, hooks,
and branch protections configured for them, and never work around a denied
action. Only a maintainer marks a pull request ready, merges it, publishes a
release, or changes repository settings.

## Enforcement

Maintainers may close a contribution that does not follow this policy without a
full review. Changes to this policy go through a pull request and are recorded
in [`CHANGELOG.md`](CHANGELOG.md).
