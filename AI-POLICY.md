# AI Policy

PolyPhys is developed with the help of AI coding assistants. This document
states how they are used, who is accountable for the result, and what they are
not permitted to do.

## 1. Scope

This policy governs the `amirhs1/poly-phys` repository and the releases cut from
it. It is a governance document written for people — users, reviewers, and
researchers deciding whether to trust or cite this software.

It is not the operating contract for the assistants themselves. That contract
lives in [`AGENTS.md`](AGENTS.md), with tool-specific additions in
[`CLAUDE.md`](CLAUDE.md). Those files describe *how* an agent works in this
repository; this file describes *what is permitted and who answers for it*. Read
this one first; the contract is subordinate to it.

This policy does not govern the writing of manuscripts or the conduct of the
research that PolyPhys supports. Those follow the policies of the relevant
journal and institution.

## 2. Human accountability

The maintainer is accountable for every line merged into `main`, regardless of
which tool drafted it. AI assistance changes how the work is produced; it does
not transfer responsibility for the result.

Concretely:

- Every change is reviewed by the maintainer before it is merged. Assistants
  open draft pull requests; they do not mark them ready, approve them, or merge
  them.
- Only the maintainer publishes a release, moves a tag, or changes repository
  protections.
- **AI systems are not authors.** No assistant appears in
  [`CITATION.cff`](CITATION.cff) or [`AUTHORS.rst`](AUTHORS.rst), and none is
  credited as a contributor to a release. A tool cannot take responsibility for
  the work, which is what authorship means.

## 3. How AI is used here

Assistants are used to draft and revise code, tests, docstrings, and
documentation; to refactor existing code; to review diffs; to triage and draft
issues; and to write commit messages and pull-request descriptions.

They are not used to produce research results. Specifically, no assistant is
used to generate simulation data, to decide a scientific question, to invent a
numerical value that a test then enshrines, or to supply a citation that is
carried into the repository unchecked.

The assistants currently in use are Claude Code and OpenAI Codex. This list will
change; the rules in this document do not depend on which tool is in use.

## 4. Scientific integrity

PolyPhys produces numbers that end up in published work, so the ordinary risks
of generated code are compounded by a specific one: a plausible-looking value,
formula, or reference that no one ever verified. The following rules exist to
close that gap.

- **Expected values are derived, not asserted.** Every numerical fixture and
  expected result in the test suite must come from an analytic derivation, a
  published reference, or a reproducible computation. A value is never accepted
  because a model produced it, and a failing test is never "fixed" by replacing
  the expectation with the observed output.
- **Doctest output comes from running the code.** Examples in docstrings and in
  the README are executed by the test suite. Their output is copied from a real
  run, never predicted.
- **Citations are verified.** Physical models, statistical methods, algorithms,
  and nontrivial formulas carry a citation to a paper, textbook, standard, or
  official library document that has been confirmed to exist *and* to support
  the claim being made. Language models fabricate plausible references; an
  unverified citation is treated as a defect, not a formatting detail.
- **Units, shapes, and assumptions are preserved.** Physical units, array-shape
  contracts, numerical meaning, and established scientific assumptions are not
  changed by a refactor. A changed expected value must be explained, with its
  scientific basis, in the pull request that changes it.
- **Verification is observed, not assumed.** No check, test, build, or benchmark
  is reported as passing unless its successful result was actually seen.

## 5. Provenance and attribution

Commits and pull requests follow [`docs/NAMING-CONVENTION.md`](docs/NAMING-CONVENTION.md).
AI-assisted commits carry the assistant's own configured attribution trailer, so
the history records where a change came from without duplicating that
information in the subject line. Agent and tool names do not appear as prefixes
in commit subjects or pull-request titles.

## 6. Security and data handling

- Credentials, tokens, private keys, and configuration secrets are never
  provided to an AI tool, in any form.
- Unpublished simulation data, private datasets, and draft manuscripts are not
  pasted into prompts or otherwise sent to a third-party service. Assistants
  work against the repository, not against research data.
- Assistants operate within the permissions, sandboxing, hooks, and branch
  protections configured for them. Those controls are not to be weakened or
  circumvented, and a denied action is not to be retried by another route.
- Suspected vulnerabilities follow [`SECURITY.md`](SECURITY.md) and are reported
  privately. Exploit details do not go into a public issue or pull request,
  whether written by a person or a tool.

## 7. Licensing

PolyPhys is MIT-licensed, and everything merged into it must be distributable
under that license. Generated output is treated like any other contribution: it
must not reproduce substantial portions of code whose license is unknown or
incompatible. Verbatim reproduction of a recognizable third-party implementation
is a reason to stop and check the source, not a shortcut.

## 8. Contributors

Contributions are welcome, and using an AI assistant to prepare one is fine.
Two conditions apply:

1. **Disclose material assistance.** If an assistant drafted a substantial part
   of the change, say so in the pull-request description. A one-line note is
   enough.
2. **You remain responsible for what you submit.** You should be able to explain
   what the change does and why it is correct. A pull request consisting of
   unreviewed model output is not a contribution, and will be closed.

Everything in section 4 applies to contributed changes as well.

## 9. Revision

This policy will be revised as the tools, the project, and the surrounding norms
change. Substantive changes go through a pull request like any other change, and
are recorded in [`CHANGELOG.md`](CHANGELOG.md).
