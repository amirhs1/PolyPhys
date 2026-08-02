# PolyPhys Naming Conventions

This file defines the naming rules for issues, labels, branches, commits, pull
requests, milestones, and releases in the `poly-phys` repository.

The goal is to make the GitHub Project easy to scan and to keep the history
consistent. These conventions are shared with `amirhs1/CareerDossierTeX`, so a
reader who knows one repository already knows the other.

## 1. Core rule

Use different naming styles for different GitHub objects:

| Object | Convention | Example |
|---|---|---|
| Issue title | `[area] Verb object` | `[organizer] Raise test coverage to 85 percent` |
| Epic issue title | `[epic] Release version goal` | `[epic] Release v0.5.0 trustworthy ensemble statistics` |
| Branch name | `type/short-description` | `test/organizer-coverage` |
| Commit message | `type(scope): imperative summary` | `test(organizer): cover the ensemble reducers` |
| Pull request title | `type(scope): imperative summary` | `fix(packaging): drop the broken console script` |
| Label | `type:*` or `area:*` | `type:test`, `area:manage` |
| Milestone | `vX.Y.Z — Release Name` | `v0.5.0 — Trustworthy Ensemble Statistics` |
| Tag | `vX.Y.Z` | `v0.4.0` |
| GitHub Release title | `PolyPhys vX.Y.Z — Release Name` | `PolyPhys v0.4.0 — Release Readiness` |

Do not force one convention onto every object. Issues, branches, commits, pull
requests, and labels serve different purposes.

---

## 2. Issue title convention

Use:

```text
[area] Verb object
```

Examples:

```text
[organizer] Replace skipped doctests with executable examples
[organizer] Raise test coverage from 13 percent to at least 85 percent
[analyze] Add autocorrelation-aware error estimation for correlated frames
[parser] Move per-project constants onto parser classes
[docs] Document time and space complexity for the organizer reducers
[packaging] Remove the unused runtime dependencies
[ci] Fail the build on Sphinx warnings
[release] Tag and archive v0.4.0
```

### Epic issue titles

Use lowercase `[epic]` for visual consistency with the other bracket prefixes:

```text
[epic] Release v0.5.0 trustworthy ensemble statistics
```

---

## 3. Branch naming convention

Use:

```text
type/short-description
```

Allowed branch types:

```text
feat/      fix/      docs/      test/
ci/        refactor/ chore/     release/
```

Examples:

```text
fix/packaging-metadata
docs/release-readiness
docs/naming-convention
test/organizer-coverage
feat/block-averaging
refactor/parser-dump-keyword
release/v0.4.0
```

Rules:

- Use lowercase.
- Use hyphens, not spaces or underscores.
- Keep the name short but specific.
- Match the branch type to the main purpose of the work.
- Do not include issue numbers unless you find them useful later.

---

## 4. Commit message convention

Use a lightweight Conventional Commits style:

```text
type(scope): imperative summary
```

Examples:

```text
fix(packaging): remove the broken console script
feat(analyze): add block-averaged error estimation
test(organizer): cover the ensemble reducers
docs(sphinx): generate the API reference
refactor(parser): move the dump keyword onto the parser class
ci(build): fail the docs job on Sphinx warnings
release: prepare v0.4.0
```

### Commit types

```text
feat      New user-facing or maintainer-facing feature
fix       Bug fix
docs      Documentation-only change
test      Tests, doctests, examples, or regression coverage
ci        GitHub Actions or automation changes
refactor  Code restructuring without changing public behavior
chore     Maintenance that does not fit another type
release   Version, changelog, tag, or release preparation
```

### Scopes

Use the scope to identify the part of the repository affected:

```text
parser      organizer   utils       types
measurer    analyze     manage
packaging   sphinx      docs
ci          github      release     agents
```

Rules:

- Write the summary in the imperative mood: `add`, `define`, `fix`, `prepare`.
- Keep the first line short.
- Do not combine unrelated changes in one commit.
- Prefer one coherent change per commit.

### Commit body

The body explains *why* before *what*. State the problem the change addresses,
then the change, then its boundaries. Report only checks that actually ran, with
their real outcomes. Put trailers in one final block separated by a blank line;
see the attribution rules in `AGENTS.md`.

---

## 5. Pull request title convention

Use the same style as commit messages:

```text
type(scope): imperative summary
```

Rules:

- A pull-request title describes the whole branch, not every small commit.
- Use `Closes #N` in the pull-request body, never in the title.
- Do not close a large epic from an early implementation pull request. Close the
  focused sub-issue instead.
- Use draft pull requests for unfinished branches that need CI or notes.
- Do not add agent or tool prefixes to pull-request titles.

---

## 6. Label naming convention

Use labels as metadata, not as titles.

### Type labels

Apply exactly one primary type label when possible:

```text
type:feature   type:bug       type:docs     type:test
type:ci        type:refactor  type:release  type:deps
type:question
```

The branch prefix determines the type label; see `AGENTS.md`.

### Area labels

Use one or more area labels when useful:

```text
area:manage    area:analyze   area:packaging
area:documentation            area:ci        area:agents
```

### State and contributor labels

Use only when needed:

```text
blocked   technical-debt   breaking-change   help-wanted
```

Rules:

- Do not use labels for status. Use the Project `Status` field.
- Do not use labels for priority. Use the Project `Priority` field.
- Do not use labels for release numbers. Use GitHub milestones.
- Do not duplicate information already shown by GitHub fields.

---

## 7. Milestone, tag, and release naming

Milestones represent releases:

```text
vX.Y.Z — Release Name
```

Examples:

```text
v0.4.0 — Release Readiness
v0.5.0 — Trustworthy Ensemble Statistics
v0.6.0 — Project Decoupling
v0.7.0 — Observables
v1.0.0 — Stable API
```

Tags use semantic versions (`vX.Y.Z`). GitHub Release titles use
`PolyPhys vX.Y.Z — Release Name`.

Rules:

- Issues and pull requests can belong to milestones.
- Do not create labels like `v0.4.0`; the milestone already tracks this.
- Only the maintainer publishes a release or moves a tag.
