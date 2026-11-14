# Contributing to PolyPhys

## Reporting a problem

Bug reports, feature requests, and questions are welcome via
[issues](https://github.com/amirhs1/PolyPhys/issues); each has a template.
Report a suspected vulnerability privately instead, as
[`SECURITY.md`](SECURITY.md) describes; do not open a public issue for it.

## Development setup

```bash
git clone https://github.com/amirhs1/PolyPhys.git
cd PolyPhys
python -m pip install -e ".[dev,viz]"
```

The `viz` extra installs matplotlib, which the doctests of `polyphys/viz/`
need. Add the `mda` extra to run the probe tests, which are skipped without
MDAnalysis, and the `docs` extra to build the documentation.

### Checklist for each clone

- [ ] `git config core.hooksPath .githooks` has been run, which turns on the
      `commit-msg` hook.
- [ ] `git config fetch.prune true` has been run.
- [ ] `git config commit.template .gitmessage` has been run, for commits
      written in an editor.

## Tests and documentation checks

Run these from the repository root before opening a pull request:

| Check                | Command                                                                                                                                               |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| Lint                 | `python -m flake8 polyphys`                                                                                                                           |
| Type-check           | `python -m mypy polyphys/analyze polyphys/manage polyphys/probe polyphys/simulate`                                                                    |
| All tests            | `python -m pytest polyphys README.md --cov=polyphys --cov-report=term-missing --cov-report=xml --doctest-modules --doctest-glob='README.md'`          |
| One test             | `python -m pytest polyphys/tests/manage/test_parser.py::<test name>`                                                                                  |
| Package build        | `python -m build`, after a packaging change                                                                                                           |
| Documentation build  | `python -m sphinx -W -b html docs/source docs/_build`, when `docs/source/` changed                                                                    |

CI runs lint, type-check, and the tests on Python 3.11 to 3.14, and builds
the package and the documentation, on every pull request.

## Submitting a pull request

- Branch from `main`, and name the branch and the pull request title as
  "Names" below sets them.
- Pull requests merge with a merge commit, so each commit lands on `main`
  unchanged; keep every commit coherent.
- Fill in every section of the pull request template; it is the pull
  request's full report.

### Changelog

A pull request adds one line per change under `## [Unreleased]` in
`CHANGELOG.md`, in the form the file's header names, unless its only area is
`agents`. A pull request that only edits the changelog, and a release pull
request, add none.

Write the entry in the commit that makes the change, citing the issue the
pull request closes. A pull request that closes no issue cites itself in a
second commit, pushed once it is open. The release adds the commit
references.

## AI assistance

AI tools are welcome. Read the [AI policy](AI-POLICY.md) before opening a pull
request. The rules most often missed:

- Understand, and be able to explain, everything you submit.
- Say in the pull request which AI tools you used and for what, with
  `not recorded` for a model you cannot name.
- Write issues, pull request descriptions, and replies in your own words, or
  let an AI agent write them. Then you answer for every word, and its text
  ends with its `Assisted-by:` trailer.
- For text you write yourself, the pull-request statement is enough; the
  commit trailers are the maintainers' record, and optional for you.

Pull requests that don't follow the policy may be closed without review.

### Checking an AI tool's commits

- [ ] A test commit made by each tool in use ends with one trailer block that
      includes `Assisted-by:` and no AI `Co-authored-by:`
      (`git log -1 --format=%B`).

### Repository checks for the maintainer

- [ ] `.claude/settings.json` is committed with the attribution block.
- [ ] `.gitmessage` is committed.
- [ ] `.gitignore` has the agent-files block:

  ```gitignore
  # Agent files: local only
  CLAUDE.local.md
  AGENTS.local.md
  .claude/settings.local.json
  .claude/worktrees/
  .claude/.cc-writes/
  ```

## Names

Use a type to describe the kind of change and an area to identify the affected
part of the repository. A commit's scope is its area, and the labels are
`type:<type>` and `area:<area>`.

```text
types: feat fix docs test refactor chore release question
areas: analyze manage simulate probe viz packaging guide ci agents deps
```

| Type       | Means                                                                                              |
| ---------- | -------------------------------------------------------------------------------------------------- |
| `feat`     | Adds a capability to the package or its tooling.                                                   |
| `fix`      | Corrects behavior that was wrong.                                                                 |
| `docs`     | Changes documentation only, in any area: the README, a `docs/` page, the changelog, a docstring.   |
| `test`     | Adds or corrects tests without changing package code.                                              |
| `refactor` | Restructures code without changing what it does.                                                   |
| `chore`    | Maintenance that fits no other type, such as a dependency bump or a CI runner pin.                 |
| `release`  | Prepares or merges a release: the version, the changelog, and the citation metadata.               |
| `question` | Asks for information instead of proposing a change; an issue type only, never a branch or commit.  |

| Area        | Covers                                                                                                                            |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `analyze`   | `polyphys/analyze/`, `polyphys/tests/analyze/`                                                                                    |
| `manage`    | `polyphys/manage/`, `polyphys/tests/manage/`                                                                                      |
| `simulate`  | `polyphys/simulate/`, `polyphys/tests/simulate/`                                                                                  |
| `probe`     | `polyphys/probe/`, `polyphys/tests/probe/`                                                                                        |
| `viz`       | `polyphys/viz/`, `polyphys/tests/viz/`                                                                                            |
| `packaging` | `pyproject.toml`, `polyphys/__version__.py`, `polyphys/tests/test_packaging.py`, `CITATION.cff`, `AUTHORS.rst`                    |
| `guide`     | `README.md`, `CONTRIBUTING.md`, `CHANGELOG.md`, `SECURITY.md`, `docs/`, `.github/ISSUE_TEMPLATE/`                                 |
| `ci`        | `.github/workflows/`, `.github/dependabot.yml`                                                                                    |
| `agents`    | `AGENTS.md`, `CLAUDE.md`, `AI-POLICY.md`, `.agents/`, `.claude/`, `.githooks/`, `.gitmessage`, `.github/pull_request_template.md` |
| `deps`      | Dependency versions: the dependency lists in `pyproject.toml` and the action versions in `.github/workflows/`                     |

A change that only moves a dependency version takes the area `deps`, not
`packaging` or `ci`.

| Name                          | Form                                    |
| ----------------------------- | --------------------------------------- |
| Branch                        | `<type>/<short-name>`                   |
| Commit and pull request title | `<type>(<area>): <imperative summary>`  |
| Issue title                   | `[<area>] <imperative summary>`         |
| Epic                          | A parent issue, titled as an issue      |
| Milestone                     | `v<MAJOR>.<MINOR>.<PATCH>`, the release |
| Tag and release title         | `v<MAJOR>.<MINOR>.<PATCH>`              |

Labels combine: one `type:` plus one or more `area:`. There are no
combination labels. Status, priority, and release are not labels.

Issue templates set a default type: `type:fix` for a bug report, `type:feat`
for a feature request, and `type:question` for a question.

### Where these names are used

- `.gitmessage`
- `.github/pull_request_template.md`
- `.github/ISSUE_TEMPLATE/*`
- `.github/dependabot.yml`
- `.agents/skills/create-branch/SKILL.md`,
  `.agents/skills/draft-release/SKILL.md`,
  `.agents/skills/open-issue/SKILL.md`,
  `.agents/skills/open-pull-request/SKILL.md`,
  `.agents/skills/write-commit/SKILL.md`
- `AGENTS.md`, "Git"

Each of these cites this section.

## Releases

This section is the maintainer's release procedure.

A milestone is released from `main` as soon as its last issue closes, and
before any pull request that closes an issue of a later milestone merges.
Each release is tagged on `main`; never tag an older commit to leave later
work out. The `draft-release` skill checks this, and stops when `main`
already holds later work.

From v1.0.0, publishing the GitHub release also starts
`.github/workflows/publish.yaml`, which uploads the release to TestPyPI and
then to PyPI; a release marked as a pre-release goes to TestPyPI only.
After the upload, the `draft-release` skill installs the released version
from PyPI into a fresh environment outside any checkout and checks that
`import polyphys` gives that version from `site-packages`.

Once PolyPhys has a conda-forge feedstock, conda-forge's bot opens a
version-update pull request on it for each release; the maintainer, as the
feedstock's maintainer, reviews and merges it. The `draft-release` skill
checks that the pull request was opened.

CI runs neither the `private_data` tests nor the `lammps` tests, so before
each release the maintainer runs the full test command on their own machine
with `POLYPHYS_DATA_DIR` set and LAMMPS available, and the release commit
records the result in a `Checks-run:` line.
