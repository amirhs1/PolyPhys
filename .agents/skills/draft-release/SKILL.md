---
name: draft-release
description: Draft a GitHub release and its changelog section; never publish it. Use only when the person running you asks for a release.
disable-model-invocation: true
---

# Draft a release

A release is tagged on `main`, once its milestone is complete and before work
of a later milestone merges, as CONTRIBUTING.md, "Releases", sets it. Never
tag an older commit to leave later work out; stop and report instead.

1. Check the branch, the working tree, and `HEAD`: on `main`, clean, and
   level with `origin/main`.
2. Tag and release title: both `v<MAJOR>.<MINOR>.<PATCH>`, as
   CONTRIBUTING.md, "Names", sets them. The person running you chooses the
   version.
3. Check that `main` holds this release and nothing later: the milestone
   `v<version>` has no open issue
   (`gh issue list --milestone v<version> --state open`), and each pull
   request merged since the last tag
   (`git log --first-parent --format='%h %s' <last tag>..HEAD`) closes no
   issue of another open milestone
   (`gh pr view <n> --json closingIssuesReferences`,
   `gh issue view <n> --json milestone`). If either fails, stop and report
   the issues and pull requests to the person running you.
4. List the changes since the last tag:
   `git log --no-merges --format='%h %s' <last tag>..HEAD`.
5. Ask the person running you to run the full test command of AGENTS.md,
   "Commands", on their own machine, with `-rs` added, `POLYPHYS_DATA_DIR`
   set to their data directory, and `lmp_serial` on `PATH` or `LMP` naming
   LAMMPS, and to report the summary line. Do not run it yourself: a failing
   `private_data` test can print data values. Stop if a test fails or any
   test is skipped. The bookkeeping commit of step 6 carries the result as
   `Checks-run: full test command with POLYPHYS_DATA_DIR and LMP set — <summary line>`.
6. On a new branch, make one bookkeeping commit: the version in
   `polyphys/__version__.py`, `version` and `date-released` in
   `CITATION.cff`, and the changelog section of step 7. Open a pull request;
   the release commit is its merge commit.
7. CHANGELOG.md: rename `## [Unreleased]` to a version heading whose version
   links to the release,
   `## [<MAJOR>.<MINOR>.<PATCH>](<release URL>) - <YYYY-MM-DD>`, and add an
   empty `## [Unreleased]` above it. In the same edit, make the
   `[Unreleased]` reference at the end of the file compare from the new tag:
   `[Unreleased]: https://github.com/amirhs1/PolyPhys/compare/v<version>...HEAD`.
   The section has no trailer; its commit carries one. Check its entries
   against step 4 in Common Changelog form:
   `### Changed`, `### Added`, `### Removed`, `### Fixed`, in that order,
   each only when it has a change; one imperative line per change,
   `**Breaking:**` lines first, citing its commit, and its pull request when
   there is one. Add a missing change, and the commit references the entries
   lack.
8. Release notes: the same categories without the heading, then your trailer
   block (AGENTS.md, "Provenance"). Once CI passes on the release commit,
   create the draft there, and never publish it:
   `gh release create v<version> --draft --target <full SHA> --title v<version> --notes-file <file>`.
   `--target` fixes the commit that the tag lands on when the release is
   published, even if `main` moves first.
9. To revise the draft: `gh release edit v<version> --notes-file <file>`.
10. After the maintainer publishes, confirm the new version under the
   concept record in `CITATION.cff`, `doi`, through the Zenodo API, not the
   webhook's response:
   `curl -s 'https://zenodo.org/api/records?q=conceptrecid:15858407&all_versions=true&sort=mostrecent&size=1'`.
11. After the maintainer publishes, confirm the upload. A release marked as a
   pre-release goes to TestPyPI only, and any other release to TestPyPI and
   then PyPI (CONTRIBUTING.md, "Releases"); see which it is with
   `gh release view v<version> --json isPrerelease`. The
   `.github/workflows/publish.yaml` run for this release, not merely the
   latest run, passed:
   `gh run list --workflow publish.yaml --event release --commit <full SHA>`,
   with the SHA of step 8, shows `success` for it. The registry serves the
   version: for a pre-release,
   `curl -s -o /dev/null -w '%{http_code}' https://test.pypi.org/pypi/polyphys/<version>/json`,
   and for a release,
   `curl -s -o /dev/null -w '%{http_code}' https://pypi.org/pypi/polyphys/<version>/json`,
   prints `200`. Then install it from that registry into a fresh environment
   in a temporary directory outside any checkout, where `import polyphys`
   cannot load a source tree instead. For a release:
   `cd "$(mktemp -d)" && python3 -m venv .venv && .venv/bin/python -m pip install -q polyphys==<version> && .venv/bin/python -c "import polyphys; print(polyphys.__version__, polyphys.__file__)"`;
   for a pre-release, which takes its dependencies from PyPI, add
   `--index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/`
   to the `pip install`. Either prints `<version>` and a path under
   `site-packages`.
12. After the PyPI upload of a release, not a pre-release, check
   conda-forge. If
   `curl -s -o /dev/null -w '%{http_code}' https://api.anaconda.org/package/conda-forge/polyphys`
   prints `404`, there is no feedstock yet; say so and stop. Otherwise find
   the bot's version-update pull request for `<version>`:
   `gh pr list --repo conda-forge/polyphys-feedstock --search <version>`.
   Report it to the maintainer, who reviews and merges it; never merge or
   comment on it yourself.
