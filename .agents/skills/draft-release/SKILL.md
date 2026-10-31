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
5. On a new branch, make one bookkeeping commit: the version in
   `polyphys/__version__.py`, `version` and `date-released` in
   `CITATION.cff`, and the changelog section of step 6. Open a pull request;
   the release commit is its merge commit.
6. CHANGELOG.md: rename `## [Unreleased]` to a version heading whose version
   links to the release,
   `## [<MAJOR>.<MINOR>.<PATCH>](<release URL>) - <YYYY-MM-DD>`, and add an
   empty `## [Unreleased]` above it. The section has no trailer; its commit
   carries one. Check its entries against step 4 in Common Changelog form:
   `### Changed`, `### Added`, `### Removed`, `### Fixed`, in that order,
   each only when it has a change; one imperative line per change,
   `**Breaking:**` lines first, citing its commit, and its pull request when
   there is one. Add a missing change, and the commit references the entries
   lack.
7. Release notes: the same categories without the heading, then your trailer
   block (AGENTS.md, "Provenance"). Once CI passes on the release commit,
   create the draft there, and never publish it:
   `gh release create v<version> --draft --target <full SHA> --title v<version> --notes-file <file>`.
   `--target` fixes the commit that the tag lands on when the release is
   published, even if `main` moves first.
8. To revise the draft: `gh release edit v<version> --notes-file <file>`.
9. After the maintainer publishes, confirm the new version under the
   concept record in `CITATION.cff`, `doi`, through the Zenodo API, not the
   webhook's response:
   `curl -s 'https://zenodo.org/api/records?q=conceptrecid:15858407&all_versions=true&sort=mostrecent&size=1'`.
