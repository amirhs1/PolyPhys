---
name: draft-release
description: Draft a GitHub release and its changelog section; never publish it. Use only when the person running you asks for a release.
disable-model-invocation: true
---

# Draft a release

1. Check the branch, the working tree, and `HEAD`: on `main`, clean, and
   level with `origin/main`.
2. Tag and release title: both `v<MAJOR>.<MINOR>.<PATCH>`, as
   CONTRIBUTING.md, "Names", sets them. The person running you chooses the
   version.
3. List the changes since the last tag:
   `git log --no-merges --format='%h %s' <last tag>..HEAD`.
4. Write the changes in Common Changelog form: `### Changed`, `### Added`,
   `### Removed`, `### Fixed`, in that order, each only when it has a change;
   one imperative line per change, citing its commit, and its pull request
   when there is one; no Unreleased section.
5. CHANGELOG.md: add them at the top under a version heading whose version
   links to the release,
   `## [<MAJOR>.<MINOR>.<PATCH>](<release URL>) - <YYYY-MM-DD>`, with no
   trailer; its commit carries one. Do this on a new branch, and open a pull
   request.
6. Release notes: the same categories without the heading, then your trailer
   block (AGENTS.md, "Provenance"). Create the draft, and never publish it:
   `gh release create v<version> --draft --title v<version> --notes-file <file>`.
7. To revise the draft: `gh release edit v<version> --notes-file <file>`.
