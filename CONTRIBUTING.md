# Contributing to PolyPhys

Bug reports and questions are welcome via
[issues](https://github.com/amirhs1/PolyPhys/issues). Please see
[`SECURITY.md`](SECURITY.md) for reporting suspected vulnerabilities privately.

## AI-assisted contributions

AI tools are welcome. Read the [AI policy](AI-POLICY.md) before opening a pull
request. The rules most often missed:

- Understand, and be able to explain, everything you submit.
- Say in the pull request which AI tools you used and for what.
- Write issues, pull request descriptions, and replies in your own words. If
  an AI agent writes them for you, you answer for every word, and its text
  ends with its `Assisted-by:` trailer, with `not recorded` for an unknown
  model.

Pull requests that don't follow the policy may be closed without review.

### Setup

- [ ] `.claude/settings.json` is committed with the attribution block.
- [ ] `git config core.hooksPath .githooks` has been run in this clone.
- [ ] `.gitmessage` is committed.
- [ ] `git config commit.template .gitmessage` has been run, for commits
      written in an editor.
- [ ] A test commit made by each tool in use ends with one trailer block that
      includes `Assisted-by:` and no AI `Co-authored-by:`
      (`git log -1 --format=%B`).
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

Two axes name every change. _Type_, the kind of change, is one word list, used
unchanged as the branch prefix, the commit type, the pull request title's
type, and the label `type:<word>`. _Area_, the part of the repository, is one
list, used unchanged as the commit scope, the issue title's bracket, and the
label `area:<word>`. Scope and area are one list.

```text
types: feat fix docs test refactor chore release question
areas: analyze manage packaging guide ci agents deps
```

| Type       | Means                                                                                              |
| ---------- | -------------------------------------------------------------------------------------------------- |
| `feat`     | Adds a capability to the package or its tooling.                                                   |
| `fix`      | Corrects behaviour that was wrong.                                                                 |
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
