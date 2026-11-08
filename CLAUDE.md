@AGENTS.md

## Local and scoped instructions

- Claude Code does not read `AGENTS.local.md` on its own. On a machine that
  has one, a gitignored `CLAUDE.local.md` containing `@AGENTS.local.md` loads
  it at launch.
- Put genuinely path-specific Claude guidance in `.claude/rules/` so it loads
  only for matching files.
- `.claude/skills` is a symlink to `.agents/skills/`, the canonical home of
  this repository's skills, so Claude Code loads them on demand.

## Closing out a merged pull request

- The Claude desktop app's pull request monitor wakes a session on CI
  failures, merge conflicts, and review comments, not on a merge; on a merge
  it can only archive the session, when its auto-archive switch is on. So
  nothing runs the `close-out` skill by itself: run it when the maintainer
  says a pull request merged, or at the next session start when the hook
  below names a merged branch.
- The `SessionStart` hook in `.claude/settings.json` runs
  `.claude/hooks/list-merged-branches.sh`, which lists the local branches
  already merged into `origin/main`. It lists; the skill deletes. Branches it
  lists as having no commits of their own, such as a worktree branch before
  its first commit, are not merged pull requests and are not for close-out.
- In the Claude desktop app, archiving a session is reversible and is the
  app's action; whether it also removes the session's worktree depends on
  the worktree that `git worktree list` shows (the `close-out` skill, step
  2). Never delete a session or its conversation. End the skill's report
  with whether the session is safe to archive, and leave the archiving to
  the maintainer.

## Git attribution

Claude Code's own commit and pull-request attribution is turned off in the
committed `.claude/settings.json`. Do not rely on the deprecated
`includeCoAuthoredBy` setting.
