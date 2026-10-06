@AGENTS.md

## Local and scoped instructions

- Claude Code does not read `AGENTS.local.md` on its own. On a machine that
  has one, a gitignored `CLAUDE.local.md` containing `@AGENTS.local.md` loads
  it at launch.
- Put genuinely path-specific Claude guidance in `.claude/rules/` so it loads
  only for matching files.
- Put recurring Claude-only procedures in `.claude/skills/` so they load on
  demand rather than expanding this file.

## Git attribution

Claude Code's own commit and pull-request attribution is turned off in the
committed `.claude/settings.json`. Do not rely on the deprecated
`includeCoAuthoredBy` setting.
