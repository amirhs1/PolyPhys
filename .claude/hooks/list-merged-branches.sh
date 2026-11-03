#!/bin/sh
# Claude Code SessionStart hook: list the local branches already merged into
# origin/main, so the session knows what the close-out skill can delete.
# It lists and never deletes. Its stdout goes into the session's context, so
# it prints nothing when there is nothing to report.
set -u

cd "${CLAUDE_PROJECT_DIR:-.}" || exit 0
git rev-parse --git-dir >/dev/null 2>&1 || exit 0

# Refresh origin/main and drop stale remote-tracking refs; carry on offline.
git fetch --quiet --prune origin main >/dev/null 2>&1 || true
git show-ref --verify --quiet refs/remotes/origin/main || exit 0

current=$(git branch --show-current 2>/dev/null)
# Skip main and the branch this session is on; the first field is the tip and
# the second the name.
branches=$(
  git branch --merged origin/main \
    --format='%(objectname) %(refname:short)%(if)%(worktreepath)%(then) (checked out in %(worktreepath))%(end)' \
    | awk -v cur="${current:-main}" '$2 != "main" && $2 != cur'
)
[ -n "$branches" ] || exit 0

# git branch --merged also matches a branch with no commits of its own, such
# as a worktree branch before its first commit. Pull requests merge with a
# merge commit (AGENTS.md, "Git"), so a merged branch's tip is never on the
# first-parent line of origin/main, and a branch with no commits always is.
mainline=$(git rev-list --first-parent origin/main)

# Print the branches whose tip is (1) or is not (0) on that line, without the
# tip.
pick() {
  printf '%s\n' "$branches" | MAINLINE=$mainline awk -v on_line="$1" '
    BEGIN { n = split(ENVIRON["MAINLINE"], c, "\n"); for (i = 1; i <= n; i++) line[c[i]] = 1 }
    (($1 in line) ? 1 : 0) == on_line { sub(/^[^ ]+ /, ""); print }'
}
merged=$(pick 0)
empty=$(pick 1)

if [ -n "$merged" ]; then
  echo "Local branches merged into origin/main by a merge commit; the close-out skill deletes them (git branch -d):"
  echo "$merged"
fi
if [ -n "$empty" ]; then
  echo "Local branches with no commits of their own; not merged pull requests, so not for close-out:"
  echo "$empty"
fi
