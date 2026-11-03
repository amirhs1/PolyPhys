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
# Skip main and the branch this session is on; the first field is the name.
merged=$(
  git branch --merged origin/main \
    --format='%(refname:short)%(if)%(worktreepath)%(then) (checked out in %(worktreepath))%(end)' \
    | awk -v cur="${current:-main}" '$1 != "main" && $1 != cur'
)

[ -n "$merged" ] || exit 0
echo "Local branches already merged into origin/main; the close-out skill deletes them (git branch -d):"
echo "$merged"
