---
name: close-out
description: Clean up the clone after a pull request merges and its remote branch is deleted, and say whether the session that made it is safe to archive. Use when a person says a pull request merged, or when the session-start listing names a merged branch; not for an open pull request, or for a branch the listing gives as having no commits of its own.
---

# Close out a merged pull request

1. Confirm the merge: `gh pr view <n> --json state,mergedAt,headRefName`.
   `state` must be `MERGED`. An open pull request, or one closed without a
   merge, is not closed out: stop and say so.
2. Find where you are: `git branch --show-current` and `git worktree list`.
   - In the main checkout: `git switch main`. Stop if `git status --short`
     prints anything; never stash or discard changes.
   - In a worktree that has the merged branch checked out: the branch cannot
     be deleted from here, and `main` is checked out elsewhere. Do steps 3
     and 5, skip step 4, and say in the report that archiving the session
     removes the worktree and a later close-out from the main checkout
     deletes the branch; the session-start listing will name it.
3. Refresh and prune: `git fetch --prune origin main` drops the stale
   remote-tracking ref. In the main checkout, `git merge --ff-only
   origin/main` then brings `main` to the merge; if that refuses, report it
   and stop.
4. Delete the local branch with `git branch -d <branch>`, never `-D`. Git
   refuses a branch whose commits are not all in `main`, and a branch checked
   out in another worktree; report the refusal with its reason instead of
   forcing the deletion.
5. Remove the pull-request and issue body files written outside the
   repository for this change, by the absolute paths used when they were
   written (AGENTS.md, "How to work here"). Leave anything else there.
6. Report a checklist, each line with the command that showed it:
   - pull request `#n` merged on `<mergedAt>`;
   - local branch deleted, or the reason it is not;
   - `git status --short` prints nothing;
   - `git stash list` prints nothing;
   - the worktree state, from `git worktree list`.

   End with `Safe to archive the session` when every line holds, otherwise
   with what blocks it. Never delete a session or its conversation, and never
   archive one yourself: archiving is reversible and deleting is not, and
   both are the person's action.
