---
name: open-pull-request
description: Open a draft GitHub pull request whose body is the full report. Use when a change is committed and ready for review, or to update an open pull request's body; not for reviewing someone else's pull request.
---

# Open a pull request

1. Check the branch, the working tree, and `HEAD`. Run the full gate in
   AGENTS.md, "Commands".
2. Before the push, review the diff: `git status --short`, then the whole
   branch diff, `git diff main...HEAD`, for unrelated files, secrets,
   private data, and accidental deletions. A change on the list in AGENTS.md,
   "Ask first", needs that approval before you push it.
3. Write the body from `.github/pull_request_template.md`: every section, in
   order; a section that does not apply says `None`.
   - Summary: the reason only from the linked issue, the maintainer, or an
     outside report.
   - Related issues: one `Closes #n` per issue the pull request completes;
     `Refs #n` for one it covers only in part.
   - Checks run: commands you ran in this session, with their actual output;
     tick only the listed checks that ran.
   - Scientific correctness: tick only what you checked, and say why for
     anything left unticked.
   - Notes for review: mark every wording or design you proposed.
   - Last, your own trailer block, as AGENTS.md, "Provenance", gives it. Do
     not list the commits' trailers. Keep each trailer line in the body
     within 72 characters: GitHub wraps the body there when it builds the
     merge commit, and a wrapped line breaks the trailer block.
4. Title and labels: follow "Names" in CONTRIBUTING.md. The type comes from
   the branch name; the areas, from the parts the diff changes. `Closes #n`
   goes in the body, never in the title, and the title has no agent or tool
   prefix. If the branch prefix and the intended type label disagree, stop
   and ask.
5. Push the branch, open a draft, then read it back with `gh pr view`:
   `gh pr create --draft --base main --title "<title>" --label <labels> --body-file <file>`.
   Never mark it ready or merge it.
6. Before pushing to a branch that has a pull request, run
   `gh pr view <n> --json state`; if it has merged, stop and tell the person
   running you.
7. If the branch already has an open pull request, push, then replace its
   body and read it back: `gh pr edit <n> --body-file <file>`.
