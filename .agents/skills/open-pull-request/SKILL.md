---
name: open-pull-request
description: Open a draft GitHub pull request whose body is the full report. Use when a change is committed and ready for review, or to update an open pull request's body; not for reviewing someone else's pull request.
---

# Open a pull request

1. Check the branch, the working tree, and `HEAD`. Run the pull request gate
   in AGENTS.md, "Commands".
2. Before the push, review the diff: `git status --short`, then the whole
   branch diff, `git diff main...HEAD`, for unrelated files, secrets,
   private data, and accidental deletions. A change on the list in AGENTS.md,
   "Ask first", needs that approval before you push it. If the change is
   breaking or deprecates a public name, check that the commit body and the
   pull request body say so (CONTRIBUTING.md, "Changelog").
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
     not list the commits' trailers. Keep each trailer line within 72
     characters, since GitHub wraps the body into the merge commit; give a
     long check result in "Checks run" and a short one in the trailer.

   Then check the trailer block, the body's last paragraph:
   `awk -v RS= 'END { print }' <file> | LC_ALL=C awk 'length > 72'` prints
   nothing, and `git interpret-trailers --parse <file>` prints every trailer.
   Fix the body until both hold. The check counts bytes, and an em dash is
   three bytes, so it can flag a line of 72 characters or fewer; shorten it.
4. Title and labels: follow "Names" in CONTRIBUTING.md. The type comes from
   the branch name; the areas, from the parts the diff changes. `Closes #n`
   goes in the body, never in the title, and the title has no agent or tool
   prefix. If the branch prefix, the title's type, and the `type:` label do
   not all agree, stop and ask the person running you.
5. Before any push, check whether the branch has a pull request:
   `gh pr list --head <branch> --state all --json number,state`. If one has
   merged, stop and tell the person running you. Run this check and the push
   in one command, so the push runs only when the check passes.
6. If the branch has no pull request, push the branch, open a draft, then
   read it back with `gh pr view`:
   `gh pr create --draft --base main --title "<title>" --label <labels> --body-file <file>`.
   Never mark it ready or merge it.
7. If the branch has an open pull request `<n>`, push, check the trailer
   block as in step 3, then replace its body and read it back:
   `gh pr edit <n> --body-file <file>`.
