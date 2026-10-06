---
name: open-pull-request
description: Push the branch and open a draft pull request whose body is the full report. Use when a change is ready for review.
---

# Open a pull request

1. Check the branch, the working tree, and `HEAD`. Run the full gate in
   AGENTS.md, "Commands".
2. Write the body from `.github/pull_request_template.md`: every section, in
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
   - AI assistance, last: tool, model, role, then the branch's `Assisted-by:`
     lines from
     `git log --no-merges --format=%B main..HEAD | grep '^Assisted-by:'`.
3. Title: follow "Names" in AGENTS.md, "Git". Labels: the type label the
   branch prefix sets and one or more `area:*` labels, as "Git" says; if the
   prefix and the intended label disagree, stop and ask.
4. Push the branch and open a draft:
   `gh pr create --draft --base main --title "<title>" --label <type label> --label <area label> --body-file <file>`.
   Never mark it ready or merge it.
5. Read the body back (`gh pr view`), then give the full chat report.
