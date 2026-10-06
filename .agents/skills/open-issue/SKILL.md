---
name: open-issue
description: Open an issue that states a problem, its evidence, and the proposed change. Use when asked to file an issue.
---

# Open an issue

1. Search open and closed issues for a duplicate first:
   `gh issue list --state all --search "<terms>"`.
2. Title: follow "Names" in AGENTS.md, "Git".
3. Body: the issue template in `.github/ISSUE_TEMPLATE/` that fits
   (`bug_report.md`, `feature_request.md`, or `question.md`), every section in
   order. Give evidence as `path:line`, `command → result`, or a link, and
   mark wording you drafted as a proposal.
4. Labels: the template's `labels:` value and one or more `area:*` labels.
5. The reason comes from the person who asked, or from the evidence; never
   invent it. Include no secrets, personal data, or unpublished simulation
   data.
6. Open it with
   `gh issue create --title "<title>" --label <type label> --label <area label> --body-file <file>`,
   then give the full chat report.
