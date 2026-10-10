---
name: open-issue
description: Open a GitHub issue that states a problem, its evidence, and the proposed change. Use when asked to file or revise an issue; not for a suspected vulnerability or a reply on an existing issue.
---

# Open an issue

1. A suspected vulnerability goes where `SECURITY.md` says, never into a
   public issue: stop and tell the person running you.
2. Search open and closed issues for a duplicate first:
   `gh issue list --state all --search "<terms>"`. Without network access,
   say that duplication was not checked.
3. Title: follow "Names" in CONTRIBUTING.md.
4. Body: the issue template that fits, every section in order:
   `.github/ISSUE_TEMPLATE/bug_report.md`,
   `.github/ISSUE_TEMPLATE/feature_request.md`, or
   `.github/ISSUE_TEMPLATE/question.md`. Give evidence as `path:line`,
   `command → result`, or a link, and mark wording you drafted as a
   proposal. Last, your trailer block, as AGENTS.md, "Provenance", gives it.
5. The reason comes only from a person or from the source of the evidence;
   never invent it. Include no secrets, personal data, or unpublished
   simulation data.
6. Labels, from "Names": the type that fits the work, with the issue
   template's type as the default, and the area of each part it changes.
7. Open it with
   `gh issue create --title "<title>" --label <labels> --body-file <file>`.
8. To revise an issue you opened, rewrite its body file and run
   `gh issue edit <n> --body-file <file>`.
