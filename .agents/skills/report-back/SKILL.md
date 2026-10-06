---
name: report-back
description: Report back in chat at the end of a task, in the short or the full form. Use at the end of every task and whenever you stop for a decision.
---

# Report back

1. Choose the form, as AGENTS.md, "Reporting", says. Full: this session
   changed a file, opened or updated a pull request or issue, or needs a
   decision from the maintainer. Otherwise short. Posting a comment gets the
   short form.
2. Fill every section; a section that does not apply says `None`. Give the
   verdict once, at the top.
3. Report actual output, not expected output. Say what was not run, and why.
4. List each decision you made that was the maintainer's.
5. In the full form, name under "5 Close-out" the scientific or compatibility
   assumptions the change makes, or say `None`.

Short form:

```text
<Answer in one or two sentences.>
Based on: <files read or commands run; "memory only" if nothing was checked>
Open: <anything unverified, or None>
```

Full form:

```text
## <title>
**Verdict: COMPLETE | NOT COMPLETE — <one line; anything remaining goes here>**
**End product:** <code change | design | issue #n | PR #n | decision for you> — <path or link>

1 What changed — files as path:line, or the issue or PR created
2 Checks run — command → result; anything not run → why
3 Decisions I made that were yours — choice, rejected alternative, cost to reverse
4 What I need from you — Action Needed / Decision Needed, blocking items first; or None
5 Close-out — what to review, branch state, what to keep
```
