---
name: write-commit
description: Write and make a commit in this project's format. Use for every commit.
---

# Write a commit

The rules live in AGENTS.md, "Commit format"; this skill is the procedure.

1. Read the staged diff (`git diff --cached`). One commit holds one coherent
   change.
2. Subject: `<type>(<scope>): <imperative subject>`, with the types and scopes
   in AGENTS.md, "Commit format".
3. Body: bullets of what changed, including which wording or code you were
   given and which you wrote. This is where the detail of your role goes.
4. `Why:` only for a reason from the linked issue, from the maintainer (in the
   pull request or this session), or from an outside report the change
   answers, such as a bug report, a security alert, or a CI failure.
   Otherwise leave it out.
5. End with one trailer block, after a blank line, with no blank line in it
   and nothing after it:
   `Assisted-by: <tool>, <model id or not recorded> (<role>)`, then
   `Checks-run:` for each check you ran, then `Ground-truth-source:` if a
   reference value changed. Pick the role as AGENTS.md, "Commit format",
   defines it.
6. Commit from a file: `git commit -F <message file>`. Never add an AI
   `Co-authored-by:` line, and never use `--no-verify`.
7. Check that git reads every trailer: `git log -1 --format='%(trailers)'`.
