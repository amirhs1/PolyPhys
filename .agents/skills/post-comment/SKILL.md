---
name: post-comment
description: Post a short comment on an issue or pull request. Use when asked to reply or comment.
---

# Post a comment

1. Comment only on the issue or pull request the person running you named.
2. Write the comment in this shape:

   ```text
   <Answer in one or two sentences.>
   Based on: <files read or commands run>
   Open: <anything unverified, or None>
   ```

3. Cite evidence as `path:line` or `command → result`. Include no secrets or
   personal data.
4. Post it with `gh issue comment <n> --body-file <file>` or
   `gh pr comment <n> --body-file <file>`, then give the short chat report.
