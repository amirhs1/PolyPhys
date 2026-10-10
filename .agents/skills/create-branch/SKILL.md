---
name: create-branch
description: Create a branch for a change from an up-to-date base. Use when starting a change, before its first commit; not for switching to a branch that exists.
---

# Create a branch

1. Check that the working tree is clean: `git status --short` prints nothing.
   Otherwise stop and ask the person running you; never stash or discard
   their changes.
2. Base: the branch AGENTS.md, "Git", names. Bring it up to date:
   `git fetch origin main`.
3. Name: `<type>/<short-name>`, lowercase and hyphenated, with a type from
   CONTRIBUTING.md, "Names".
4. Create the branch from the fetched base, without tracking it:
   `git switch --no-track -c <type>/<short-name> origin/main`.
5. Check: `git status --short --branch` shows the new branch and nothing else.
