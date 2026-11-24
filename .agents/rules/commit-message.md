---
trigger: always_on
---

# Commit Message Standard

Whenever generating a commit message or running `git commit`, inspect the effective commit diff — including staged changes (`git diff --staged`), tracked modifications included by `git commit -a` (`git diff HEAD`), and existing commit changes retained when amending (`git show` or `git diff HEAD~1`) — and follow both the **Conventional Commits** specification and the **50/72** rule:

1. **Conventional Commits Format**
   - Header: `<type>(<optional scope>): <subject>` (append `!` before `:` for breaking changes, e.g. `feat!: ...`).
   - Allowed types: `feat`, `fix`, `docs`, `style`, `refactor`, `perf`, `test`, `build`, `ci`, `chore`, `revert`.

2. **50/72 Formatting**
   - **Subject line (header):** Maximum **50 characters** total (including `<type>(<scope>): `). Use imperative mood (e.g., `add`, not `added` or `adds`) and do not end with a period.
   - **Blank line:** Leave a single blank line between the subject line and the body.
   - **Body:** Wrap every line at **72 characters**. Explain *why* the change was made rather than *what* changed in the diff.
   - **Footer (optional):** Include `BREAKING CHANGE: ...` or issue references when relevant, wrapped at 72 characters.
