# CLAUDE.md

The agent instructions for this repository are in `AGENTS.md`, which every AI tool shares. It is
imported here:

@AGENTS.md

## Claude Code only

- **Path-scoped rules.** `.claude/rules/` holds short reminders that load when you read or edit
  matching files. They point to the folder READMEs, which hold the design. When a rule and a README
  disagree, the README is right.
- **Permissions.** `.claude/settings.json` makes some commands ask first and refuses others. If a
  command is refused or needs approval, do not reach the same result another way (a different
  command, a script, the web interface): give the person the command to run themselves.
- **Worktrees** go in `.claude/worktrees/<name>`, which git ignores. Remove yours once its branch
  has merged.
- **Subagents** get a self-contained brief and return their findings as text, not as files. Check
  what they report before you repeat it, literature references and "bit-for-bit" claims first.
- **Machine-specific notes** (toolchains, paths, the scheduler) go in `CLAUDE.local.md`, which git
  ignores, or in your user memory, never in a tracked file.
