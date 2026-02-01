# CLAUDE.md - Project Context for Claude Code

## Project Overview

**fmrilatent** is an R package for latent space representations of fMRI neuroimaging data. It provides S4 classes (`LatentNeuroVec`) for efficient storage and manipulation of dimensionality-reduced brain data.

## Quick Commands

```bash
# Development
Rscript -e "devtools::load_all()"      # Load package
Rscript -e "devtools::test()"          # Run tests
Rscript -e "devtools::document()"      # Regenerate docs
Rscript -e "devtools::check()"         # Full R CMD check

# Issue Tracking (beads)
bd ready                               # Find next task
bd create "title" -p 1                 # Create issue (P0-P3)
bd show <id>                           # View issue
bd close <id>                          # Close issue
bd sync                                # Sync to git
```

## Issue Tracking

This project uses **beads** (`bd`) for git-backed issue tracking. See AGENTS.md for full workflow.

**Key commands:**
- `bd ready` - Tasks without blockers (your next work)
- `bd create "title" -p 1` - Create task with priority
- `bd close <id> --reason "text"` - Close completed task
- `bd sync` - Force sync to git

## Multi-Agent Coordination

This project uses **MCP Agent Mail** for multi-agent coordination.

**Project key:** `/Users/bbuchsbaum/code/fmrilatent`

**Before editing files**, reserve them to prevent conflicts:
```
mcp__mcp-agent-mail__file_reservation_paths
  project_key: /Users/bbuchsbaum/code/fmrilatent
  agent_name: <your-agent-name>
  paths: ["R/*.R"]
  ttl_seconds: 3600
  exclusive: true
```

**Check inbox** for messages from other agents:
```
mcp__mcp-agent-mail__fetch_inbox
  project_key: /Users/bbuchsbaum/code/fmrilatent
  agent_name: <your-agent-name>
```

## Session Management

When ending a work session, always "land the plane":

1. File issues for remaining work (`bd create`)
2. Run quality gates (`devtools::test()`, `devtools::check()`)
3. Update issue status (`bd close`, `bd update`)
4. Release file reservations (`release_file_reservations`)
5. Sync and push: `bd sync && git push`
6. Provide handoff context for next session

## Code Style

- 2-space indent, snake_case for functions/objects
- S4 methods follow `generic,Class` pattern
- roxygen2 comments for user-facing functions
- Prefer pure functions; avoid side effects
