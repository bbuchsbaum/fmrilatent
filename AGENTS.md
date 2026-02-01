# Repository Guidelines

## Project Structure & Modules
- `R/` – Core S4 classes and generics (`LatentNeuroVec`, accessors) plus helpers in `all_class.R`, `all_generic.R`, `latent_neurovector.R`.
- `man/` – Auto-generated Rd docs (roxygen2). Do not hand-edit.
- HRBF source to port from `~/code/neuroarchive`: analytic helpers in `R/hrbf_helpers.R`, `R/hrbf_core.R`; optional Rcpp kernel in `src/lna_hrbf_rcpp.cpp`. Only the math/algorithms are relevant—ignore neuroarchive-specific pipeline/LNA layers or HDF5 paths.

## Build, Test, and Development
- Load locally: `Rscript -e "devtools::load_all()"`.
- Run tests (once added): `Rscript -e "devtools::test()"`.
- Full check: `R CMD check .` (or `Rscript -e "devtools::check()"`).
- Refresh docs/namespace after new functions: `Rscript -e "devtools::document()"`.

## Coding Style & Naming
- R: 2-space indent, snake_case for functions/objects; S4 methods follow `generic,Class` with registrations in `all_generic.R`.
- Keep matrices as `Matrix` classes where sparse/dense matters; validate dimensions with existing helpers (`check_same_dims`, `validate_same_dims`).
- Keep files ASCII; roxygen comments for user-facing functions. Prefer pure functions; avoid side effects or hidden globals.

## Testing Guidelines
- Framework: testthat (edition 3). Place new specs under `tests/testthat/test-*.R`.
- Use small deterministic masks (`set.seed`) when porting HRBF sampling to keep tests reproducible.
- Aim for round-trip checks: HRBF basis -> coefficients -> `LatentNeuroVec` reconstruction back to dense voxels for a toy mask; include encode/decode and component access paths.
- Add light performance smoke tests (microbenchmark) for basis generation and projection to guard against regressions; prioritize speed and memory use.

## HRBF Porting Plan (from neuroarchive → fmrilatent)
- Phase 1 (pure R port): bring over analytic HRBF generators (`hrbf_basis_from_params`, `generate_hrbf_atom`, `poisson_disk_sample_neuroim2`) adapted to `LatentNeuroVec` inputs. Output: spatial loadings (voxels × atoms) and temporal coefficients (time × atoms) consumable by `LatentNeuroVec`.
- Phase 2 (feature parity): support extra levels and kernel aliases (`gaussian`, `wendland_c6`); keep parameter lists simple and in-R (no descriptor/HDF5 plumbing).
- Phase 3 (performance): optionally port `hrbf_atoms_rcpp` with Eigen/OpenMP for faster basis assembly; guard with an in-package option flag to choose R vs C++.
- Source of truth during port: `R/hrbf_helpers.R` and `src/lna_hrbf_rcpp.cpp` in neuroarchive; keep function names close for easier diffing and testing.

## Commit & PR Guidelines
- Commits: imperative, concise (e.g., “Add HRBF basis generator”, “Wire HRBF to LatentNeuroVec”).
- PRs: describe scope, note tests run, and call out performance implications (sampling determinism, C++ toggle, memory). Link issues when available.

## Issue Tracking with Beads

This project uses **beads** (`bd`) for git-backed issue tracking. See https://github.com/steveyegge/beads

### Essential Commands

| Command | Purpose |
|---------|---------|
| `bd ready` | List tasks without blockers (your next work) |
| `bd create "title" -p 1` | Create task (P0=critical, P1=high, P2=medium, P3=low) |
| `bd show <id>` | View issue details and history |
| `bd update <id> --status in_progress` | Mark task as in progress |
| `bd close <id> --reason "text"` | Close completed task |
| `bd dep add <child> <parent>` | Add dependency |
| `bd list --json` | List all open issues |
| `bd sync` | Force sync to git |

### Critical Rules for Agents

1. **NEVER use `bd edit`** - it opens an interactive editor. Use flag-based updates:
   ```bash
   bd update <id> --description "new description"
   bd update <id> --title "new title"
   ```

2. **Always use `--json` flag** for programmatic access

3. **Run `bd sync` after changes** to ensure immediate git sync

### Finding Work

```bash
bd ready --json          # Tasks without blockers
bd list --status open    # All open tasks
bd stale --days 7        # Neglected tasks
```

## Multi-Agent Coordination (MCP Agent Mail)

This project uses **MCP Agent Mail** for coordination between multiple AI agents. See https://github.com/anthropics/agent-mail

### Getting Started

Before working, register your agent and start a session:

```
mcp__mcp-agent-mail__macro_start_session
  human_key: /Users/bbuchsbaum/code/fmrilatent
  program: claude-code
  model: <your-model>
  task_description: <what you're working on>
```

This returns your agent name (e.g., `GentleMountain`) - use it for all subsequent calls.

### File Reservations (Prevent Conflicts)

**Before editing files**, reserve them to prevent other agents from conflicting edits:

```
mcp__mcp-agent-mail__file_reservation_paths
  project_key: /Users/bbuchsbaum/code/fmrilatent
  agent_name: <your-agent-name>
  paths: ["R/latent_neurovector.R", "R/all_class.R"]
  ttl_seconds: 3600
  exclusive: true
  reason: "Refactoring LatentNeuroVec class"
```

**Release when done:**
```
mcp__mcp-agent-mail__release_file_reservations
  project_key: /Users/bbuchsbaum/code/fmrilatent
  agent_name: <your-agent-name>
```

### Messaging Other Agents

**Send a message:**
```
mcp__mcp-agent-mail__send_message
  project_key: /Users/bbuchsbaum/code/fmrilatent
  sender_name: <your-agent-name>
  to: ["OtherAgentName"]
  subject: "Question about HRBF implementation"
  body_md: "Message content here..."
  thread_id: "HRBF-refactor"  # optional, for threading
```

**Check inbox:**
```
mcp__mcp-agent-mail__fetch_inbox
  project_key: /Users/bbuchsbaum/code/fmrilatent
  agent_name: <your-agent-name>
```

**Acknowledge messages** (especially if `ack_required=true`):
```
mcp__mcp-agent-mail__acknowledge_message
  project_key: /Users/bbuchsbaum/code/fmrilatent
  agent_name: <your-agent-name>
  message_id: <id>
```

### Critical Rules for Agents

1. **Always register first** - Call `macro_start_session` before other operations
2. **Reserve before editing** - Prevents conflicts with other agents
3. **Check inbox periodically** - Other agents may send coordination messages
4. **Release reservations** - When done editing, release your file reservations
5. **Use thread_id** - Keep related discussions organized

### Quick Reference

| Tool | Purpose |
|------|---------|
| `macro_start_session` | Register and boot session |
| `file_reservation_paths` | Reserve files before editing |
| `release_file_reservations` | Release when done |
| `send_message` | Send async message |
| `fetch_inbox` | Check incoming messages |
| `acknowledge_message` | Confirm receipt |

## Landing the Plane (Session Completion)

**When ending a work session**, you MUST complete ALL steps below. Work is NOT complete until `git push` succeeds.

**MANDATORY WORKFLOW:**

1. **File issues for remaining work** - Create issues for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **Release file reservations** - Free any reserved files for other agents:
   ```
   mcp__mcp-agent-mail__release_file_reservations
     project_key: /Users/bbuchsbaum/code/fmrilatent
     agent_name: <your-agent-name>
   ```
5. **PUSH TO REMOTE** - This is MANDATORY:
   ```bash
   git pull --rebase
   bd sync
   git push
   git status  # MUST show "up to date with origin"
   ```
6. **Clean up** - Clear stashes, prune remote branches
7. **Verify** - All changes committed AND pushed
8. **Hand off** - Provide context for next session

**CRITICAL RULES:**
- Work is NOT complete until `git push` succeeds
- NEVER stop before pushing - that leaves work stranded locally
- NEVER say "ready to push when you are" - YOU must push
- If push fails, resolve and retry until it succeeds
