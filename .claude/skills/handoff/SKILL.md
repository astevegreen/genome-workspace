---
name: handoff
description: Write a handoff note so a fresh session can pick up the work, or resume from one. Use when the owner says handoff, wrap up or resume, when the handoff hook says the conversation is past the handoff line, or before the owner moves to an unrelated task.
argument-hint: "[resume [note-name]]"
---
# Handoff

A long conversation costs more on every turn and attends worse; a 40-line note costs almost nothing to carry.
The note lets the owner `/clear` and carry on in a fresh session, at the PC or through Remote Control.
A Stop hook (`tools/hooks/handoff-nudge.mjs`) asks for one once the conversation passes 60% of the auto-compact window.

`node tools/handoff-notes.mjs` (Windows `node`) keeps the notes in `.handoff/` at the main checkout. Every worktree's
sessions share that folder, and it is git-ignored, so notes are never committed or pushed.
Two branches: **write** (the default) and **resume** (the argument starts with `resume`).

**Notes are personal-data-free.** This repo is public and the notes outlive the session: no genotypes, variants,
findings, report contents, sample names or labels, or file paths that carry a personal name. Name code by
`path:line`, runs by step number ("Step 10 skipped: D: unplugged").

## Write

A note exists only for work still open. **Nothing outstanding means no note**: when the task is finished and nothing waits on the owner or a run, put any future work on `docs/follow-ups.md`, close this session's open notes, and tell the owner in one line that nothing is left to hand off. `save` refuses a note whose **Next** reads `Complete` with no **Decide first** or **Waiting on**.

**Close notes whose work is done.** Before writing a new note, and when this session finishes work an open note lists (`list` shows them), file that note: `node tools/handoff-notes.mjs close <name> <reason>`, with the reason in a few words ("finished in this session", "overtaken by <sha>"). It moves the note to `.handoff/done/` with a `Closed:` line; nothing is deleted. At session start the hook files open notes that leave nothing to do and lists those older than 3 days: when the owner asks, check each against the repo and close the finished ones.

1. **Stop at a green point.** Finish the step in hand until the tests pass in WSL (CLAUDE.md, Commands), or go back to the last step that did and put the rest under **Next**. Background work (a test run, a pipeline step, an agent) finishes first: end the turn saying what it waits on, and write the note when its notification arrives, with the result under **Verification**.
2. **Settle the work.** Commit what is finished, in the normal flow. Work in progress stays uncommitted and is listed by file, with its state, under **Done**; save records which folder it is in. A question whose answer changes what gets committed goes to the owner now; every other owner question goes under **Decide first**. Medical interpretation and scope are always the owner's: carry them, never settle them.
3. **Send durable facts to memory.** A fact that outlives this task (an owner preference, a tooling gotcha) goes to the memory directory; the note carries only what this task needs.
4. **Write the draft** to `.handoff/<slug>.md` (2-4 words naming the task) from the template, at most 40 lines, as a new file even when this session began by resuming a note. Write for a **cold reader**: a session with only the note, the repo and memory. Point at code by `path:line` and at history by sha.
5. **Save it**: `node tools/handoff-notes.mjs save .handoff/<slug>.md`. It writes line 2 (time, branch, HEAD, tree state, folder), checks the note, and files it under a timestamped name. Done when it prints `saved`; on `not saved`, fix what it names and save again.
6. **Tell the owner**, in plain English, in the same turn: the note's name; each **Decide first** question on one line with its recommendation; the **Model** line, if any; and the resume line: type `/clear`, then `/handoff resume`. A quick question can be answered in place; the rest wait for the fresh session.

### Template

    # Handoff: <task in a few words>
    (save writes this line)

    ## Goal
    <the owner's ask, quoted where it was quoted; the PR number; what done looks like>

    ## Decide first
    - <question in plain terms: what changes in the reports or the workflow> — recommend: <option>, because <reason>

    ## Done
    - <sha> <subject>
    - Uncommitted: <path> — <state>   (or: none)

    ## Next
    1. <concrete step, starting from a named file, function or command>

    ## Waiting on
    - <an owner action (an install script, plugging in D:) or a run, and what it unblocks>

    ## Model
    <Opus | Sonnet> at <low | medium | high> effort — <what in Next earns it>

    ## Decisions
    - <decision> — owner: "<quote>"   (or: chosen because <reason>)

    ## Dead ends
    - <what was tried> — <why it failed>

    ## Read first
    - <path:line> — <why>   (five at most)

    ## Verification
    - pytest (WSL): <n> passed on <sha> | failing: <one line> | not run
    - CI: green on PR #<n> | failing: <one line> | not pushed

Keep **Goal** and **Next** always; drop any other section that would be empty. `Complete` only while a **Decide first** or **Waiting on** remains (finished work with an owner answer or a run still pending); with neither, write no note.

**Model**: Opus for anything that decides what a result means (variant logic, a new step's design, tracing a wrong call); Sonnet for implementation with a clear spec. Effort follows the hardest step in **Next**: low for mechanical edits, medium for a specified change, high for tracing a cause or open design. Omit the section when **Next** reads `Complete`.

## Resume

1. **Claim the note**: `node tools/handoff-notes.mjs claim [name]`, with the name from the argument if there is one. With no name it takes the only open note; when several are open it claims nothing and lists them: show the list, ask which, and claim that one. Claiming moves the note to `.handoff/done/`, so no other session takes it, and prints it followed by **Since the note**. `already claimed` means another session took it first: run `list` and ask.
2. **Get to the work**, as **Since the note** says. Run `git switch <branch>` when it says so (uncommitted changes here go to the owner first). When the branch is checked out in another folder, tell the owner to open a session in that folder, and stop. Report commits or changes the note doesn't know about, and re-plan any step they invalidate.
3. **Ask the Decide first questions**, recommended option first, in plain English. The answers are settled for this session, like the note's **Decisions**.
4. **Read the Read-first files**, and only those, before starting.
5. **Restate the goal and the next step** in two lines, re-planned for the answers from step 3. When **Model** names a different model from this session's, say so in one line and let the owner switch before the step starts. A session can't see its own effort level, so never claim it differs: when the note recommends one, name it (`The note recommends high effort.`). Then start.
