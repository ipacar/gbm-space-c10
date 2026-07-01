---
name: feedback-usage-limit-handling
description: User wants long-running autonomous sessions to self-resume across usage-limit resets rather than just stopping
metadata: 
  node_type: memory
  type: feedback
  originSessionId: e2ce4505-eebd-4dd3-8fbe-be5b78067320
---

If a long-running task (e.g. multi-hour GPU job babysitting) risks running into the current session's usage limit, proactively schedule a wakeup (ScheduleWakeup) timed to resume after the limit resets, and leave enough state (TodoWrite list, memory notes, in-progress file edits) that work can continue from where it left off.

**Why:** Explicit user instruction during a long GPU-job-monitoring session ([[project-c10-gbm-space]]) — they don't want the session to just stop and lose momentum when hitting the 5-hour usage window; they want it to pick back up automatically once the limit expires.

**How to apply:** During long autonomous/background-heavy tasks, keep the todo list current and don't leave critical state only in conversation memory. If a usage-limit warning appears, schedule a wakeup for after the reset window instead of ending the task silently.
