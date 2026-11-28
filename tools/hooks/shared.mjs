// Helpers shared by the handoff tools (tools/handoff-notes.mjs, tools/hooks/handoff-nudge.mjs).
// No dependencies; Node built-ins only.
import { spawnSync } from "node:child_process";
import { existsSync, realpathSync } from "node:fs";
import path from "node:path";

export async function readInput() {
  let raw = "";
  for await (const chunk of process.stdin) raw += chunk;
  try {
    return JSON.parse(raw);
  } catch {
    return {};
  }
}

export function toPosix(p) {
  return p.replace(/\\/g, "/");
}

export function git(args, cwd) {
  const r = spawnSync("git", args, { cwd, encoding: "utf8", windowsHide: true });
  return r.status === 0 ? r.stdout.trim() : null;
}

// Resolve symlinks and Windows short names (C:\Users\RUNNER~1) in the part of `p` that exists,
// so paths compare equal to the long form git reports.
export function canonical(p) {
  let dir = path.resolve(p);
  const rest = [];
  while (!existsSync(dir)) {
    const up = path.dirname(dir);
    if (up === dir) break;
    rest.unshift(path.basename(dir));
    dir = up;
  }
  try { dir = realpathSync.native(dir); } catch {}
  return path.join(dir, ...rest);
}
