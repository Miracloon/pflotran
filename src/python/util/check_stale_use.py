#!/usr/bin/env python
"""Detect likely stale Fortran ``use`` statements in PFLOTRAN sources.

The checker is scope-aware: it inspects each subroutine/function, compares
``use`` imports against symbols referenced in that scope, and tries to avoid
flagging host-associated imports already available from the parent module.

Typical workflow after a refactor:

  ./check_stale_use.py --branch master
  ./check_stale_use.py --files realization_subsurface.F90 init_subsurface.F90

Review findings manually, remove imports, and rebuild.

Git compare modes (two-dot vs three-dot)
----------------------------------------

Given a history like::

  master:      A - B - C - D
                       \\
  your branch:          E - F - G (HEAD)

Three-dot (``master...HEAD``) — default for ``--branch``
  Question: what changed on *my branch* since it split from master?
  Git diffs merge-base ``C`` against ``G``.  Shows only commits/files
  introduced on your branch (``E``, ``F``, ``G``).  Ignores commits that
  landed on master after you branched (``D``).  Does not include uncommitted
  edits.

Two-dot (``master`` or ``master..HEAD``)
  Question: what differs between these two endpoints?
  Git diffs tip of master ``D`` against ``G`` (or the working tree when no
  range is given).  Broader than three-dot: can include master-only changes
  and, for a plain ``git diff master``, uncommitted local edits.

Script behaviour
  ``--branch master`` tries three-dot first, then falls back to two-dot if no
  changed ``.F90`` files are found.  Use ``--two-dot`` to force two-dot;
  use ``--git-diff --base REF --three-dot`` for the same logic without
  ``--branch``.

Author: Glenn Hammond (script generated for PFLOTRAN dev workflow)
"""
from __future__ import print_function

import argparse
import os
import re
import subprocess
import sys

DEFAULT_ROOT = os.path.dirname(os.path.abspath(__file__))

# Whole-module imports that are usually present only for derived-type dummy
# arguments or host typing; suppress unless no related type/name appears.
TYPE_ONLY_MODULE_SUFFIXES = ("_class", "_module")

# Modules whose symbols are often referenced only through dummy-argument typing.
TYPE_ONLY_MODULES = {
    "Realization_Subsurface_class",
    "Realization_Base_class",
    "Simulation_Base_class",
    "Simulation_Subsurface_class",
    "PM_Base_class",
    "Field_module",
    "Patch_module",
    "Grid_module",
    "Option_module",
    "Discretization_module",
}


def run_git(args, cwd):
    """Run a git command and return stripped stdout."""
    proc = subprocess.Popen(
        ["git"] + args,
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(
            "git {} failed: {}".format(" ".join(args), err.decode().strip())
        )
    return out.decode().strip()


def git_repo_root(start):
    """Return git repository root containing ``start``."""
    return run_git(["rev-parse", "--show-toplevel"], start)


def resolve_git_ref(repo, ref):
    """Resolve a branch or ref name to a verified committish."""
    candidates = [ref]
    if "/" not in ref:
        candidates.append("origin/{}".format(ref))
    last_err = None
    for candidate in candidates:
        proc = subprocess.Popen(
            ["git", "rev-parse", "--verify", candidate],
            cwd=repo,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        out, err = proc.communicate()
        if proc.returncode == 0:
            return candidate, out.decode().strip()
        last_err = err.decode().strip()
    raise SystemExit(
        "Unknown git ref '{}'. {}".format(ref, last_err or "")
    )


def git_changed_f90_files(repo, rel_root, compare_ref, three_dot):
    """Return .F90 paths changed relative to ``compare_ref``."""
    if three_dot:
        rev_range = "{}...HEAD".format(compare_ref)
    else:
        rev_range = compare_ref
    names = run_git(
        ["diff", "--name-only", rev_range, "--", rel_root], repo
    )
    files = []
    for name in names.splitlines():
        name = name.strip()
        if name.endswith(".F90"):
            files.append(os.path.join(repo, name))
    return sorted(files)


def collect_f90_files(
    root, explicit_files, git_diff, git_base, git_branch, three_dot, scan_all
):
    """Collect .F90 files to scan."""
    if scan_all:
        files = []
        for dirpath, _, fnames in os.walk(root):
            for fn in fnames:
                if fn.endswith(".F90"):
                    files.append(os.path.join(dirpath, fn))
        return sorted(files), None

    if explicit_files:
        out = []
        for path in explicit_files:
            if not os.path.isabs(path):
                path = os.path.join(root, path)
            if os.path.isfile(path):
                out.append(path)
            else:
                print("WARNING: file not found: {}".format(path), file=sys.stderr)
        return sorted(out), None

    if git_diff or git_branch:
        repo = git_repo_root(root)
        rel_root = os.path.relpath(root, repo)
        compare_ref = git_branch if git_branch else git_base
        resolved_ref, _ = resolve_git_ref(repo, compare_ref)
        if git_branch:
            used_three_dot = three_dot
            files = git_changed_f90_files(
                repo, rel_root, resolved_ref, three_dot=used_three_dot
            )
            if not files and used_three_dot:
                files = git_changed_f90_files(
                    repo, rel_root, resolved_ref, three_dot=False
                )
                used_three_dot = False
            compare_label = (
                "{}...HEAD".format(resolved_ref)
                if used_three_dot
                else resolved_ref
            )
        else:
            files = git_changed_f90_files(
                repo, rel_root, resolved_ref, three_dot=three_dot
            )
            compare_label = (
                "{}...HEAD".format(resolved_ref)
                if three_dot
                else resolved_ref
            )
        return files, compare_label

    raise SystemExit(
        "Specify --branch, --git-diff, --files, or --all.\n"
        "Examples:\n"
        "  ./check_stale_use.py --branch master\n"
        "  ./check_stale_use.py --git-diff --base HEAD"
    )


def read_text(path):
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        return f.read()


def join_continued_lines(lines):
    """Join Fortran continuation lines ending with ``&``."""
    joined = []
    buf = ""
    for line in lines:
        stripped = line.rstrip()
        if not buf:
            buf = stripped
        else:
            buf += " " + stripped.lstrip()
        if "&" not in stripped.rstrip():
            joined.append(buf.replace("&", " ").strip())
            buf = ""
    if buf:
        joined.append(buf.replace("&", " ").strip())
    return joined


def parse_public_symbols(module_text):
    """Extract public symbols from a Fortran module file."""
    symbols = set()
    lines = join_continued_lines(module_text.splitlines())

    for line in lines:
        m = re.match(r"^\s*public\s*::\s*(.+)$", line, re.I)
        if m:
            chunk = m.group(1)
            for part in chunk.split(","):
                part = part.strip()
                part = re.sub(r"=>.*", "", part).strip()
                if part and part != "&":
                    symbols.add(part)

    for m in re.finditer(
        r"^\s*type\s*,\s*public(?:\s*,\s*extends\s*\([^)]+\))?\s*::\s*(\w+)",
        module_text,
        re.I | re.M,
    ):
        symbols.add(m.group(1))

    for m in re.finditer(
        r"^\s*(?:integer|real|logical|character|class)\s*,\s*"
        r"parameter\s*,\s*public\s*::\s*(\w+)",
        module_text,
        re.I | re.M,
    ):
        symbols.add(m.group(1))

    for m in re.finditer(
        r"^\s*(subroutine|function)\s+(\w+)",
        module_text,
        re.I | re.M,
    ):
        symbols.add(m.group(2))

    return symbols


def build_module_index(root):
    """Map module name -> metadata used for stale-import detection."""
    index = {}
    for dirpath, _, fnames in os.walk(root):
        for fn in fnames:
            if not fn.endswith(".F90"):
                continue
            path = os.path.join(dirpath, fn)
            text = read_text(path)
            m = re.search(r"^\s*module\s+(\w+)\s*$", text, re.I | re.M)
            if not m:
                continue
            mod = m.group(1)
            symbols = parse_public_symbols(text)
            type_names = {s for s in symbols if s.endswith("_type")}
            prefixes = set()
            for sym in symbols:
                m2 = re.match(r"([A-Z][a-zA-Z0-9]*)", sym)
                if m2 and len(m2.group(1)) > 2:
                    prefixes.add(m2.group(1))
            index[mod] = {
                "path": path,
                "symbols": symbols,
                "type_names": type_names,
                "prefixes": sorted(prefixes, key=len, reverse=True),
            }
    return index


def parse_use_statement(line):
    """Parse one ``use`` statement line (already continuation-joined)."""
    line = line.strip()
    if not re.match(r"use\b", line, re.I):
        return None
    if re.match(r"use\s*,\s*(intrinsic|non_intrinsic)\b", line, re.I):
        return None

    m = re.match(
        r"use\s+(?:,\s*[^\s]+\s*::\s*)?(\w+)(?:\s*,\s*only\s*:\s*(.+))?",
        line,
        re.I,
    )
    if not m:
        return None

    mod = m.group(1)
    only_part = m.group(2)
    only_names = []
    if only_part:
        for part in only_part.split(","):
            part = part.strip()
            part = re.sub(r"=>.*", "", part).strip()
            if part:
                only_names.append(part)
    return mod, only_names


def symbol_used(sym, body):
    """Return True if ``sym`` appears to be referenced in ``body``."""
    return bool(re.search(r"\b" + re.escape(sym) + r"\b", body, re.I))


def type_only_module_used(mod, meta, body):
    """Heuristic for imports used only through derived-type declarations."""
    if mod in TYPE_ONLY_MODULES:
        for tname in meta.get("type_names", ()):
            if symbol_used(tname, body):
                return True
        # e.g. Realization_Subsurface_class -> realization_subsurface_type
        guess = mod.replace("_class", "_type")
        if symbol_used(guess, body):
            return True
        guess = re.sub(
            r"_module$", "_type", mod, flags=re.I
        )
        if guess != mod and symbol_used(guess, body):
            return True
    if mod.endswith(TYPE_ONLY_MODULE_SUFFIXES):
        for tname in meta.get("type_names", ()):
            if symbol_used(tname, body):
                return True
    return False


def select_type_names_used(body):
    """Return type names referenced in ``select type`` / ``class is`` blocks."""
    names = set()
    for m in re.finditer(
        r"class\s+is\s*\(\s*(\w+)", body, re.I
    ):
        names.add(m.group(1))
    for m in re.finditer(
        r"type\s+is\s*\(\s*(\w+)", body, re.I
    ):
        names.add(m.group(1))
    return names


def module_symbols_used(mod, meta, body):
    """Return True if any exported symbol from ``mod`` appears in ``body``."""
    if not meta:
        return False

    for sym in meta.get("symbols", ()):
        if symbol_used(sym, body):
            return True

    for tname in select_type_names_used(body):
        if tname in meta.get("type_names", ()):
            return True
        if tname in meta.get("symbols", ()):
            return True

    for prefix in meta.get("prefixes", ()):
        if re.search(r"\b" + re.escape(prefix) + r"\w*\b", body):
            return True

    return type_only_module_used(mod, meta, body)


def collect_module_uses(lines, start, end):
    """Collect ``use`` statements between ``start`` and ``end`` line indices."""
    uses = []
    i = start
    while i < end:
        stripped = lines[i].strip()
        if re.match(r"use\b", stripped, re.I):
            chunk = [stripped]
            j = i + 1
            while j < end and "&" in lines[j - 1].rstrip():
                chunk.append(lines[j].strip())
                j += 1
            joined = join_continued_lines(chunk)
            if joined:
                uses.append(joined[0])
            i = j
            continue
        i += 1
    return uses


def parse_scopes(text):
    """Yield scope dicts for subroutines/functions in a source file."""
    lines = text.splitlines()
    in_module = False
    module_name = None
    module_start = 0
    module_uses = []

    i = 0
    while i < len(lines):
        stripped = lines[i].strip()
        low = stripped.lower()

        if re.match(r"module\s+\w+", low) and "procedure" not in low:
            in_module = True
            module_name = re.match(r"module\s+(\w+)", low, re.I).group(1)
            module_start = i
            module_uses = []

        if in_module and re.match(r"use\b", stripped, re.I):
            chunk = [stripped]
            j = i + 1
            while j < len(lines) and "&" in lines[j - 1].rstrip():
                chunk.append(lines[j].strip())
                j += 1
            joined = join_continued_lines(chunk)
            if joined:
                module_uses.append(joined[0])
            i = j
            continue

        m = re.match(r"(subroutine|function)\s+(\w+)", low)
        if m:
            kind = m.group(1)
            name = m.group(2)
            start = i
            i += 1
            while i < len(lines):
                low2 = lines[i].strip().lower()
                if re.match(r"(subroutine|function)\s+\w+", low2):
                    break
                if re.match(r"end\s+(subroutine|function)", low2):
                    break
                i += 1
            end = i
            scope_uses = collect_module_uses(lines, start + 1, end)
            body = "\n".join(lines[start : end + 1])
            yield {
                "kind": kind,
                "name": name,
                "start_line": start + 1,
                "module": module_name,
                "module_uses": list(module_uses),
                "scope_uses": scope_uses,
                "body": body,
            }
            continue

        if re.match(r"end\s+module", low):
            in_module = False
            module_name = None
            module_uses = []

        i += 1


def host_provides_module(module_uses, mod):
    """Return True if parent module already imports ``mod``."""
    for u in module_uses:
        parsed = parse_use_statement(u)
        if parsed and parsed[0] == mod:
            return True
    return False


def check_scope(scope, module_index, verbose):
    """Return finding dicts for one subroutine/function scope."""
    findings = []
    body = scope["body"]

    for use_line in scope["scope_uses"]:
        parsed = parse_use_statement(use_line)
        if not parsed:
            continue
        mod, only_names = parsed
        meta = module_index.get(mod, {})

        if only_names:
            for sym in only_names:
                if not symbol_used(sym, body):
                    findings.append({
                        "line": scope["start_line"],
                        "scope": scope["name"],
                        "module": mod,
                        "symbol": sym,
                        "kind": "unused-only",
                        "confidence": "high",
                        "detail": (
                            "symbol '{}' from 'only:' list not referenced "
                            "in this scope".format(sym)
                        ),
                    })
            continue

        if module_symbols_used(mod, meta, body):
            continue

        if host_provides_module(scope["module_uses"], mod):
            confidence = "medium"
            detail = (
                "no direct references found; parent module already "
                "imports '{}' (likely redundant)".format(mod)
            )
        else:
            confidence = "medium"
            detail = "no exported symbols from '{}' found in this scope".format(
                mod
            )

        if mod in TYPE_ONLY_MODULES and type_only_module_used(mod, meta, body):
            continue

        findings.append({
            "line": scope["start_line"],
            "scope": scope["name"],
            "module": mod,
            "symbol": None,
            "kind": "stale-whole-module",
            "confidence": confidence,
            "detail": detail,
        })

    if verbose:
        print(
            "  scanned {} {} (line {})".format(
                scope["kind"], scope["name"], scope["start_line"]
            ),
            file=sys.stderr,
        )
    return findings


def format_finding(path, finding):
    """Format one finding for terminal output."""
    loc = "{}:{} in {}".format(
        os.path.basename(path),
        finding["line"],
        finding["scope"],
    )
    if finding["symbol"]:
        target = "{} :: {}".format(finding["module"], finding["symbol"])
    else:
        target = finding["module"]
    return "[{}:{}] {} -> {} ({})".format(
        finding["confidence"],
        finding["kind"],
        loc,
        target,
        finding["detail"],
    )


def main():
    parser = argparse.ArgumentParser(
        description="Detect likely stale Fortran use statements in PFLOTRAN.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap_epilog(),
    )
    parser.add_argument(
        "--root",
        default=DEFAULT_ROOT,
        help="PFLOTRAN Fortran source root (default: script directory)",
    )
    parser.add_argument(
        "--branch",
        metavar="NAME",
        help=(
            "Scan .F90 files changed on the current branch since diverging "
            "from branch NAME (e.g. master). Implies --git-diff."
        ),
    )
    parser.add_argument(
        "--git-diff",
        action="store_true",
        help="Scan .F90 files changed vs a git ref (--base or --branch)",
    )
    parser.add_argument(
        "--base",
        default="HEAD",
        help=(
            "Git ref for --git-diff (default: HEAD). Accepts commits, tags, "
            "and branch names. Use --branch for the common master/compare case."
        ),
    )
    parser.add_argument(
        "--three-dot",
        action="store_true",
        default=None,
        help=(
            "With --git-diff/--base, use REF...HEAD (changes on current "
            "branch). Default for --branch."
        ),
    )
    parser.add_argument(
        "--two-dot",
        action="store_true",
        help=(
            "With --git-diff/--base, compare working tree to REF (not "
            "three-dot). Ignored when --branch is set unless no commits differ."
        ),
    )
    parser.add_argument(
        "--files",
        nargs="*",
        metavar="FILE",
        help="Explicit .F90 files to scan",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Scan all .F90 files under --root (noisy)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print scanned scopes to stderr",
    )
    args = parser.parse_args()

    if args.branch and args.git_diff and args.base != "HEAD":
        print(
            "NOTE: --branch {} overrides --base {}".format(
                args.branch, args.base
            ),
            file=sys.stderr,
        )

    if args.two_dot and args.three_dot:
        raise SystemExit("Specify only one of --two-dot or --three-dot.")

    if args.branch:
        three_dot = not args.two_dot
    elif args.git_diff:
        three_dot = args.three_dot if args.three_dot is not None else False
    else:
        three_dot = False

    root = os.path.abspath(args.root)
    files, compare_label = collect_f90_files(
        root,
        args.files,
        args.git_diff or bool(args.branch),
        args.base,
        args.branch,
        three_dot,
        args.all,
    )
    if not files:
        if compare_label:
            print(
                "No changed .F90 files under {} (git diff {}).".format(
                    root, compare_label
                )
            )
        else:
            print("No .F90 files to scan.")
        return 0

    if compare_label:
        print(
            "Git compare: {} ({} file(s))".format(compare_label, len(files)),
            file=sys.stderr,
        )

    if args.verbose:
        print("Building module index from {} ...".format(root), file=sys.stderr)
    module_index = build_module_index(root)

    conf_rank = {"high": 0, "medium": 1, "low": 2}
    findings_by_file = {}
    if args.verbose:
        print("Scanning {} file(s)...".format(len(files)), file=sys.stderr)

    for path in files:
        text = read_text(path)
        file_findings = []
        for scope in parse_scopes(text):
            file_findings.extend(
                check_scope(scope, module_index, args.verbose)
            )
        if file_findings:
            findings_by_file[path] = file_findings

    total = 0
    for path in files:
        file_findings = findings_by_file.get(path, [])
        if not file_findings:
            continue
        file_findings.sort(
            key=lambda f: (
                conf_rank.get(f["confidence"], 9),
                f["line"],
                f["scope"],
            )
        )
        print(os.path.relpath(path, root))
        for finding in file_findings:
            print("  " + format_finding(path, finding))
            total += 1
        print("")

    if total == 0:
        print("No likely stale use statements found.")
    else:
        print(
            "Found {} finding(s). Review manually and rebuild after removals.".format(
                total
            )
        )
        print(
            "Tip: prioritize [high:unused-only] and recent-refactor scopes first.",
            file=sys.stderr,
        )
        return 1
    return 0


def textwrap_epilog():
    return """
Examples:
  %(prog)s --branch master
  %(prog)s --branch develop
  %(prog)s --branch master --two-dot
  %(prog)s --git-diff
  %(prog)s --git-diff --base origin/main --three-dot
  %(prog)s --files realization_subsurface.F90 init_subsurface.F90
  %(prog)s --all --verbose

Git compare: two-dot vs three-dot
---------------------------------
History::

  master:      A - B - C - D
                       \\
  your branch:          E - F - G (HEAD)

Three-dot: master...HEAD  (default for --branch)
  - Question: what changed on MY branch since it forked from master?
  - Compares merge-base C to HEAD G (commits E, F, G only).
  - Ignores commits on master after you branched (D).
  - Commit-based; uncommitted edits are not included.
  - Best for reviewing a feature branch / PR vs master.

Two-dot: master  (or --two-dot / fallback after --branch)
  - Question: what is different between master and here?
  - Compares tip of master D to HEAD G (or working tree vs master).
  - Can include master-only changes and local uncommitted edits.
  - Use when you want everything that differs from master right now.

Side-by-side:
  three-dot master...HEAD     two-dot master
  your branch commits only    full endpoint diff + dirty tree
  default: --branch master    force: --branch master --two-dot

This script:
  - --branch NAME resolves NAME (then origin/NAME), tries NAME...HEAD first.
  - If no changed .F90 files, falls back to two-dot vs NAME.
  - --git-diff --base REF --three-dot gives three-dot without --branch.

Notes:
  - Findings are heuristic; always rebuild after removing imports.
  - Whole-module stale detection is medium confidence by design.
  - use ..., only : unused symbols are high confidence.
"""


if __name__ == "__main__":
    sys.exit(main())