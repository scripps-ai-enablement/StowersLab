"""Install this skill into your own Claude Science catalog.

Run in the `repl` tool of a Claude Science session (it needs `host`, which the
python/r kernels do not have):

    exec(open("/path/to/vno-receptor-rnaseq/install.py").read())

Idempotent -- re-run after `git pull` to pick up an update.
Reads SKILL.md and kernel.py from this directory; touches nothing else.
"""
import os

SKILL = "vno-receptor-rnaseq"


def _repo_root():
    """Directory holding SKILL.md. Tries cwd, then common relative spots."""
    for c in (".", SKILL, os.path.join("..", SKILL), os.path.dirname(os.path.abspath("install.py"))):
        if os.path.isfile(os.path.join(c, "SKILL.md")):
            return c
    raise SystemExit(
        "cannot find SKILL.md.\n"
        "  cd into the cloned repo and re-run, or call:\n"
        "     install(repo_dir='/abs/path/to/vno-receptor-rnaseq')"
    )


def install(repo_dir=None):
    d = repo_dir or _repo_root()
    files = {}
    for fn in ("SKILL.md", "kernel.py"):
        p = os.path.join(d, fn)
        if os.path.isfile(p):
            files[fn] = open(p).read()
    if "SKILL.md" not in files:
        raise SystemExit(f"no SKILL.md in {d}")

    present = SKILL in {s["name"] for s in host.skills.list()}
    actions = []
    for fn, content in files.items():
        current = None
        if present:
            try:
                current = host.skills.read(SKILL, path=fn)["content"]
            except Exception:
                current = None
        if current == content:
            actions.append(f"{fn}: already current")
            continue
        r = (host.skills.edit(SKILL, fn, content, old_string=current)
             if current is not None else host.skills.edit(SKILL, fn, content))
        gate = r.get("sidecar_gate")
        if gate and not gate.get("ok"):
            raise SystemExit(f"{fn} rejected by the sidecar gate: {gate.get('error')}")
        actions.append(f"{fn}: {r.get('action', 'written')}")

    pub = host.skills.publish(SKILL, overwrite=True)
    for a in actions:
        print("  " + a)
    print("  published:", pub.get("status"))
    print(f'\nUse it: skill({{skill: "{SKILL}"}})')
    print("Then in a python cell: vno_skill_version(), vno_install_check()")
    return {"actions": actions, "published": pub.get("status"), "source": d}


_result = install()
