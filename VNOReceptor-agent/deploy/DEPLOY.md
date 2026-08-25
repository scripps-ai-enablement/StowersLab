# Deploying so terminal users install nothing

Install the pipeline once in a shared location and colleagues run it as a plain
command — no clone, no `pip install`, no environment to activate.

This requires someone with write access to that shared location. That is
inherent: making software appear for other people is a privileged act, and it is
why per-user installs are the default.

For driving the pipeline from **Claude Science or Claude Code** instead, see
[`AGENT_USAGE.md`](AGENT_USAGE.md) — that path needs no privileges at all.

---

## A shared Environment Modules module

**After setup, a user runs two commands from any directory:**

```bash
module load vno-receptor
vr_analyze.py --results <nfcore results dir> --gtf <genes.gtf>
```

Nothing to clone, no `pip install`, no environment to activate, no config to
edit — `module load` brings the pipeline's python and `samtools` with it.

### How

```bash
PREFIX=/path/to/your/site/software \
bash deploy/install_shared_module.sh vno-receptor-rnaseq-repo.tar.gz
```

Writes two things:

| path | typical permissions |
|---|---|
| `$PREFIX/vno-receptor/1.0.0` | the payload — often group-writable |
| `$PREFIX/modules/vno-receptor/1.0.0` | the modulefile — **often owner-only** |

That split is the usual snag: you can frequently stage the payload yourself but
not the modulefile. The script handles it — it installs the payload, then writes
a ready-to-copy modulefile with the paths already filled in and tells you where
it is, so the request to the modules owner is a one-line copy rather than an
explanation.

Override the dependency module names if your site's differ:

```bash
PREFIX=/path/to/software PYTHON_MODULE=python/3.11 SAMTOOLS_MODULE=samtools \
bash deploy/install_shared_module.sh vno-receptor-rnaseq-repo.tar.gz
```

The pipeline needs pandas, numpy, matplotlib, pyyaml and scipy — a site's
*default* python usually has none of them, which is why the modulefile names an
explicit one. `samtools` is required for the unique-read evidence channel;
without it `vr_analyze.py` stops at its preflight guard, and `--no-bam` runs but
degrades the result.

### Verified

`deploy/install_shared_module.sh` was executed against a scratch prefix on a
SLURM cluster, with `PYTHON_MODULE` and `SAMTOOLS_MODULE` passed explicitly, and
the scratch tree removed afterward. What that run confirmed:

- installer exit 0; 30 files staged; modulefile written
- the generated modulefile had all three site lines substituted and **zero**
  `/path/to/software` placeholders remaining
- `.version` written, so `module load vno-receptor` resolves with no version suffix
- `module load` through the generated modulefile set `VNO_RECEPTOR_HOME`, put
  `vr_analyze.py` on PATH, and pulled in samtools 1.19 and python 3.11.4 /
  pandas 2.2.2
- `module help vno-receptor` printed the usage text
- `vr_analyze.py --selftest` → 13 passed, 0 failed
- a full run invoked as a plain command from an unrelated directory: exit 0,
  18 clusters called, 1 suppressed, 25 named candidates, all
  `tentative_unconfirmed` — matching the reference run

One thing that run exposed, worth knowing if you relocate the files: the script
finds `vno-receptor.modulefile` via its own directory, so **the two must stay
together**. Copying only the `.sh` elsewhere breaks it.

### If your site does not use Environment Modules

Any shared location works. Install the payload somewhere world-readable and add
its `bin/` to `PATH` — a line in a shared profile, a Lmod modulefile, a
`/usr/local/bin` symlink to `vr_analyze.py`, or a container image. The pipeline
has no site assumptions: it resolves its own root from the script location and
reads all thresholds from `config/`.

---

## Until then

The per-user fallback works with no privileges: clone the repo and run
`python3 bin/vr_analyze.py ...`. It has one real hazard — copies drift silently,
and two people can run different conventions with no warning. `kernel.py`
carries `VNO_SKILL_VERSION` (`vno_skill_version()`) so the drift is at least
*detectable*. A single shared install removes it entirely, which is the main
argument for doing one.
