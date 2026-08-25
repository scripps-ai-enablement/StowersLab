# Driving this pipeline from Claude Science or Claude Code

Both let you ask in plain language instead of remembering flags. Neither needs
an administrator or a shared install — this is the path a single user can set up
alone, in a couple of minutes.

The two are different tools with different setups, so they are covered
separately. Pick the one you actually use.

| | Claude Science | Claude Code |
|---|---|---|
| where it runs | the Claude Science web app | your terminal |
| setup | run `install.py` once | clone the repo; the shipped `CLAUDE.md` is picked up automatically |
| how it knows the conventions | the `vno-receptor-rnaseq` skill (`SKILL.md` + `kernel.py`) | `CLAUDE.md` at the repo root |
| helper functions in a kernel | yes (`vno_*`) | no — it drives the CLI |
| good for | interpreting results, figures, follow-up analysis | running it on a cluster, scripting, batch work |

---

## Claude Science

### 1. Install the skill (once per user)

In a Claude Science session, in the **`repl`** tool — not a `python` cell, because
the installer needs the `host` object:

```python
exec(open("/abs/path/to/vno-receptor-rnaseq/install.py").read())
```

Expected output:

```
  SKILL.md: created
  kernel.py: created
  published: published
```

Re-running it after a `git pull` is safe and is how you update — it reports
`already current` when nothing changed.

If you would rather not run a script, `SKILL.md` and `kernel.py` at the repo root
are plain text; write them with `host.skills.edit(...)` and then
`host.skills.publish("vno-receptor-rnaseq")`.

### 2. Load it

```python
skill({skill: "vno-receptor-rnaseq"})
```

The tool result lists the helper functions now defined in your python kernel.
Confirm the environment before spending a run:

```python
vno_install_check()      # what is missing, with the fix for your platform
vno_skill_version()      # which copy you are on
```

### 3. Ask for the analysis

Plain language is enough. A request of this shape works:

> Run the VNO receptor pipeline on `/data/run7/results`, GTF at
> `/refs/GRCm38/genes.gtf`. Tell me which receptor clusters each library
> expresses and flag anything that failed QC.

The agent validates the inputs, runs the stages, and reports by tier. It will not
hand you a receptor name from a library that failed QC, because the pipeline
refuses to produce one.

Useful follow-ups, once results exist:

```python
res = vno_load_results(work_dir="~/vr_run", trial="trial1")
vno_clearance(res["sample_qc"])     # which libraries support any claim
vno_em_verdict(row)                 # read a flag row without re-judging it
vno_flag_vocabulary()               # the exact allowed values per column
```

### What this was tested to do

A plain-language request — two paths in a sentence, no flags, no config, no stage
names — was sent to a session running the specialist profile. It recognised the
existing skill, ran the self-tests (13/13 driver, 35/35 tier gate) before
trusting anything, and submitted the full pipeline to the cluster without being
told any command. So the plain-language entry point genuinely reaches a dispatched
analysis.

**One caveat measured in that same test:** on a busy SLURM cluster the job sat
`PENDING (Resources)` for well over an hour. The pipeline itself takes about a
minute; the queue wait is unbounded and has nothing to do with the agent. Expect
"submitted, waiting" rather than a chat-speed answer when the cluster is loaded.

### The specialist profile

There is also a **VNO Receptor Specialist** profile whose instructions encode the
biology, the five interpretation rules and the tier hierarchy, so it applies them
without being reminded. It is **personal to the account that created it** — a
colleague cannot select it. Making it available to others requires an
administrator to publish it at organization scope, exactly as for the skill. Ask
for both in one message:

> Please publish the skill `vno-receptor-rnaseq` and the agent profile
> `VNO_RECEPTOR_SPECIALIST` at organization scope, so lab members see them
> without installing personal copies. `SKILL.md` and `kernel.py` at the
> repository root are the two files that constitute the skill.

Until then the skill install above gives any session the procedure — you lose the
persona, not the capability.

---

## Claude Code

Claude Code runs in your terminal and drives the pipeline through its
command-line interface. There is no skill to install: clone the repo and the
`CLAUDE.md` at its root is picked up as project context, so the agent starts out
knowing the entry point, the exit-code contract and the interpretation rules.

```bash
git clone <this repo> && cd vno-receptor-rnaseq
claude
```

Then ask:

> Run the VNO receptor pipeline on ~/data/run7/results with the GTF at
> ~/refs/GRCm38/genes.gtf, then tell me which libraries failed QC and why.

### Why the CLI is safe for an agent to drive

Three properties of `vr_analyze.py` matter more here than in interactive use,
because an agent reacts to output automatically rather than reading it:

**`--json-only` puts structured JSON on stdout** and all progress on stderr, so
parsing never has to scrape prose. The JSON carries `libraries[]` (each with
`highest_tier` and `stopped_because`), `clusters_called`, `clusters_suppressed`,
`called_clusters[]`, `candidates`, `top_candidates[]`, `flags{}` and
`unique_read_channel`.

**Exit codes distinguish "no result" from "broken."**

| code | meaning |
|---|---|
| `0` | ran; at least one library reached a cluster-level call |
| `1` | a pipeline stage failed — a genuine error; the log path is printed |
| `3` | ran; **nothing survived QC** — a result, not a crash |
| `4` | sample roles could not be inferred; needs `--target`/`--nontarget` |
| `5` | missing dependency |
| `6` | bad input |

**Check order matters when diagnosing an exit code.** Dependencies are checked
*before* inputs, so on a machine missing `samtools` a bad `--results` path still
reports `5` (missing dependency), not `6`. Resolve `5` first, then re-run — the
input error surfaces after. Verified: the same bad path returns `6` once the
dependency gate passes.

**Never retry on a non-zero exit.** `3` will never succeed on a retry, and an
agent looping on it both wastes cluster time and hides a real finding — the
correct response to `3` is to read the stop reasons and report them.

**`clusters_called` already excludes QC-suppressed libraries.** The withheld
count is reported separately as `clusters_suppressed`, so nothing is silently
dropped. Do not recompute a call count from the raw tables to "check" it; that
bypasses the gate. (The first version of this summary leaked a failed library's
called cluster, which is why `--selftest` now asserts against it.)

### Before a real run

```bash
python3 bin/vr_analyze.py --selftest    # 13 assertions; expect 13 passed, 0 failed
python3 bin/vr_analyze.py --dry-run --results <dir> --gtf <gtf>   # scaffold + validate only
```

`--dry-run` is the cheap way to have an agent confirm it understood the inputs
before committing cluster time.

### On a cluster

Submit rather than running on a login node:

```bash
sbatch --partition=<partition> --time=60 --mem=16G --cpus-per-task=8 \
  --wrap="python3 bin/vr_analyze.py --results <dir> --gtf <gtf> --trial myrun"
```

Not tested from Claude Code specifically — the CLI contract above is what was
tested, and it is the whole interface Claude Code uses.

---

## What neither path changes

Both give you a faster way to ask. Neither loosens the science:

- a wrong-tissue or QC-failed library still yields **zero** receptor statements
- every individual-receptor call still comes back `tentative_unconfirmed`
- `--no-bam` still degrades the result (20 named candidates instead of 25 on the
  reference dataset) and the agent should say so rather than quietly use it

If an agent ever reports a confirmed receptor identity from this pipeline, it is
wrong — no evidence channel here can confirm one at these read lengths.
