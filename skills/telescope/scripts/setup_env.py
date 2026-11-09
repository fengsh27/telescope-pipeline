#!/usr/bin/env python3
"""Step 0 of the telescope skill: make sure the conda environment exists.

Every job script, and check_strandedness.py, starts with

    module load <conda_module>
    conda activate <conda_env>

and dies on the spot when the env is missing. This script creates it from
skills/telescope/environment.yml if it does not exist, then verifies, inside
the env, every tool and Python module the pipeline shells out to or imports.

    python setup_env.py --config telescope.config.json            # create if missing, then verify
    python setup_env.py --config telescope.config.json --check    # verify only, never create
    python setup_env.py --config telescope.config.json --dry-run  # say what would happen

Creating the env takes 10-30 min (conda solve + compiling Telescope's Cython
extension and RSeQC's dependencies) and needs a C compiler on PATH; the
system gcc on OSC is enough. Run it on a login node once per cluster account,
not inside a job.

Telescope's C extension: the historical env was compiled with icc and has
undefined `__intel_*` symbols, which is why the step-2 template LD_PRELOADs
the Intel runtime. A fresh gcc build does not need that; the verification
reports which case the env is in. The preload in the template is harmless
either way, so it stays.

Python 3.6 compatible (runs under the env's interpreter and the system 3.9).
"""

import argparse
import os
import shlex
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (add_config_arg, load_config, conda_prelude, run_bash,  # noqa: E402
                     conda_env_path, ENV_FILE, CONDA_UNAVAILABLE_RC)

# Executables the pipeline runs (worker_patch.py, bam_input.py, strandedness.py).
TOOLS = ("telescope", "bowtie2", "samtools", "trim_galore", "cutadapt", "pigz",
         "kallisto", "infer_experiment.py", "check_strandedness")
# Modules code/lib imports under the env.
MODULES = ("telescope", "pysam", "pandas", "numpy", "scipy", "intervaltree", "yaml")

PROBE = r'''
import importlib, sys
py = sys.version_info
print("PYTHON %%d.%%d.%%d" %% (py.major, py.minor, py.micro))
for m in %(modules)s:
    try:
        mod = importlib.import_module(m)
        print("MODULE %%s OK %%s" %% (m, getattr(mod, "__version__", "")))
    except Exception as e:
        print("MODULE %%s FAIL %%s" %% (m, str(e).replace("\n", " ")[:200]))
# Telescope's Cython extension, loaded the way step 2 checks it.
import ctypes, glob, os
p = glob.glob(os.path.join(os.environ["CONDA_PREFIX"], "lib", "python*",
                           "site-packages", "telescope", "utils", "calignment*.so"))
if not p:
    print("EXT MISSING")
else:
    try:
        ctypes.CDLL(p[0])
        print("EXT OK")
    except OSError as e:
        print("EXT FAIL %%s" %% str(e).replace("\n", " ")[:200])
'''


def activate(cfg):
    return conda_prelude(cfg) + f"conda activate {shlex.quote(cfg['conda_env'])} || exit 96\n"


def intel_preload(cfg):
    """The step-2 template's LD_PRELOAD, for retrying the extension probe."""
    root = cfg["intel_cmplr_root"]
    intel64 = f"{root}/linux/compiler/lib/intel64_lin"
    return (f"module load {cfg['intel_module']} >/dev/null 2>&1 || true\n"
            f'export LD_LIBRARY_PATH="{intel64}:{root}/linux/lib:$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"\n'
            f'export LD_PRELOAD="{intel64}/libimf.so:{intel64}/libsvml.so:{intel64}/libirc.so"\n')


def create_env(cfg, env_file, stream=True):
    target = cfg["conda_env"]
    flag = "-p" if "/" in target else "-n"
    script = (conda_prelude(cfg)
              + f"set -o pipefail\n"
              + f"echo \"creating conda env {target} from {env_file}\"\n"
              + f"conda env create {flag} {shlex.quote(target)} -f {shlex.quote(env_file)}\n")
    return run_bash(script, stream=stream)


def verify(cfg):
    """-> (ok: bool, lines: [str]) after activating the env."""
    lines, ok = [], True
    which = " ".join(f"echo \"TOOL {t} $(command -v {t} || echo MISSING)\";" for t in TOOLS)
    probe = PROBE % {"modules": repr(list(MODULES))}
    script = (activate(cfg)
              + f'echo "PREFIX $CONDA_PREFIX"\n'
              + which + "\n"
              + "python - <<'PY'\n" + probe + "PY\n"
              + "samtools --version 2>&1 | head -1 | sed 's/^/VERSION /'\n"
              + "bowtie2 --version 2>&1 | head -1 | sed 's/^/VERSION bowtie2 /'\n")
    rc, out = run_bash(script, timeout=600)
    if rc == CONDA_UNAVAILABLE_RC:
        return False, ["conda is not available here (module %s)" % cfg["conda_module"]]
    if rc == 96:
        return False, [f"conda activate {cfg['conda_env']} failed"]

    ext_needs_intel = None
    for ln in out.splitlines():
        parts = ln.split(None, 2)
        if not parts:
            continue
        tag = parts[0]
        if tag == "PREFIX":
            lines.append(f"env      : {parts[1] if len(parts) > 1 else '?'}")
        elif tag == "PYTHON":
            v = parts[1]
            good = v.startswith("3.6.")
            ok &= good
            lines.append(f"python   : {v}{'' if good else '  <- expected 3.6.x'}")
        elif tag == "TOOL":
            name, where = parts[1], parts[2] if len(parts) > 2 else "MISSING"
            good = where != "MISSING"
            ok &= good
            lines.append(f"tool     : {name:<22} {'ok  ' + where if good else 'MISSING'}")
        elif tag == "MODULE":
            name, rest = parts[1], parts[2] if len(parts) > 2 else "FAIL"
            good = rest.startswith("OK")
            ok &= good
            lines.append(f"module   : {name:<22} {'ok  ' + rest[3:] if good else rest}")
        elif tag == "EXT":
            if ln.startswith("EXT OK"):
                ext_needs_intel = False
            elif ln.startswith("EXT MISSING"):
                ok = False
                lines.append("telescope C extension: MISSING (install failed?)")
            else:
                ext_needs_intel = True
        elif tag == "VERSION":
            lines.append(f"version  : {ln[8:]}")

    if ext_needs_intel:
        # Same as the step-2 template: does it load once the Intel runtime is preloaded?
        rc2, out2 = run_bash(activate(cfg) + intel_preload(cfg)
                             + "python - <<'PY'\n" + probe + "PY\n", timeout=600)
        if any(l.startswith("EXT OK") for l in out2.splitlines()):
            lines.append("telescope C extension: loads only with the Intel runtime preloaded "
                         "(icc build) -- the step-2 template does that")
        else:
            ok = False
            fail = next((l for l in out2.splitlines() if l.startswith("EXT")), "EXT FAIL")
            lines.append(f"telescope C extension: FAILS to load even with the Intel preload: {fail[9:]}")
    elif ext_needs_intel is False:
        lines.append("telescope C extension: loads without the Intel runtime (gcc build)")
    return ok, lines


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    add_config_arg(ap)
    ap.add_argument("--env-file", default=ENV_FILE,
                    help="environment.yml to create from (default: the skill's)")
    ap.add_argument("--check", action="store_true", help="verify only; never create")
    ap.add_argument("--dry-run", action="store_true", help="report what would be done")
    ap.add_argument("--recreate", action="store_true",
                    help="remove an existing env first (only while no job that uses it is queued)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    name = cfg["conda_env"]
    env_file = os.path.abspath(args.env_file)
    if not os.path.isfile(env_file):
        sys.exit(f"environment file not found: {env_file}")

    path, err = conda_env_path(cfg)
    if err:
        sys.exit(f"cannot reach conda: {err}\n  module load {cfg['conda_module']} and retry")
    print(f"conda module : {cfg['conda_module']}")
    print(f"conda env    : {name}  ({'exists at ' + path if path else 'NOT FOUND'})")
    print(f"env file     : {env_file}")

    if path and args.recreate:
        if args.check:
            sys.exit("--recreate and --check contradict each other")
        print(f"\n--recreate: removing {path}")
        if not args.dry_run:
            rc, out = run_bash(conda_prelude(cfg) + f"conda env remove -y -p {shlex.quote(path)}", stream=True)
            if rc != 0:
                return 1
        path = None

    if not path:
        if args.check:
            print("\nenv missing; --check given, so not creating it. To create:")
            print(f"  python {os.path.abspath(__file__)} --config {args.config}")
            return 1
        if not shutil.which("gcc") and not shutil.which("cc"):
            print("WARNING: no C compiler on PATH; installing telescope-ngs will fail", file=sys.stderr)
        if args.dry_run:
            print("\n--dry-run: would run")
            print(f"  conda env create {'-p' if '/' in name else '-n'} {name} -f {env_file}")
            return 0
        print("\nthis takes 10-30 min (solve + building Telescope's extension) ...\n")
        rc, _ = create_env(cfg, env_file)
        if rc != 0:
            print(f"\nconda env create failed (rc={rc}); see the output above", file=sys.stderr)
            return 1
        path, _ = conda_env_path(cfg)
        if not path:
            print("\nconda env create returned 0 but the env is still not listed", file=sys.stderr)
            return 1
        print(f"\ncreated {path}")
    elif args.dry_run:
        print("\n--dry-run: env exists, would only verify it")
        return 0

    print("\nverifying inside the env ...")
    ok, lines = verify(cfg)
    for ln in lines:
        print("  " + ln)
    if not ok:
        print("\nenv is NOT usable by the pipeline. Fix it, or rebuild with --recreate "
              "(only when nothing that uses it is queued).", file=sys.stderr)
        return 1
    print(f"\nenv {name} is ready. next: check_fastq.py (FASTQ) / check_strandedness.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
