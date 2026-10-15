#!/usr/bin/env python3
"""Offline Quantum2 reproduction; never launches Gaussian or contacts a cluster."""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def run(command):
    print("+ " + " ".join(str(v) for v in command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-only", action="store_true", help="Rebuild numeric tables only; no R or image fonts needed")
    parser.add_argument("--verify-raw", action="store_true", help="First verify all 16 .out/.fchk/.gjf bundles against the snapshots")
    parser.add_argument("--orbitals", action="store_true", help="Also reconstruct the orbital PNG from conf-02/uv.fchk (slower)")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable or full path")
    args = parser.parse_args()
    if sys.flags.optimize:
        parser.error("Do not use Python -O; scientific assertions must remain enabled.")
    if args.data_only and args.orbitals:
        parser.error("--orbitals and --data-only cannot be combined")
    if not args.data_only and not shutil.which(args.rscript):
        parser.error("Rscript not found; see ENVIRONMENT.md or use --data-only")
    def python(script):
        run([sys.executable, str(ROOT/"scripts"/script)])
    if args.verify_raw:
        python("verify_public_data.py")
    for script in ["export_tables.py", "compare_properties.py", "compare_experiment.py"]:
        python(script)
    if not args.data_only:
        python("render_stage1_structures.py")
        for script in ["plot_spectra.R", "plot_comparison.R", "plot_experiment.R"]:
            run([args.rscript, str(ROOT/"scripts"/script)])
        if args.orbitals:
            python("render_orbitals.py")
    run([sys.executable, "-m", "unittest", "discover", "-s", "scripts", "-p", "test_experiment.py", "-v"])
    python("check_report.py")
    print("Done. Derived analysis files updated; calculation inputs/results and report source/PDF unchanged.")


if __name__ == "__main__":
    main()
