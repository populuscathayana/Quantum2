#!/usr/bin/env python3
"""Regenerate the base transition, shielding and coupling tables from snapshots."""
import json
from gaussian_io import ROOT, save_csv


def main():
    paths = [ROOT/"analysis"/name/"verified-results.json" for name in ["stage2","conf01"]]
    paths += sorted((ROOT/"analysis/dimer-properties").glob("*/*/verified-results.json"))
    assert len(paths) == 8
    for path in paths:
        data = json.loads(path.read_text(encoding="utf-8"))
        if "uv_states" in data:
            save_csv(path.parent, "uv-transitions.csv", [{k:v for k,v in s.items() if k != "transitions"} for s in data["uv_states"]])
        if "nmr_shieldings" in data:
            save_csv(path.parent, "nmr-shieldings.csv", data["nmr_shieldings"])
            j = data["nmr_j_hz"]
            save_csv(path.parent, "nmr-couplings.csv", [{"atom_i":i+1,"atom_j":k+1,"J_hz":j[i][k]}
                     for i in range(len(j)) for k in range(i)])
    print("Exported transition, shielding and J-coupling tables for five systems.")


if __name__ == "__main__":
    main()
