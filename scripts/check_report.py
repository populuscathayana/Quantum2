#!/usr/bin/env python3
"""Check report tables against reviewed raw-result extractions."""
import argparse
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--latex-log", type=Path, help="Optional log from a confirmed current compilation; never assume an older log is current")
args = parser.parse_args()
data = json.loads((ROOT / "analysis/stage2/verified-results.json").read_text(encoding="utf-8"))
text = (ROOT / "experiment-report/sorbic_acid_uv_nmr.tex").read_text(encoding="utf-8")
gs = json.loads((ROOT / "analysis/stage1/audit.json").read_text(encoding="utf-8"))
for row in gs["systems"]:
    expected = (f'{row["id"]} & ${row["E_hartree"]:.9f}$ & '
                f'{row["relative_E_kJmol"]:.4f} & {row["relative_G_kJmol"]:.4f} & {row["min_frequency_cm1"]:.4f}')
    assert expected in text, expected
for state in data["uv_states"]:
    expected = (f'{state["state"]} & {state["energy_ev"]:.4f} & '
                f'{state["wavelength_nm"]:.2f} & {state["oscillator_strength"]:.4f}')
    assert expected in text, expected
for row in data["nmr_shieldings"]:
    atom = f'{row["element"]}{row["atom"]}'
    pattern = rf'{atom}\s*&[^&]+&\s*\$?{re.escape(format(row["shielding_ppm"],".4f"))}\$?'
    assert re.search(pattern,text), atom
for i,j in [(12,13),(14,15),(13,14),(12,14),(13,15),(2,12),(5,15)]:
    ai = next(r["element"] for r in data["nmr_shieldings"] if r["atom"]==i)+str(i)
    aj = next(r["element"] for r in data["nmr_shieldings"] if r["atom"]==j)+str(j)
    value = format(data["nmr_j_hz"][i-1][j-1], ".4f")
    assert re.search(rf'{ai}--{aj}\s*&\s*\$?{re.escape(value)}\$?',text), (ai,aj,value)
for figure in re.findall(r'\\includegraphics\[[^\]]*\]\{([^}]+)\}',text):
    if figure == "image/zju.png" and r"\IfFileExists{image/zju.png}" in text:
        continue  # Optional third-party cover emblem is not distributed.
    assert (ROOT / "experiment-report" / figure).exists(), figure
comparison = json.loads((ROOT / "analysis/comparison/summary.json").read_text(encoding="utf-8"))
for row in comparison["dimers"]:
    label = row["system"].removeprefix("dimer-").replace("-", "/")
    expected = f'{label} & ${row["energy_hartree"]:.9f}$ & {row["relative_energy_kjmol"]:.4f} & {row["min_frequency_cm1"]:.4f}'
    assert expected in text, expected
for row in comparison["systems"]:
    label = row["system"] if row["units"] == 1 else row["system"].removeprefix("dimer-").replace("-", "/")
    expected = f'{label} & {row["bright_state"]} & {row["bright_wavelength_nm"]:.2f} & {row["bright_f_per_molecule"]:.4f} & {row["broadened_peak_nm"]:.2f}'
    assert expected in text, expected
    values = []
    for key in ["OH_shieldings_ppm", "carbonyl_C_shieldings_ppm"]:
        unique = list(dict.fromkeys(format(value,".4f") for value in row[key]))
        values.append(" / ".join(unique))
    assert f'{label} & {values[0]} & {values[1]}' in text, label
experimental = json.loads((ROOT / "analysis/experimental/summary.json").read_text(encoding="utf-8"))
assert "255.56" in text and "256.41" in text
assert f'{experimental["uv"]["fwhm"]["width_nm"]:.1f}' in text
for row in experimental["nmr_H"]["metrics"]:
    label = row["system"] if row["system"].startswith("conf") else row["system"].removeprefix("dimer-").replace("-","/")
    oh = " / ".join(dict.fromkeys(f'{v:.3f}' for v in row["OH_relative_shifts_ppm"]))
    assert f'{label} & {row["alkene_MAE_ppm"]:.3f} & {oh}' in text
assert "DSS-$d_6$" in text and "\\label{fig:uvexp}" in text and "\\label{fig:hexp}" in text
for row in experimental["nmr_C"]["metrics"]:
    label = row["system"] if row["system"].startswith("conf") else row["system"].removeprefix("dimer-").replace("-","/")
    carbonyl = " / ".join(dict.fromkeys(f'{v:.2f}' for v in row["carbonyl_relative_shifts_ppm"]))
    assert f'{label} & {row["alkene_assigned_MAE_ppm"]:.2f} & {row["alkene_peak_set_MAE_ppm"]:.2f} & {carbonyl}' in text
assert "\\label{fig:cexp}" in text and "10.1248/jhs1956.37.89" in text
for value in experimental["nmr_C"]["source"]["peaks"]:
    assert f'{value["delta_ppm"]:.2f}' in text
assert "\\label{tab:cassignment}" in text and "逐碳 MAE" in text and "峰列 MAE" in text
import csv
assigned = list(csv.DictReader((ROOT / "analysis/experimental/nmr-C-assigned-relative-shifts.csv").read_text(encoding="utf-8").splitlines()))
for p in experimental["nmr_C"]["source"]["peaks"]:
    atom = 7-p["source_carbon_number"]
    calc = next(r for r in assigned if r["system"]=="conf-02" and int(r["atom"])==atom)
    assert (f'C{atom} & {p["source_carbon_number"]} & {p["delta_ppm"]:.2f} & '
            f'{float(calc["reference_relative_shift_ppm"]):.2f} & {float(calc["relative_shift_ppm"]):.2f}') in text
for pair in experimental["nmr_J"]["comparisons"]:
    assert f'{pair["pair"]} & {pair["calculated_J_hz"]:.4f} & {pair["experimental_magnitude_hz"]:.1f}' in text
    assert f'{pair["difference_hz"]:.2f}' in text
stack = []
for action, environment in re.findall(r'\\(begin|end)\{([^}]+)\}', text):
    if action == "begin":
        stack.append(environment)
    else:
        assert stack.pop() == environment
assert not stack
labels = re.findall(r'\\label\{([^}]+)\}',text)
assert len(labels)==len(set(labels)), "Duplicate label"
assert set(re.findall(r'\\(?:eqref|ref)\{([^}]+)\}',text)) <= set(labels)
bibkeys = set(re.findall(r'\\bibitem\{([^}]+)\}',text))
assert {key for group in re.findall(r'\\cite\{([^}]+)\}',text) for key in group.split(',')} <= bibkeys
assert "SHA256" not in text and "Mac" not in text and "本次只计算了" not in text
if args.latex_log:
    log = args.latex_log.read_text(encoding="utf-8")
    for forbidden in ("Overfull", "Missing character", "undefined", "! Unable", "Emergency stop"):
        assert forbidden not in log, forbidden
    assert "Output written" in log
print("PASS: baseline data, 3 dimer rows, 5 UV/NMR comparison rows, experimental UV/H/C comparisons, figure references and TeX environments. "
      + ("Supplied compile log passed." if args.latex_log else "PDF compilation/layout not asserted by this source check."))
