#!/usr/bin/env python3
"""Compare audited calculations with traceable experimental reference data."""
import json
from pathlib import Path

import numpy as np

from gaussian_io import save_csv

ROOT = Path(__file__).resolve().parents[1]
DEST = ROOT / "analysis/experimental"


def read_uv_reference():
    """Read the attributed numeric extract distributed with the repository."""
    values = np.genfromtxt(ROOT / "references/experimental/uv-reference-points.csv",
                          delimiter=",", names=True)
    x, y = values["wavelength_nm"], values["log10_epsilon"]
    assert len(x) == 126 and np.all(np.diff(x) > 0)
    assert np.isfinite(x).all() and np.isfinite(y).all()
    return x, y


def read_uv_jcamp(path):
    text = path.read_text(encoding="utf-8")
    assert "##CAS REGISTRY NO=110-44-1" in text
    assert "##XUNITS=Wavelength (nm)" in text
    assert "##YUNITS=Logarithm epsilon" in text
    assert "##XFACTOR=1.0" in text and "##YFACTOR=1.0" in text
    block = text.split("##XYPOINTS=(XY..XY)\n")[1].split("##END=")[0]
    values = np.array([[float(v) for v in line.split(",")] for line in block.splitlines() if line])
    assert values.shape == (126, 2) and np.all(np.diff(values[:, 0]) > 0)
    assert np.isfinite(values).all()
    return values[:, 0], values[:, 1]


def fwhm(x, y):
    """Linear interpolation of both half-height crossings around the maximum."""
    p = int(np.argmax(y))
    level = y[p] / 2
    left = np.flatnonzero(y[:p] < level)
    right = np.flatnonzero(y[p+1:] < level)
    if not len(left) or not len(right):
        return None
    a, b = int(left[-1]), p + 1 + int(right[0])
    low = float(np.interp(level, y[a:a+2], x[a:a+2]))
    high = float(np.interp(level, y[b-1:b+1][::-1], x[b-1:b+1][::-1]))
    return {"left_nm": low, "right_nm": high, "width_nm": high-low}


def compare_protons():
    import csv
    reference = json.loads((ROOT / "references/experimental/nmr-h-reference.json").read_text(encoding="utf-8"))
    shields = list(csv.DictReader((ROOT / "analysis/comparison/nmr-shieldings.csv").read_text(encoding="utf-8").splitlines()))
    anchor = next(p["delta_ppm"] for p in reference["peaks"] if p["site"] == "CH3")
    rows, metrics = [], []
    for p in reference["peaks"]:
        rows.append({"system":"reference","fragment":1,"site":p["site"],
                     "relative_shift_ppm":p["delta_ppm"]-anchor,
                     "reference_relative_shift_ppm":p["delta_ppm"]-anchor,"error_ppm":0.0})
    for name in dict.fromkeys(r["system"] for r in shields):
        nuclei = {int(r["atom"]): float(r["shielding_ppm"]) for r in shields if r["system"] == name}
        units = 1 if name.startswith("conf") else 2
        methyl = float(np.mean([nuclei[i+16*f] for f in range(units) for i in [9,10,11]]))
        # One methyl anchor for the whole species; do not independently align fragments.
        for fragment in range(units):
            for p in reference["peaks"]:
                sigma = float(np.mean([nuclei[i+16*fragment] for i in p["atoms"]]))
                relative = methyl - sigma
                ref = p["delta_ppm"]-anchor
                rows.append({"system":name,"fragment":fragment+1,"site":p["site"],
                             "relative_shift_ppm":relative,"reference_relative_shift_ppm":ref,
                             "error_ppm":relative-ref})
        eligible = [r for r in rows if r["system"]==name and r["site"] in ["H12","H13","H14","H15"]]
        errors = np.array([r["error_ppm"] for r in eligible])
        metrics.append({"system":name,"alkene_site_count":len(errors),
                        "alkene_MAE_ppm":float(np.mean(abs(errors))),
                        "alkene_RMSE_ppm":float(np.sqrt(np.mean(errors**2))),
                        "OH_relative_shifts_ppm":[r["relative_shift_ppm"] for r in rows if r["system"]==name and r["site"]=="OH"]})
    save_csv(DEST,"nmr-H-relative-shifts.csv",rows)
    return {"source":reference,"metrics":metrics,
            "method":"Compare delta_i - delta_methyl with mean(sigma_methyl) - sigma_i; methyl protons averaged as a fast-rotation approximation. Dimers use one mean over both methyl groups but retain individual fragment sites. No empirical scale fit, TMS calculation, population weights, chemical-exchange model or J-split experimental trace reconstruction. MAE excludes methyl anchor and exchangeable OH, and gives each of four alkene sites in each fragment equal weight.",
            "limitations":"Gas-phase fixed-structure predictions versus DMSO-d6 experiment. DMSO hydrogen bonding and exchange can affect OH, so agreement of a cyclic dimer OH spacing does not establish its solution population. Experimental H12/H13 overlap at the reported 0.01 ppm precision; not proof of exact equivalence."}


def compare_carbons():
    import csv
    reference = json.loads((ROOT / "references/experimental/nmr-c-reference.json").read_text(encoding="utf-8"))
    shields = list(csv.DictReader((ROOT / "analysis/comparison/nmr-shieldings.csv").read_text(encoding="utf-8").splitlines()))
    methyl = next(p["delta_ppm"] for p in reference["peaks"] if p["group"] == "CH3")
    ordered = sorted(reference["peaks"], key=lambda p: p["delta_ppm"])
    sites = ["CH3", "alkene-1", "alkene-2", "alkene-3", "alkene-4", "COOH"]
    refs = {site: peak["delta_ppm"]-methyl for site, peak in zip(sites, ordered)}
    rows = [{"system":"reference", "fragment":1, "site":site, "atom":"",
             "source_carbon_number":peak["source_carbon_number"],
             "relative_shift_ppm":refs[site], "reference_relative_shift_ppm":refs[site],
             "error_ppm":0.0} for site, peak in zip(sites, ordered)]
    assigned_refs = {7-p["source_carbon_number"]: p["delta_ppm"]-methyl for p in reference["peaks"]}
    assigned_rows = [{"system":"reference", "fragment":1, "site":f"C{a}", "atom":a,
                      "source_carbon_number":7-a, "relative_shift_ppm":v,
                      "reference_relative_shift_ppm":v, "error_ppm":0.0}
                     for a,v in sorted(assigned_refs.items())]
    metrics = []
    for name in dict.fromkeys(r["system"] for r in shields):
        nuclei = {int(r["atom"]):float(r["shielding_ppm"]) for r in shields
                  if r["system"] == name and r["element"] == "C"}
        units = len(nuclei)//6
        anchor = float(np.mean([nuclei[1+16*f] for f in range(units)]))
        for fragment in range(units):
            offset = 16*fragment
            for site_atom, ref in sorted(assigned_refs.items()):
                atom = site_atom + offset
                relative = anchor-nuclei[atom]
                assigned_rows.append({"system":name,"fragment":fragment+1,"site":f"C{site_atom}",
                                      "atom":atom,"source_carbon_number":7-site_atom,
                                      "relative_shift_ppm":relative,"reference_relative_shift_ppm":ref,
                                      "error_ppm":relative-ref})
            alkene = sorted(range(2+offset, 6+offset), key=lambda a: anchor-nuclei[a])
            # Peak-set matching only, not an experimental atom assignment.
            atoms = [1+offset, *alkene, 6+offset]
            for site, atom in zip(sites, atoms):
                relative = anchor-nuclei[atom]
                rows.append({"system":name,"fragment":fragment+1,"site":site,"atom":atom,
                             "source_carbon_number":"", "relative_shift_ppm":relative,
                             "reference_relative_shift_ppm":refs[site],"error_ppm":relative-refs[site]})
        selected = [r for r in rows if r["system"] == name]
        errors = np.array([r["error_ppm"] for r in selected if r["site"].startswith("alkene-")])
        assigned_errors = np.array([r["error_ppm"] for r in assigned_rows
                                    if r["system"]==name and r["site"] in ["C2","C3","C4","C5"]])
        metrics.append({"system":name, "alkene_peak_count":len(errors),
                        "alkene_assigned_MAE_ppm":float(np.mean(abs(assigned_errors))),
                        "alkene_assigned_RMSE_ppm":float(np.sqrt(np.mean(assigned_errors**2))),
                        "alkene_peak_set_MAE_ppm":float(np.mean(abs(errors))),
                        "alkene_peak_set_RMSE_ppm":float(np.sqrt(np.mean(errors**2))),
                        "carbonyl_relative_shifts_ppm":[r["relative_shift_ppm"] for r in selected if r["site"]=="COOH"],
                        "carbonyl_residuals_ppm":[r["error_ppm"] for r in selected if r["site"]=="COOH"]})
    save_csv(DEST, "nmr-C-relative-shifts.csv", rows)
    save_csv(DEST, "nmr-C-assigned-relative-shifts.csv", assigned_rows)
    return {"source":reference,"metrics":metrics,
            "method":"Experimental delta_i - delta_methyl versus calculated mean(sigma_methyl_C) - sigma_i. Each dimer uses one average of its two methyl carbon shieldings, retaining both fragments. Assigned metrics preserve Table II carbon labels (project carbon = 7 - source carbon). Separate peak-set metrics match ascending positions within each fragment. Both MAEs/RMSEs equally weight the four/eight alkene sites, excluding methyl and carboxyl carbon. No empirical slope fitting or spectral mixture.",
            "limitations":"Gas-phase fixed-structure calculations versus CDCl3 solution; temperature and concentration not specified. Peak-set matching cannot establish site assignments, populations or exchange rates. The source contains peak positions, not a digitized carbon spectrum; equal-height sticks carry no measured intensity information."}


def compare_couplings():
    reference = json.loads((ROOT / "references/experimental/nmr-j-reference.json").read_text(encoding="utf-8"))
    data = json.loads((ROOT / "analysis/stage2/verified-results.json").read_text(encoding="utf-8"))
    rows = []
    for pair in reference["couplings"]:
        i, j = pair["atoms"]
        calculated = data["nmr_j_hz"][i-1][j-1]
        rows.append({"system":"conf-02", "pair":f"H{i}--H{j}", "calculated_J_hz":calculated,
                     "experimental_magnitude_hz":pair["J_hz"], "difference_hz":abs(calculated)-pair["J_hz"]})
    save_csv(DEST, "nmr-J-comparison.csv", rows)
    return {"source":reference,"comparisons":rows}


def build():
    DEST.mkdir(exist_ok=True)
    x, logeps = read_uv_reference()
    eps = 10**logeps
    norm = eps / max(eps)
    plateau = x[logeps == max(logeps)]
    midpoint = float((plateau.min()+plateau.max())/2)
    save_csv(DEST, "uv-reference.csv", [
        {"wavelength_nm": float(a), "log10_epsilon": float(b), "epsilon": float(c), "normalized_intensity": float(d)}
        for a,b,c,d in zip(x,logeps,eps,norm)])
    states = json.loads((ROOT / "analysis/comparison/summary.json").read_text(encoding="utf-8"))
    import csv
    transitions = list(csv.DictReader((ROOT / "analysis/comparison/uv-transitions.csv").read_text(encoding="utf-8").splitlines()))
    curves, peaks = [], []
    grid = np.arange(200.07,292.23001,.01)
    for system in states["systems"]:
        s = [r for r in transitions if r["system"] == system["system"]]
        y = sum(float(r["oscillator_strength"])*np.exp(-4*np.log(2)*((1239.841984/grid-float(r["energy_ev"]))/0.30)**2) for r in s)
        y /= max(y)
        for a,b in zip(grid,y):
            curves.append({"system":system["system"],"wavelength_nm":round(float(a),2),"normalized_intensity":float(b)})
        peak = round(float(grid[np.argmax(y)]),2)
        peaks.append({"system":system["system"],"peak_nm":peak,
                      "difference_from_reference_midpoint_nm":peak-midpoint,"fwhm":fwhm(grid,y)})
    save_csv(DEST,"uv-normalized-calculations.csv",curves)
    result = {"uv":{
        "reference_url":"https://webbook.nist.gov/cgi/cbook.cgi?ID=C110441&Mask=400",
        "download_url":"https://webbook.nist.gov/cgi/cbook.cgi?Index=0&JCAMP=C110441&Type=UVVis",
        "source":"Timmons, C. J., UV atlas of organic compounds, 1966, 1, B7/1; NIST/RAS UV No. 19947",
        "solvent":None,"concentration":None,"temperature":None,
        "n_points":len(x),"range_nm":[float(x.min()),float(x.max())],
        "maximum_log10_epsilon":float(max(logeps)),
        "maximum_plateau_nm":[float(plateau.min()),float(plateau.max())],
        "peak_midpoint_nm":midpoint,"fwhm":fwhm(x,eps),"calculations":peaks,
        "method":"Convert log10(epsilon) to linear epsilon. Normalize each experimental/calculated curve by its own maximum over 200.07-292.23 nm. No horizontal shift, fitted width, baseline subtraction, mixture weights, or extrapolation. Calculations retain energy Gaussian FWHM 0.30 eV.",
        "limitations":"The NIST entry does not state solvent, concentration or temperature; this is not a matched-solvent validation. Plateau midpoint is descriptive, not a measured precision or confidence interval. Experimental spectrum digitized from a historical atlas; FWHM linearly interpolated from tabulated points."
    }}
    result["nmr_H"] = compare_protons()
    result["nmr_C"] = compare_carbons()
    result["nmr_J"] = compare_couplings()
    (DEST/"summary.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n", encoding="utf-8")
    print(f"Saved experimental comparisons: {len(x)} UV points, H/C relative shifts and 2 J couplings.")
    return result


if __name__ == "__main__":
    build()
