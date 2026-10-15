#!/usr/bin/env python3
"""Assemble audited monomer/dimer comparisons; no solution population model."""
import json
from pathlib import Path

import numpy as np

from gaussian_io import BONDS, ROOT, file_hash, save_csv


def build():
    destination = ROOT / "analysis/comparison"
    destination.mkdir(exist_ok=True)
    states, shields, curves, summaries, geometries, bonds, dimers = [], [], [], [], [], [], []
    provenance = {}
    names = ["conf-01", "conf-02", "dimer-01-01", "dimer-01-02", "dimer-02-02"]
    optimized = {name: json.loads((ROOT / "analysis/dimers" / name / "audit.json").read_text(encoding="utf-8")) for name in names[2:]}
    minimum = min(d["E_hartree"] for d in optimized.values())
    for name in names:
        units = 1 if name.startswith("conf") else 2
        if units == 1:
            uvpath = ROOT / ("analysis/conf01" if name == "conf-01" else "analysis/stage2") / "verified-results.json"
            nmrpath = uvpath
        else:
            uvpath = ROOT / "analysis/dimer-properties" / name / "uv/verified-results.json"
            nmrpath = uvpath.parent.parent / "nmr/verified-results.json"
        uv, nmr = json.loads(uvpath.read_text(encoding="utf-8")), json.loads(nmrpath.read_text(encoding="utf-8"))
        for path in [uvpath, nmrpath]:
            provenance[str(path.relative_to(ROOT))] = file_hash(path)
        assert len(uv["uv_states"]) == 20 * units and len(nmr["nmr_shieldings"]) == 16 * units
        for state in uv["uv_states"]:
            states.append({"system": name, "units": units, **{k: v for k, v in state.items() if k != "transitions"}})
        for row in nmr["nmr_shieldings"]:
            shields.append({"system": name, "units": units, **row})
        wavelengths = np.arange(160, 350.0001, .01)
        energy = 1239.841984 / wavelengths
        signal = sum(s["oscillator_strength"] * np.exp(-4*np.log(2)*((energy-s["energy_ev"])/.30)**2)
                     for s in uv["uv_states"]) / units
        assert 1239.841984 / 160 + .90 < uv["uv_states"][-1]["energy_ev"]
        assert all(s["energy_ev"] > 1239.841984 / 350 for s in uv["uv_states"])
        for x, y in zip(wavelengths, signal):
            curves.append({"system": name, "wavelength_nm": round(float(x), 2), "amplitude_per_acid_unit": float(y)})
        bright = max(uv["uv_states"], key=lambda s: s["oscillator_strength"])
        oh = [s["shielding_ppm"] for s in nmr["nmr_shieldings"] if s["atom"] in [16, 32]]
        carbonyl = [s["shielding_ppm"] for s in nmr["nmr_shieldings"] if s["atom"] in [6, 22]]
        summaries.append({"system": name, "units": units, "bright_state": bright["state"],
                          "bright_energy_ev": bright["energy_ev"], "bright_wavelength_nm": bright["wavelength_nm"],
                          "bright_f_per_molecule": bright["oscillator_strength"],
                          "broadened_peak_nm": round(float(wavelengths[np.argmax(signal)]), 2),
                          "highest_excitation_ev": uv["uv_states"][-1]["energy_ev"],
                          "OH_shieldings_ppm": oh, "carbonyl_C_shieldings_ppm": carbonyl,
                          "alkene_J_hz": [nmr["nmr_j_hz"][i][j] for i, j in ([(11, 12), (13, 14)] if units == 1 else [(11, 12), (13, 14), (27, 28), (29, 30)])]})
        if units == 2:
            d = optimized[name]
            assert d["verified"] and d["imaginary_count"] == 0 and not d["same_minimum_as"]
            dimers.append({"system": name, "energy_hartree": d["E_hartree"],
                           "relative_energy_kjmol": (d["E_hartree"]-minimum)*2625.4996394799,
                           "min_frequency_cm1": d["min_frequency_cm1"],
                           "OO_A": [h["O_O_angstrom"] for h in d["hydrogen_bonds"]],
                           "HO_A": [h["H_O_angstrom"] for h in d["hydrogen_bonds"]]})
            xyz = np.array(d["coordinates_angstrom"])
            xyz -= xyz.mean(axis=0)
            heavy = [i for i, number in enumerate(d["atom_numbers"]) if number != 1]
            _, _, axes = np.linalg.svd(xyz[heavy] - xyz[heavy].mean(axis=0), full_matrices=False)
            xyz = xyz @ axes.T
            if xyz[0, 0] > 0:
                xyz[:, 0] *= -1
            if xyz[7, 1] < 0:
                xyz[:, 1] *= -1
            for i, (number, point) in enumerate(zip(d["atom_numbers"], xyz)):
                geometries.append({"system": name, "atom": i+1, "element": {1:"H",6:"C",8:"O"}[number],
                                   "x": float(point[0]), "y": float(point[1])})
            for offset in [0, 16]:
                for a, b, order in BONDS:
                    a, b = a-1+offset, b-1+offset
                    bonds.append({"system": name,"x":float(xyz[a,0]),"y":float(xyz[a,1]),"xend":float(xyz[b,0]),"yend":float(xyz[b,1]),"kind":"covalent"})
            for a,b in [(15,22),(31,6)]:
                bonds.append({"system": name,"x":float(xyz[a,0]),"y":float(xyz[a,1]),"xend":float(xyz[b,0]),"yend":float(xyz[b,1]),"kind":"hydrogen"})
    for filename, rows in [("uv-transitions.csv", states),("nmr-shieldings.csv", shields),("uv-curves.csv", curves),
                           ("dimer-atoms.csv", geometries),("dimer-bonds.csv", bonds)]:
        save_csv(destination, filename, rows)
    result = {"systems": summaries, "dimers": dimers, "provenance": provenance,
              "uv_method": "Energy Gaussian FWHM 0.30 eV, amplitude sum(f*Gaussian)/acid units; no Jacobian, no population weights. Plot 160-350 nm; all state ceilings exceed high-energy edge by >3 FWHM.",
              "nmr_method": "Equal-area Lorentzian per nucleus divided by acid units; H FWHM 0.08 ppm, C 0.80 ppm. One shared maximum per element across all systems for display. Static shieldings, no TMS, J splitting or rotational averaging.",
              "scope": "Two monomer conformers and three verified gas-phase cyclic-dimer local minima; no solution equilibrium prediction."}
    (destination / "summary.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print("Saved theoretical comparison for 2 monomers and 3 dimers; UV peaks (nm): "
          + ", ".join(f'{s["system"]}={s["broadened_peak_nm"]:.2f}' for s in summaries))
    return result


if __name__ == "__main__":
    build()
