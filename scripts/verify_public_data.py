#!/usr/bin/env python3
"""Read-only scientific verification of the 16 published calculation bundles.

This checks calculation data, not the original cluster submission/archive history.
Run without Python's -O flag so that scientific assertions remain active.
"""
import itertools
import json
import re
import sys

import numpy as np

from gaussian_io import (BONDS, BOHR_TO_ANGSTROM, ROOT, coupling_matrix,
                         difference, file_hash, internal_hessian_frequencies,
                         last_orientation, parse_states, read_fchk, torsion)


def snapshot(system, kind):
    if kind == "opt":
        if system.startswith("conf"):
            return next(s for s in json.loads((ROOT / "analysis/stage1/audit.json").read_text(encoding="utf-8"))["systems"]
                        if s["id"] == system)
        path = ROOT / "analysis/dimers" / system / "audit.json"
    elif system.startswith("conf"):
        path = ROOT / "analysis" / ("conf01" if system == "conf-01" else "stage2") / "verified-results.json"
    else:
        path = ROOT / "analysis/dimer-properties" / system / kind / "verified-results.json"
    return json.loads(path.read_text(encoding="utf-8"))


def check_geometry(fields, text, system):
    units = 1 if system.startswith("conf") else 2
    numbers = ([6]*6 + [8]*2 + [1]*8)*units
    assert fields["Atomic numbers"] == numbers
    assert (fields["Number of atoms"], fields["Number of electrons"],
            fields["Charge"], fields["Multiplicity"]) == (16*units, 60*units, 0, 1)
    xyz = np.array(fields["Current cartesian coordinates"]).reshape(-1, 3)*BOHR_TO_ANGSTROM
    out_numbers, out_xyz = last_orientation(text)
    assert numbers == out_numbers and np.max(abs(xyz-out_xyz)) < 1e-6
    optimized = snapshot(system, "opt")
    assert np.max(abs(xyz-optimized["coordinates_angstrom"])) < 1e-7
    radii = {1: .31, 6: .76, 8: .66}
    bonds = set()
    for a, b in itertools.combinations(range(len(numbers)), 2):
        distance = np.linalg.norm(xyz[a]-xyz[b])
        assert distance > .7
        if distance < 1.2*(radii[numbers[a]]+radii[numbers[b]]):
            bonds.add((a+1, b+1))
    assert bonds == {(a+16*u, b+16*u) for u in range(units) for a,b,_ in BONDS}
    for offset in range(0, 16*units, 16):
        for indices in [(1,2,3,4), (3,4,5,6)]:
            angle = torsion(xyz.tolist(), [i+offset for i in indices])
            assert abs(difference(angle, 180)) < 5
    if units == 2:
        for donor, hydrogen, acceptor in [(7,15,22), (23,31,6)]:
            assert .9 < np.linalg.norm(xyz[donor]-xyz[hydrogen]) < 1.1
            assert 1.3 < np.linalg.norm(xyz[acceptor]-xyz[hydrogen]) < 2.1
    return xyz


def verify_job(job):
    system, kind = job["system"], job["kind"]
    saved = snapshot(system, kind)
    for ext, metadata in job["files"].items():
        path = ROOT / metadata["path"]
        assert path.stat().st_size == metadata["bytes"]
        assert file_hash(path) == metadata["published_sha256"], path
        original_name = metadata["original_path"].rsplit("/", 1)[1]
        if kind == "opt":
            expected = saved.get("sha256", saved.get("source_hashes"))[original_name]
        elif system.startswith("conf"):
            expected = saved["jobs"][kind][{"gjf":"input_sha256", "out":"output_sha256", "fchk":"fchk_sha256"}[ext]]
        else:
            expected = saved["source_hashes"][original_name]
        assert metadata["original_sha256"] == expected
        if ext != "out":
            assert metadata["original_sha256"] == metadata["published_sha256"]
    text = (ROOT / job["files"]["out"]["path"]).read_text(encoding="utf-8")
    fields = read_fchk(ROOT / job["files"]["fchk"]["path"])
    input_text = (ROOT / job["files"]["gjf"]["path"]).read_text(encoding="utf-8")
    assert "Error termination" not in text
    assert "Normal termination of Gaussian" in "\n".join(text.splitlines()[-10:])
    route = "".join(input_text.split()).upper()
    assert "%NPROCSHARED=4" in route and "%MEM=8GB" in route
    assert "PBE1PBE/6-31G(D)" in route
    xyz = check_geometry(fields, text, system)
    units = len(xyz)//16
    if kind == "opt":
        assert "Optimization completed" in text and "FREQ" in route
        assert "SCRF=" not in route
        frequencies = [float(v) for line in re.findall(r"Frequencies --([^\n]+)", text) for v in line.split()]
        count = 3*len(xyz)-6
        assert len(frequencies) == count and min(frequencies) > 0
        assert frequencies == saved["frequencies_cm1"]
        assert np.max(abs(np.array(fields["Vib-E2"][:count])-frequencies)) < 6e-5
        derived = internal_hessian_frequencies(fields, xyz)
        assert min(derived) > 0 and np.max(abs(derived-frequencies)) < .05
        assert abs(fields["SCF Energy"]-saved["E_hartree"]) < 1e-9
        for label, key in [("zero-point Energies", "E_ZPE"), ("thermal Enthalpies", "H"), ("thermal Free Energies", "G")]:
            value = float(re.findall(r"Sum of electronic and "+label+r"=\s*(-?\d+\.\d+)", text)[-1])
            assert value == saved["thermochemistry_hartree"][key]
    elif kind == "uv":
        count = 20*units
        assert f"TD=(SINGLETS,NSTATES={count})" in route
        assert "SCRF=(SMD,SOLVENT=DIETHYLETHER)" in route
        assert "Solvent              : DiethylEther" in text
        states = parse_states(text, count)
        assert states == saved["uv_states"]
        assert fields["ETran scalars"][:2] == [count, 16]
        assert fields["ETran spin"] == [0]*count
        response = np.array(fields["ETran state values"]).reshape(count, 16)
        excitation = response[:,0]-fields["SCF Energy"]
        strengths = 2/3*excitation*np.sum(response[:,1:4]**2, axis=1)
        assert np.max(abs(excitation*27.211386245988-[s["energy_ev"] for s in states])) < .0001
        assert np.max(abs(strengths-[s["oscillator_strength"] for s in states])) < .0001
        for orbital in saved.get("frontier_orbitals", []):
            assert fields["Alpha Orbital Energies"][orbital["mo"]-1] == orbital["energy_hartree"]
    else:
        assert kind == "nmr" and "NMR=(GIAO,SPINSPIN,MIXED)" in route and "SCRF" not in route
        assert abs(fields["SCF Energy"]-snapshot(system,"opt")["E_hartree"]) < 1e-8
        shield = [{"atom":int(n), "element":element, "shielding_ppm":float(value), "anisotropy_ppm":float(anisotropy)}
                  for n,element,value,anisotropy in re.findall(
                      r"(\d+)\s+([CHO])\s+Isotropic =\s*([-\d.]+)\s+Anisotropy =\s*([-\d.]+)", text)]
        assert len(shield) == len(xyz) and shield == saved["nmr_shieldings"]
        tensor = np.array(fields["NMR shielding"]).reshape(len(xyz),3,3)*7.2973525643e-3**2*1e6
        assert np.max(abs(tensor.trace(axis1=1,axis2=2)/3-[s["shielding_ppm"] for s in shield])) < .0005
        total = coupling_matrix(text, "Total nuclear spin-spin coupling J (Hz):", len(xyz))
        assert np.array_equal(total, saved["nmr_j_hz"])
        parts = []
        for key, heading in {"FC":"Fermi Contact (FC)", "SD":"Spin-dipolar (SD)",
                             "PSO":"Paramagnetic spin-orbit (PSO)", "DSO":"Diamagnetic spin-orbit (DSO)"}.items():
            part = coupling_matrix(text, heading+" contribution to J (Hz):", len(xyz))
            assert np.array_equal(part, saved["nmr_j_components_hz"][key])
            parts.append(part)
        assert np.max(abs(sum(parts)-total)) < .001
    print(f"PASS {system}/{kind}")


def main():
    if sys.flags.optimize:
        raise SystemExit("Run without -O: scientific assertions must remain enabled.")
    manifest = json.loads((ROOT / "data/manifest.json").read_text(encoding="utf-8"))
    expected = {(f"conf-{i:02}", "opt") for i in [1,2,3]}
    expected |= {(f"conf-{i:02}", kind) for i in [1,2] for kind in ["uv","nmr"]}
    expected |= {(f"dimer-{pair}", kind) for pair in ["01-01","01-02","02-02"] for kind in ["opt","uv","nmr"]}
    assert len(manifest["jobs"]) == 16
    assert {(j["system"],j["kind"]) for j in manifest["jobs"]} == expected
    for job in manifest["jobs"]:
        verify_job(job)
    print("PASS: all 16 public calculation bundles match the scientific snapshots.")


if __name__ == "__main__":
    main()
