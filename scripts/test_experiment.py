"""Offline checks for experimental reference transforms; no external requests."""
import csv
import json
import unittest

import numpy as np

from compare_experiment import ROOT, fwhm, read_uv_reference


class ExperimentalComparisonTests(unittest.TestCase):
    def test_raw_uv_identity_and_linearization(self):
        x, y = read_uv_reference()
        self.assertEqual(len(x), 126)
        self.assertEqual(list(x[y == max(y)]), [255.56, 256.41])
        rows = list(csv.DictReader((ROOT / "analysis/experimental/uv-reference.csv").read_text(encoding="utf-8").splitlines()))
        for row, xx, yy in zip(rows, x, y):
            self.assertEqual(float(row["wavelength_nm"]), xx)
            self.assertAlmostEqual(float(row["normalized_intensity"]), 10**(yy-4.443), places=12)

    def test_half_height_interpolation(self):
        x = np.linspace(-3, 3, 12001)
        y = np.exp(-4*np.log(2)*x*x/1.2**2)
        self.assertAlmostEqual(fwhm(x,y)["width_nm"],1.2,places=6)
        self.assertIsNone(fwhm(np.array([0,1,2]),np.array([.7,1,.7])))

    def test_proton_transcription_and_atom_mapping(self):
        ref = json.loads((ROOT / "references/experimental/nmr-h-reference.json").read_text(encoding="utf-8"))
        self.assertEqual([p["delta_ppm"] for p in ref["peaks"]], [1.84,6.26,6.26,7.18,5.79,12.2])
        self.assertEqual([a for p in ref["peaks"] for a in p["atoms"]], list(range(9,17)))
        self.assertEqual(sum(p["count"] for p in ref["peaks"]),8)
        self.assertEqual(ref["reference"],"DSS-d6 internal standard, delta=0 ppm")

    def test_relative_shift_sign_and_one_species_anchor(self):
        raw = list(csv.DictReader((ROOT / "analysis/comparison/nmr-shieldings.csv").read_text(encoding="utf-8").splitlines()))
        calc = list(csv.DictReader((ROOT / "analysis/experimental/nmr-H-relative-shifts.csv").read_text(encoding="utf-8").splitlines()))
        self.assertEqual(len(calc),54)
        for r in calc:
            if r["system"] == "reference":
                continue
            s = {int(p["atom"]):float(p["shielding_ppm"]) for p in raw if p["system"]==r["system"]}
            units = 1 if len(s)==16 else 2
            anchor = sum(s[a+16*f] for f in range(units) for a in [9,10,11])/(3*units)
            f = int(r["fragment"])-1
            atoms = [9,10,11] if r["site"]=="CH3" else [16 if r["site"]=="OH" else int(r["site"][1:])]
            sigma = sum(s[a+16*f] for a in atoms)/len(atoms)
            self.assertAlmostEqual(float(r["relative_shift_ppm"]),anchor-sigma,places=10)
            self.assertAlmostEqual(float(r["error_ppm"]),anchor-sigma-float(r["reference_relative_shift_ppm"]),places=10)

    def test_unshifted_uv_peaks_and_metrics(self):
        summary = json.loads((ROOT / "analysis/experimental/summary.json").read_text(encoding="utf-8"))
        theory = json.loads((ROOT / "analysis/comparison/summary.json").read_text(encoding="utf-8"))
        for a,b in zip(summary["uv"]["calculations"],theory["systems"]):
            self.assertAlmostEqual(a["peak_nm"],b["broadened_peak_nm"],places=2)
        self.assertIsNone(summary["uv"]["solvent"])
        for row in summary["nmr_H"]["metrics"]:
            self.assertGreater(row["alkene_RMSE_ppm"],row["alkene_MAE_ppm"])
            self.assertLess(row["alkene_MAE_ppm"],.16)

    def test_carbon_transcription_and_identity(self):
        ref = json.loads((ROOT / "references/experimental/nmr-c-reference.json").read_text(encoding="utf-8"))
        self.assertEqual(ref["cas"],"110-44-1")
        self.assertEqual(ref["solvent"],"CDCl3")
        self.assertEqual(ref["reference"],"TMS internal standard")
        self.assertEqual([p["source_carbon_number"] for p in ref["peaks"]],list(range(1,7)))
        self.assertEqual([p["delta_ppm"] for p in ref["peaks"]],[173.00,118.25,147.34,140.74,129.73,18.70])
        self.assertEqual([p["group"] for p in ref["peaks"]],["COOH",*(4*["alkene"]),"CH3"])
        self.assertIsNone(ref["temperature_C"])
        self.assertEqual(ref["doi"], "10.1248/jhs1956.37.89")
        self.assertTrue(ref["url"].startswith("https://www.jstage.jst.go.jp/"))

    def test_carbon_relative_shifts_and_fragment_anchor(self):
        raw = list(csv.DictReader((ROOT / "analysis/comparison/nmr-shieldings.csv").read_text(encoding="utf-8").splitlines()))
        rows = list(csv.DictReader((ROOT / "analysis/experimental/nmr-C-relative-shifts.csv").read_text(encoding="utf-8").splitlines()))
        self.assertEqual(len(rows),54)
        ref = {r["site"]:float(r["relative_shift_ppm"]) for r in rows if r["system"]=="reference"}
        for name in dict.fromkeys(r["system"] for r in raw):
            s = {int(r["atom"]):float(r["shielding_ppm"]) for r in raw if r["system"]==name and r["element"]=="C"}
            units = len(s)//6
            anchor = sum(s[1+16*f] for f in range(units))/units
            for f in range(units):
                selected = [r for r in rows if r["system"]==name and int(r["fragment"])==f+1]
                self.assertEqual(sorted(int(r["atom"]) for r in selected),list(range(1+16*f,7+16*f)))
                alkene = [float(r["relative_shift_ppm"]) for r in selected if r["site"].startswith("alkene-")]
                self.assertEqual(alkene,sorted(alkene))
                for r in selected:
                    position = anchor-s[int(r["atom"])]
                    self.assertAlmostEqual(float(r["relative_shift_ppm"]),position,places=10)
                    self.assertAlmostEqual(float(r["error_ppm"]),position-ref[r["site"]],places=10)
            # Independent MAE from sorted raw shielding-derived spacings.
            expected = []
            for f in range(units):
                calc = sorted(anchor-s[a+16*f] for a in range(2,6))
                target = [99.55,111.03,122.04,128.64]
                expected.extend(abs(c-t) for c,t in zip(calc,target))
            summary = json.loads((ROOT / "analysis/experimental/summary.json").read_text(encoding="utf-8"))
            metric = next(m for m in summary["nmr_C"]["metrics"] if m["system"]==name)
            self.assertAlmostEqual(metric["alkene_peak_set_MAE_ppm"],sum(expected)/len(expected),places=10)
            self.assertEqual(metric["alkene_peak_count"],4*units)

    def test_carbonyl_comparison_keeps_nonzero_model_error(self):
        summary = json.loads((ROOT / "analysis/experimental/summary.json").read_text(encoding="utf-8"))
        metrics = summary["nmr_C"]["metrics"]
        self.assertAlmostEqual(metrics[0]["carbonyl_relative_shifts_ppm"][0],135.0296,places=6)
        self.assertAlmostEqual(metrics[-1]["carbonyl_residuals_ppm"][0],-10.9466,places=6)
        # A methyl anchor is not a fit that forces every carbon onto experiment.
        self.assertTrue(all(v < -10 for m in metrics for v in m["carbonyl_residuals_ppm"]))

    def test_assigned_carbons_keep_source_numbering(self):
        rows = list(csv.DictReader((ROOT / "analysis/experimental/nmr-C-assigned-relative-shifts.csv").read_text(encoding="utf-8").splitlines()))
        raw = list(csv.DictReader((ROOT / "analysis/comparison/nmr-shieldings.csv").read_text(encoding="utf-8").splitlines()))
        self.assertEqual(len(rows),54)
        expected = {1:0.0,2:111.03,3:122.04,4:128.64,5:99.55,6:154.30}
        summary = json.loads((ROOT / "analysis/experimental/summary.json").read_text(encoding="utf-8"))
        for name in dict.fromkeys(r["system"] for r in raw):
            s = {int(r["atom"]):float(r["shielding_ppm"]) for r in raw if r["system"]==name and r["element"]=="C"}
            units = len(s)//6
            anchor = sum(s[1+16*f] for f in range(units))/units
            errors = []
            selected = [r for r in rows if r["system"]==name]
            self.assertEqual(len(selected),6*units)
            for r in selected:
                atom = int(r["atom"])
                site = (atom-1)%16+1
                self.assertEqual(int(r["source_carbon_number"]),7-site)
                self.assertAlmostEqual(float(r["reference_relative_shift_ppm"]),expected[site],places=10)
                self.assertAlmostEqual(float(r["relative_shift_ppm"]),anchor-s[atom],places=10)
                if site in range(2,6):
                    errors.append(abs(anchor-s[atom]-expected[site]))
            metric = next(m for m in summary["nmr_C"]["metrics"] if m["system"]==name)
            self.assertAlmostEqual(metric["alkene_assigned_MAE_ppm"],sum(errors)/len(errors),places=10)
        conf02 = next(m for m in summary["nmr_C"]["metrics"] if m["system"]=="conf-02")
        self.assertAlmostEqual(conf02["alkene_assigned_MAE_ppm"],9.081275,places=6)

    def test_experimental_coupling_mapping(self):
        summary = json.loads((ROOT / "analysis/experimental/summary.json").read_text(encoding="utf-8"))
        rows = summary["nmr_J"]["comparisons"]
        self.assertEqual([r["pair"] for r in rows],["H14--H15","H13--H14"])
        self.assertEqual([r["experimental_magnitude_hz"] for r in rows],[15.2,10.0])
        self.assertAlmostEqual(rows[0]["difference_hz"],.4217,places=6)
        self.assertAlmostEqual(rows[1]["difference_hz"],1.5586,places=6)


if __name__ == "__main__":
    unittest.main()
