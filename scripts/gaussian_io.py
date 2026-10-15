"""Shared scientific readers and geometry helpers for offline reproduction."""
import csv
import hashlib
import math
import re
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
BOHR_TO_ANGSTROM = 0.529177210903
BONDS = [(1,2,1),(2,3,2),(3,4,1),(4,5,2),(5,6,1),(6,7,2),(6,8,1),
         (1,9,1),(1,10,1),(1,11,1),(2,12,1),(3,13,1),(4,14,1),(5,15,1),(8,16,1)]

def read_fchk(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    fields = {}
    for i, line in enumerate(lines):
        match = re.match(r"^(.{40,43}?)\s+([IR])\s+(?:N=\s*(\d+)|([-+\d.EeDd]+))\s*$", line)
        if not match:
            continue
        name, kind, count, value = match.groups()
        convert = int if kind == "I" else lambda s: float(s.replace("D", "E"))
        if count:
            values, j = [], i + 1
            while len(values) < int(count):
                values.extend(lines[j].split())
                j += 1
            assert len(values) == int(count), name
            fields[name.strip()] = [convert(s) for s in values]
        else:
            fields[name.strip()] = convert(value)
    return fields


def last_orientation(text):
    pattern = r"(?:Input|Standard) orientation:\s*\n\s*-+\n.*?\n.*?\n\s*-+\n(.*?)\n\s*-+"
    rows = re.findall(pattern, text, re.S)[-1].strip().splitlines()
    atoms = [int(row.split()[1]) for row in rows]
    xyz = np.array([[float(x) for x in row.split()[-3:]] for row in rows])
    return atoms, xyz


def fitted_rmsd(first, second):
    first = first - first.mean(axis=0)
    second = second - second.mean(axis=0)
    u, _, vt = np.linalg.svd(first.T @ second)
    correction = np.diag([1, 1, np.linalg.det(u @ vt)])
    rotation = u @ correction @ vt
    return float(np.sqrt(np.mean(np.sum((first @ rotation - second) ** 2, axis=1))))


def internal_hessian_frequencies(fields, xyz):
    # Project out the six mass-weighted rigid translations/rotations, then
    # independently diagonalize the fchk Cartesian Hessian (Hartree/bohr^2).
    masses = np.array(fields["Real atomic weights"])
    size = len(masses) * 3
    hessian = np.zeros((size, size))
    hessian[np.tril_indices(size)] = fields["Cartesian Force Constants"]
    hessian = hessian + hessian.T - np.diag(hessian.diagonal())
    weighted = hessian / np.sqrt(np.repeat(masses, 3)[:, None] * np.repeat(masses, 3)[None, :])
    centered = xyz - np.average(xyz, axis=0, weights=masses)
    external = []
    for axis in np.eye(3):
        external.append((np.sqrt(masses)[:, None] * np.tile(axis, (len(masses), 1))).ravel())
        external.append((np.sqrt(masses)[:, None] * np.cross(axis, centered)).ravel())
    u, singular, _ = np.linalg.svd(np.array(external).T, full_matrices=True)
    assert np.sum(singular > 1e-8) == 6
    internal = u[:, 6:]
    values = np.linalg.eigvalsh(internal.T @ weighted @ internal)
    factor = math.sqrt(4.3597447222071e-18 / (1.66053906660e-27 * (BOHR_TO_ANGSTROM * 1e-10) ** 2)) / (2 * math.pi * 299792458 * 100)
    return np.sign(values) * np.sqrt(np.abs(values)) * factor


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def coupling_matrix(text, heading, size=16):
    block = text.rsplit(heading, 1)[1]
    matrix = np.full((size, size), np.nan)
    columns = []
    for line in block.splitlines():
        words = line.split()
        if not words:
            continue
        if all(w.isdigit() for w in words):
            columns = [int(w)-1 for w in words]
            assert columns == sorted(set(columns)) and all(0 <= col < size for col in columns)
        elif words[0].isdigit() and all(re.fullmatch(r"[-+]?\d+\.\d+D[-+]\d+", w) for w in words[1:]):
            row = int(words[0])-1
            assert 0 <= row < size
            expected_columns = [col for col in columns if col <= row]
            assert len(words) - 1 == len(expected_columns)
            for col, value in zip(columns, words[1:]):
                assert np.isnan(matrix[row, col]), (heading, row, col)
                matrix[row, col] = matrix[col, row] = float(value.replace("D", "E"))
            if np.isfinite(matrix).all():
                break
        elif columns:
            break
    assert np.isfinite(matrix).all(), heading
    return matrix


def save_csv(destination, name, rows):
    with (destination / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_states(text, expected=40):
    pattern = (r"Excited State\s+(\d+):\s+Singlet-\?Sym\s+([\d.]+) eV\s+"
               r"([\d.]+) nm\s+f=([\d.]+).*?\n"
               r"((?:\s*\d+\s*->\s*\d+\s+[-\d.]+\s*\n)+)")
    states = []
    for match in re.finditer(pattern, text):
        n, ev, nm, strength, coefficients = match.groups()
        states.append({"state": int(n), "energy_ev": float(ev),
                       "wavelength_nm": float(nm), "oscillator_strength": float(strength),
                       "transitions": [{"occupied": int(a), "virtual": int(b),
                                        "printed_amplitude": float(c)}
                                       for a, b, c in re.findall(
                                           r"(\d+)\s*->\s*(\d+)\s+([-\d.]+)", coefficients)]})
    assert [s["state"] for s in states] == list(range(1, expected + 1))
    assert np.all(np.diff([s["energy_ev"] for s in states]) >= 0)
    assert all(s["energy_ev"] > 0 and s["oscillator_strength"] >= 0 for s in states)
    assert max(abs(1239.841984 / s["energy_ev"] - s["wavelength_nm"]) for s in states) < 0.01
    return states


def sub(a, b):
    return [x-y for x, y in zip(a, b)]


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def unit(a):
    norm = math.sqrt(dot(a, a))
    return [x/norm for x in a]


def torsion(xyz, indices):
    a, b, c, d = [xyz[i-1] for i in indices]
    axis = unit(sub(c, b))
    left, right = sub(a, b), sub(d, c)
    left = [x-dot(left, axis)*y for x, y in zip(left, axis)]
    right = [x-dot(right, axis)*y for x, y in zip(right, axis)]
    return math.degrees(math.atan2(dot(cross(axis, left), right), dot(left, right)))


def difference(a, b):
    return (a-b+180) % 360 - 180
