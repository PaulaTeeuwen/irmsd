"""Regression tests for the bond-preserving assignment in get_irmsd
(centre walk plus terminal-atom refinement).

Fixture: RDKit conformers of 1-octanol, ibuprofen and decane. Each comment
line holds the exact bond-preserving RMSD to conformer 0 of the same
molecule (minimum over all graph automorphisms, proper rotations only).
"""

from pathlib import Path

import numpy as np
import pytest

from irmsd.api.irmsd_exposed import check_bond_preservation, get_irmsd

DATA = Path(__file__).resolve().parents[1] / "data" / "flexible_conformers.xyz"
SYMBOLS = {"H": 1, "C": 6, "O": 8}


def _read_blocks():
    lines = DATA.read_text().splitlines()
    blocks = {}
    i = 0
    while i < len(lines):
        nat = int(lines[i])
        name, _, k, exact = lines[i + 1].split()
        atoms = [line.split() for line in lines[i + 2 : i + 2 + nat]]
        Z = np.array([SYMBOLS[a[0]] for a in atoms], dtype=np.int32)
        P = np.array([[float(x) for x in a[1:4]] for a in atoms])
        blocks.setdefault(name, []).append((int(k), Z, P, float(exact.split("=")[1])))
        i += nat + 2
    return blocks


BLOCKS = _read_blocks()
PAIRS = [(name, k) for name, confs in BLOCKS.items() for k, *_ in confs if k > 0]


@pytest.mark.parametrize("name,k", PAIRS)
def test_mapping_is_bond_preserving_and_not_below_exact(name, k):
    confs = {c[0]: c for c in BLOCKS[name]}
    _, Z, P0, _ = confs[0]
    _, _, Pk, exact = confs[k]
    perm = np.random.default_rng(k).permutation(len(Z))

    r, Z3, P3, Z4, P4 = get_irmsd(Z, P0, Z[perm], Pk[perm], iinversion=2)
    nmis, _ = check_bond_preservation(Z3, P3, Z4, P4)

    assert nmis == 0
    assert r >= exact - 1e-5


@pytest.mark.parametrize("k", [k for name, k in PAIRS if name == "1-octanol"])
def test_octanol_reaches_exact_minimum(k):
    """All heavy atoms of 1-octanol are unique, so only the H assignment is
    open; the refinement must find the global minimum."""
    confs = {c[0]: c for c in BLOCKS["1-octanol"]}
    _, Z, P0, _ = confs[0]
    _, _, Pk, exact = confs[k]

    r, *_ = get_irmsd(Z, P0, Z, Pk, iinversion=2)

    assert r == pytest.approx(exact, abs=1e-4)


@pytest.mark.parametrize("name,k", PAIRS)
def test_close_to_exact(name, k):
    confs = {c[0]: c for c in BLOCKS[name]}
    _, Z, P0, _ = confs[0]
    _, _, Pk, exact = confs[k]
    perm = np.random.default_rng(k).permutation(len(Z))

    r, *_ = get_irmsd(Z, P0, Z[perm], Pk[perm], iinversion=2)

    assert r <= 1.01 * exact


@pytest.mark.parametrize(
    "name,k", [("decane", 1), ("decane", 3), ("ibuprofen", 4), ("ibuprofen", 7)]
)
def test_heavy_atom_symmetry_choice(name, k):
    """Pairs where the inertia/LSAP path picked the wrong heavy-atom
    symmetry class (decane chain direction, ibuprofen ring flip and
    isopropyl methyls); the centre walk must find the exact minimum."""
    confs = {c[0]: c for c in BLOCKS[name]}
    _, Z, P0, _ = confs[0]
    _, _, Pk, exact = confs[k]
    perm = np.random.default_rng(k).permutation(len(Z))

    r, *_ = get_irmsd(Z, P0, Z[perm], Pk[perm], iinversion=2)

    assert r == pytest.approx(exact, abs=1e-4)
