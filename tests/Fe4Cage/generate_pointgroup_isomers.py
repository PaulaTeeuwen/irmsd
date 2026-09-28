"""Regenerates pointgroup_isomers.json from Fe4Cage.xyz.

Self-contained: uses only irmsd's own get_symmetry_operations (an
independent, from-scratch geometric symmetry search, unrelated to
min_rmsd/get_irmsd's own alignment code) to find the real point group and
its operations directly from the structure -- no external tool or
reference data needed. Each non-identity operation IS a permutation that
maps the structure exactly onto itself (RMSD = 0.00 A by construction),
which is what makes these reliable ground truth for testing get_irmsd's
separate alignment pathway against.

    python generate_pointgroup_isomers.py
"""

import json
from pathlib import Path

from ase.io import read
from irmsd.api.symmetry_exposed import get_symmetry_operations

HERE = Path(__file__).parent

atoms = read(HERE / "Fe4Cage.xyz")
Z = atoms.get_atomic_numbers()
pos = atoms.get_positions()

symbol, ops = get_symmetry_operations(Z, pos, max_atoms=None)
print(f"point group: {symbol}, {len(ops)} operations found (incl. identity)")

isomers = [
    {"operation": op.label, "perm": [int(v) for v in op.permutation]}
    for op in ops
    if op.label != "E"
]

with open(HERE / "pointgroup_isomers.json", "w") as f:
    json.dump(isomers, f, indent=1)

print(f"wrote {len(isomers)} non-identity operations to pointgroup_isomers.json")
