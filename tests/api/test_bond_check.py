import numpy as np
import pytest
from helpers.utils import get_atom_num_and_pos_from_xyz

from irmsd.api.irmsd_exposed import check_bond_preservation, get_irmsd


def _parents_of_h(atom_numbers, positions):
    """Nearest heavy atom for every H atom."""
    heavy = np.where(atom_numbers > 1)[0]
    parents = {}
    for i in np.where(atom_numbers == 1)[0]:
        d = np.linalg.norm(positions[heavy] - positions[i], axis=1)
        parents[int(i)] = int(heavy[np.argmin(d)])
    return parents


def test_identity_mapping_is_bond_preserving(caffeine_irmsd_test_data):
    conformer1, _, _, _ = caffeine_irmsd_test_data
    Z, P = get_atom_num_and_pos_from_xyz(conformer1)

    nmis, per_atom = check_bond_preservation(Z, P, Z, P)

    assert nmis == 0
    assert not per_atom.any()


def test_cross_parent_h_swap_is_detected(caffeine_irmsd_test_data):
    conformer1, _, _, _ = caffeine_irmsd_test_data
    Z, P = get_atom_num_and_pos_from_xyz(conformer1)
    parents = _parents_of_h(Z, P)
    ha = next(iter(parents))
    hb = next(h for h, p in parents.items() if p != parents[ha])

    P2 = P.copy()
    P2[[ha, hb]] = P2[[hb, ha]]
    nmis, per_atom = check_bond_preservation(Z, P, Z, P2)

    # ha-ca and hb-cb are lost, ha-cb and hb-ca appear
    expected = np.zeros(len(Z), dtype=np.int32)
    expected[[ha, hb, parents[ha], parents[hb]]] = 2
    assert nmis == 4
    np.testing.assert_array_equal(per_atom, expected)


def test_irmsd_mapping_is_bond_preserving(caffeine_irmsd_test_data):
    """Caffeine has no same-rank atoms on different parents, so the iRMSD
    mapping of a shuffled conformer must preserve all bonds."""
    conformer1, conformer2, _, _ = caffeine_irmsd_test_data
    Z1, P1 = get_atom_num_and_pos_from_xyz(conformer1)
    Z2, P2 = get_atom_num_and_pos_from_xyz(conformer2)

    _, Z3, P3, Z4, P4 = get_irmsd(Z1, P1, Z2, P2, iinversion=1)
    nmis, _ = check_bond_preservation(Z3, P3, Z4, P4)

    assert nmis == 0


def test_shape_mismatch_raises():
    Z = np.array([6, 1, 1, 1, 1], dtype=np.int32)
    P = np.zeros((5, 3))
    with pytest.raises(ValueError):
        check_bond_preservation(Z, P, Z[:4], P[:4])
