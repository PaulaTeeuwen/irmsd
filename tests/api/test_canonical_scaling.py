"""Regression tests for the canonical rank setup on larger graphs.

The rank setup must be a function of the molecular graph only: permuting
the input atoms permutes the ranks, without renumbering them. The former
implementation broke this for long chains, where its floating-point
distance sums depended on the atom order.
"""

import numpy as np
import pytest

from irmsd.api.canonical_exposed import (
    get_canonical_fortran,
    get_canonical_from_connect_fortran,
)
from irmsd.api.irmsd_exposed import get_irmsd


def alkane(n):
    """All-trans n-alkane, carbons first, then the H atoms of each carbon."""
    dx, dy = 1.27, 0.87
    C = np.array([[k * dx, (k % 2) * dy, 0.0] for k in range(n)])
    H = []
    for k in range(n):
        out = -1.0 if k % 2 == 0 else 1.0  # away from the chain
        for z in (0.89, -0.89):
            H.append(C[k] + [0.0, out * 0.63, z])
        if k in (0, n - 1):
            H.append(C[k] + [(-1.0 if k == 0 else 1.0) * 1.03, out * 0.36, 0.0])
    Z = np.array([6] * n + [1] * len(H), dtype=np.int32)
    return Z, np.vstack([C, np.array(H)])


@pytest.mark.parametrize("n", [10, 45, 80])
def test_alkane_ranks_follow_atom_permutation(n):
    Z, P = alkane(n)
    r0 = get_canonical_fortran(Z, P)
    for seed in range(3):
        perm = np.random.default_rng(seed).permutation(len(Z))
        r = get_canonical_fortran(Z[perm], P[perm])
        np.testing.assert_array_equal(r, r0[perm])


@pytest.mark.parametrize("n", [10, 45, 80])
def test_alkane_carbon_classes(n):
    """C_k and C_(n+1-k) are equivalent, and nothing else among carbons."""
    Z, P = alkane(n)
    rc = get_canonical_fortran(Z, P)[:n]
    np.testing.assert_array_equal(rc, rc[::-1])
    assert len(set(rc)) == (n + 1) // 2


@pytest.mark.parametrize(
    "Z,P",
    [
        ([1, 1], [[0, 0, 0], [0.74, 0, 0]]),
        ([6, 1, 1, 1, 1], [[0, 0, 0], [0.63, 0.63, 0.63], [-0.63, -0.63, 0.63],
                           [-0.63, 0.63, -0.63], [0.63, -0.63, -0.63]]),
        ([8, 1, 1], [[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]]),
        ([18], [[0, 0, 0]]),
    ],
    ids=["H2", "methane", "water", "Ar"],
)
def test_tiny_systems(Z, P):
    """Graphs without any heavy-heavy bond used to crash the CANGEN path."""
    Z = np.array(Z, dtype=np.int32)
    P = np.array(P, dtype=float)
    r = get_canonical_fortran(Z, P)
    assert (r > 0).all()
    C = (np.linalg.norm(P[:, None] - P[None], axis=-1) < 1.2).astype(np.int32)
    np.fill_diagonal(C, 0)
    r = get_canonical_from_connect_fortran(Z, C)
    assert (r > 0).all()
    assert len(set(r[Z == 1])) <= 1


def test_prime_product_overflow():
    """Chain of 300 SF4 units: inner S atoms have six heavy neighbours, whose
    rank primes multiply beyond int64. Each unit and its mirror image must
    still form one class."""
    N = 300
    nat = 5 * N
    Z = np.array([16] * N + [9] * (4 * N), dtype=np.int32)
    C = np.zeros((nat, nat), dtype=np.int32)
    for i in range(N - 1):
        C[i, i + 1] = C[i + 1, i] = 1
    for i in range(N):
        for k in range(4):
            C[i, N + 4 * i + k] = C[N + 4 * i + k, i] = 1
    r = get_canonical_from_connect_fortran(Z, C, heavy=True)
    S, F = r[:N], r[N:].reshape(N, 4)
    np.testing.assert_array_equal(S, S[::-1])
    np.testing.assert_array_equal(F, F[::-1])
    assert (F == F[:, :1]).all()
    assert len(set(r)) == N


def test_rank_reuse_matches_recomputation():
    """Same atom order (ranks copied from the reference) and shuffled order
    (ranks recomputed) must give the same iRMSD."""
    Z, P = alkane(30)
    rng = np.random.default_rng(5)
    P2 = P + rng.normal(scale=0.05, size=P.shape)
    perm = rng.permutation(len(Z))
    r_same = get_irmsd(Z, P, Z, P2, iinversion=2)[0]
    r_perm = get_irmsd(Z, P, Z[perm], P2[perm], iinversion=2)[0]
    assert r_same == pytest.approx(r_perm, abs=1e-8)
