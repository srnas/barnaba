"""Inter-residue torsions must not bridge a residue skipped by Nucleic."""
from __future__ import absolute_import, division, print_function

import os
import tempfile

import barnaba as bb
import mdtraj as md
import numpy as np


DATA = os.path.join(os.path.dirname(__file__), "data", "UUCG.pdb")
EPSILON_COUPLINGS = ["C4Pe", "H3P"]


def _trajectory():
    return md.load(DATA)


def _skip_unknown_residue():
    traj = _trajectory()
    traj.topology.residue(3).name = "UNK"
    return traj


def _skip_incomplete_residue():
    traj = _trajectory()
    missing = traj.topology.residue(3).atom("C2").index
    return traj.atom_slice([atom.index for atom in traj.topology.atoms
                            if atom.index != missing])


def _check_skipped_backbone(traj):
    original, names = bb.backbone_angles_traj(_trajectory())
    actual, recognized = bb.backbone_angles_traj(traj)
    kept = [i for i in range(len(names)) if i != 3]
    assert recognized == [names[i] for i in kept]
    expected = original[:, kept, :].copy()
    # U3 cannot obtain epsilon/zeta from C5 across the skipped U4.
    expected[:, 2, 4:6] = np.nan
    # C5 cannot obtain alpha from U3 across the skipped U4.
    expected[:, 3, 0] = np.nan
    np.testing.assert_array_equal(actual, expected)


def _check_skipped_couplings(traj):
    for raw in [False, True]:
        original, names = bb.jcouplings_traj(
            _trajectory(), couplings=EPSILON_COUPLINGS, raw=raw)
        actual, recognized = bb.jcouplings_traj(
            traj, couplings=EPSILON_COUPLINGS, raw=raw)
        kept = [i for i in range(len(names)) if i != 3]
        assert recognized == [names[i] for i in kept]
        expected = original[:, kept, :].copy()
        expected[:, 2, :] = np.nan
        # Every unaffected value and pre-existing NaN must be unchanged.
        np.testing.assert_array_equal(actual, expected)


def test_unknown_residue_does_not_bridge_backbone_torsions():
    _check_skipped_backbone(_skip_unknown_residue())


def test_incomplete_residue_does_not_bridge_backbone_torsions():
    _check_skipped_backbone(_skip_incomplete_residue())


def test_unknown_residue_does_not_bridge_epsilon_couplings():
    _check_skipped_couplings(_skip_unknown_residue())


def test_incomplete_residue_does_not_bridge_epsilon_couplings():
    _check_skipped_couplings(_skip_incomplete_residue())


def test_skipped_residue_subset_preserves_requested_order_and_missing_angles():
    traj = _skip_unknown_residue()
    full, recognized = bb.backbone_angles_traj(traj)
    requested = [recognized[3], recognized[2]]
    selected, names = bb.backbone_angles_traj(
        traj, residues=requested, angles=["zeta", "alpha", "eps"])
    assert names == requested
    np.testing.assert_array_equal(selected, full[:, [3, 2], :][:, :, [5, 0, 4]])
    assert np.isnan(selected[0, 0, 1])
    assert np.isnan(selected[0, 1, [0, 2]]).all()
    couplings, names = bb.jcouplings_traj(
        traj, residues=requested, couplings=EPSILON_COUPLINGS)
    assert names == requested
    assert np.isfinite(couplings[0, 0]).all()
    assert np.isnan(couplings[0, 1]).all()


def test_contiguous_fixture_has_only_existing_terminal_missing_angles():
    traj = _trajectory()
    angles, names = bb.backbone_angles_traj(
        traj, angles=["alpha", "eps", "zeta"])
    assert len(names) == 8
    expected_missing = np.zeros((1, 8, 3), dtype=bool)
    expected_missing[:, 0, 0] = True
    expected_missing[:, -1, 1:] = True
    np.testing.assert_array_equal(np.isnan(angles), expected_missing)
    couplings, _ = bb.jcouplings_traj(traj, couplings=EPSILON_COUPLINGS)
    assert np.isfinite(couplings[:, :-1, :]).all()
    assert np.isnan(couplings[:, -1, :]).all()


def test_residue_numbers_and_insertion_codes_do_not_break_topology_neighbors():
    original = _trajectory()
    renamed = md.Trajectory(original.xyz.copy(), original.topology.copy())
    for residue, number in zip(renamed.topology.residues,
                               [10, 20, 40, 41, 75, 100, 150, 200]):
        residue.resSeq = number
    expected, _ = bb.backbone_angles_traj(original)
    actual, _ = bb.backbone_angles_traj(renamed)
    np.testing.assert_array_equal(actual, expected)
    expected, _ = bb.jcouplings_traj(original, couplings=EPSILON_COUPLINGS)
    actual, _ = bb.jcouplings_traj(renamed, couplings=EPSILON_COUPLINGS)
    np.testing.assert_array_equal(actual, expected)
    with open(DATA) as source:
        lines = source.readlines()
    renamed = []
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            index = int(line[22:26]) - 1
            number = 10 + index // 2
            insertion = "A" if index % 2 else " "
            line = line[:22] + "%4d%s" % (number, insertion) + line[27:]
        renamed.append(line)
    with tempfile.TemporaryDirectory() as directory:
        path = os.path.join(directory, "insertions.pdb")
        with open(path, "w") as target:
            target.writelines(renamed)
        inserted = md.load(path)
    original = _trajectory()
    assert inserted.n_residues == original.n_residues
    np.testing.assert_array_equal(inserted.xyz, original.xyz)
    expected, _ = bb.backbone_angles_traj(original)
    actual, _ = bb.backbone_angles_traj(inserted)
    np.testing.assert_array_equal(actual, expected)
    expected, _ = bb.jcouplings_traj(original, couplings=EPSILON_COUPLINGS)
    actual, _ = bb.jcouplings_traj(inserted, couplings=EPSILON_COUPLINGS)
    np.testing.assert_array_equal(actual, expected)


def test_chain_boundaries_still_leave_inter_residue_angles_missing():
    original = _trajectory()
    topology = md.Topology()
    chains = [topology.add_chain(), topology.add_chain()]
    atoms = {}
    for residue in original.topology.residues:
        target = topology.add_residue(residue.name, chains[int(residue.index >= 4)],
                                      resSeq=residue.resSeq)
        for atom in residue.atoms:
            atoms[atom.index] = topology.add_atom(atom.name, atom.element, target)
    for a, b in original.topology.bonds:
        if atoms[a.index].residue.chain == atoms[b.index].residue.chain:
            topology.add_bond(atoms[a.index], atoms[b.index])
    split = md.Trajectory(original.xyz.copy(), topology)
    expected, _ = bb.backbone_angles_traj(original)
    actual, names = bb.backbone_angles_traj(split)
    assert names[3].endswith("_0") and names[4].endswith("_1")
    expected[:, 3, 4:6] = np.nan
    expected[:, 4, 0] = np.nan
    np.testing.assert_array_equal(actual, expected)
    expected, _ = bb.jcouplings_traj(original, couplings=EPSILON_COUPLINGS)
    actual, _ = bb.jcouplings_traj(split, couplings=EPSILON_COUPLINGS)
    expected[:, 3, :] = np.nan
    np.testing.assert_array_equal(actual, expected)
