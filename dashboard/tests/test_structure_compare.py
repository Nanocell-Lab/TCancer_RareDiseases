from __future__ import annotations

import unittest

from rare_oncology_dashboard.structure_compare import (
    PDBComparisonError,
    compare_pdb_models,
    parse_pdb,
)


def pdb(points: list[tuple[str, int, float, float, float]]) -> bytes:
    lines = []
    for serial, (residue, number, x, y, z) in enumerate(points, 1):
        lines.append(
            f"ATOM  {serial:5d}  CA  {residue:>3s} A{number:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 90.00           C"
        )
    return ("\n".join(lines) + "\nEND\n").encode("ascii")


class PDBComparisonTests(unittest.TestCase):
    def test_kabsch_removes_rigid_rotation_and_translation(self) -> None:
        wildtype = pdb(
            [("ALA", 1, 0, 0, 0), ("GLY", 2, 1, 0, 0), ("SER", 3, 0, 1, 0)]
        )
        variant = pdb(
            [("ALA", 1, 5, 2, 1), ("GLY", 2, 5, 3, 1), ("VAL", 3, 4, 2, 1)]
        )
        result = compare_pdb_models(wildtype, variant)
        self.assertAlmostEqual(result["global_ca_rmsd_angstrom"], 0.0, places=6)
        self.assertEqual(result["changed_residue_count"], 1)
        self.assertEqual(result["common_ca_count"], 3)

    def test_parser_rejects_missing_ca_and_comparison_requires_common_ids(self) -> None:
        with self.assertRaisesRegex(PDBComparisonError, "no readable C-alpha"):
            parse_pdb(b"HEADER empty\n", "empty.pdb")
        wildtype = pdb([("ALA", 1, 0, 0, 0), ("GLY", 2, 1, 0, 0), ("SER", 3, 0, 1, 0)])
        variant = pdb([("ALA", 10, 0, 0, 0), ("GLY", 11, 1, 0, 0), ("SER", 12, 0, 1, 0)])
        with self.assertRaisesRegex(PDBComparisonError, "At least three"):
            compare_pdb_models(wildtype, variant)


if __name__ == "__main__":
    unittest.main()
