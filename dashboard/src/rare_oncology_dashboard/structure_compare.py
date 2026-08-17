"""Safe C-alpha comparison for uploaded wild-type and variant PDB models."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np


MAX_PDB_BYTES = 25 * 1024 * 1024
MAX_CA_ATOMS = 20_000


class PDBComparisonError(ValueError):
    """Raised when uploaded structures cannot be compared reliably."""


@dataclass(frozen=True)
class CAlphaAtom:
    chain: str
    residue_number: int
    insertion_code: str
    residue_name: str
    x: float
    y: float
    z: float

    @property
    def key(self) -> tuple[str, int, str]:
        """Role: Return the stable residue key used to align protein structures."""
        return (self.chain, self.residue_number, self.insertion_code)


@dataclass(frozen=True)
class ParsedPDB:
    name: str
    atoms: tuple[CAlphaAtom, ...]


def compare_pdb_models(
    wildtype_data: bytes,
    variant_data: bytes,
    *,
    wildtype_name: str = "wildtype.pdb",
    variant_name: str = "variant.pdb",
) -> dict[str, Any]:
    """Superpose common C-alpha atoms with Kabsch and return residue shifts."""

    wildtype = parse_pdb(wildtype_data, wildtype_name)
    variant = parse_pdb(variant_data, variant_name)
    wildtype_by_key = {atom.key: atom for atom in wildtype.atoms}
    variant_by_key = {atom.key: atom for atom in variant.atoms}
    common_keys = sorted(wildtype_by_key.keys() & variant_by_key.keys())
    if len(common_keys) < 3:
        raise PDBComparisonError(
            "At least three C-alpha residues with matching chain and residue IDs are required"
        )
    reference = np.asarray(
        [
            (wildtype_by_key[key].x, wildtype_by_key[key].y, wildtype_by_key[key].z)
            for key in common_keys
        ],
        dtype=float,
    )
    mobile = np.asarray(
        [
            (variant_by_key[key].x, variant_by_key[key].y, variant_by_key[key].z)
            for key in common_keys
        ],
        dtype=float,
    )
    aligned = _kabsch_align(mobile, reference)
    displacements = np.linalg.norm(aligned - reference, axis=1)
    rmsd = float(np.sqrt(np.mean(np.square(displacements))))
    rows: list[dict[str, Any]] = []
    for index, key in enumerate(common_keys):
        reference_atom = wildtype_by_key[key]
        variant_atom = variant_by_key[key]
        rows.append(
            {
                "chain": key[0] or "_",
                "residue_number": key[1],
                "insertion_code": key[2],
                "wildtype_residue": reference_atom.residue_name,
                "variant_residue": variant_atom.residue_name,
                "residue_changed": (
                    reference_atom.residue_name != variant_atom.residue_name
                ),
                "ca_displacement_angstrom": float(displacements[index]),
            }
        )
    coverage = len(common_keys) / max(len(wildtype.atoms), len(variant.atoms))
    return {
        "wildtype_name": wildtype.name,
        "variant_name": variant.name,
        "wildtype_ca_count": len(wildtype.atoms),
        "variant_ca_count": len(variant.atoms),
        "common_ca_count": len(common_keys),
        "coverage_fraction": coverage,
        "global_ca_rmsd_angstrom": rmsd,
        "mean_ca_displacement_angstrom": float(np.mean(displacements)),
        "max_ca_displacement_angstrom": float(np.max(displacements)),
        "changed_residue_count": sum(bool(row["residue_changed"]) for row in rows),
        "residues": rows,
        "wildtype_coordinates": reference.tolist(),
        "aligned_variant_coordinates": aligned.tolist(),
        "caveat": (
            "C-alpha comparison after a global Kabsch superposition. Matching uses "
            "chain, residue number and insertion code; different numbering, missing "
            "segments, multimers and model confidence can change the interpretation."
        ),
    }


def parse_pdb(data: bytes, name: str) -> ParsedPDB:
    """Role: Parse and validate pdb from its source representation."""
    if not data:
        raise PDBComparisonError(f"{name} is empty")
    if len(data) > MAX_PDB_BYTES:
        raise PDBComparisonError(f"{name} exceeds the 25 MB upload limit")
    text = data.decode("utf-8", errors="replace")
    atoms: dict[tuple[str, int, str], CAlphaAtom] = {}
    saw_model = False
    for line in text.splitlines():
        record = line[:6].strip()
        if record == "MODEL":
            if saw_model:
                break
            saw_model = True
            continue
        if record == "ENDMDL" and atoms:
            break
        if record != "ATOM" or len(line) < 54:
            continue
        if line[12:16].strip() != "CA" or line[16:17] not in {" ", "A"}:
            continue
        try:
            residue_number = int(line[22:26].strip())
            coordinates = tuple(
                float(line[start:end]) for start, end in ((30, 38), (38, 46), (46, 54))
            )
        except ValueError:
            continue
        if not all(math.isfinite(value) for value in coordinates):
            continue
        atom = CAlphaAtom(
            chain=line[21:22].strip(),
            residue_number=residue_number,
            insertion_code=line[26:27].strip(),
            residue_name=line[17:20].strip().upper() or "UNK",
            x=coordinates[0],
            y=coordinates[1],
            z=coordinates[2],
        )
        atoms.setdefault(atom.key, atom)
        if len(atoms) > MAX_CA_ATOMS:
            raise PDBComparisonError(f"{name} contains too many C-alpha atoms")
    if not atoms:
        raise PDBComparisonError(f"{name} contains no readable C-alpha ATOM records")
    return ParsedPDB(name=name, atoms=tuple(atoms.values()))


def _kabsch_align(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Role: Compute kabsch align for the structure_compare workflow."""
    mobile_center = np.mean(mobile, axis=0)
    reference_center = np.mean(reference, axis=0)
    centered_mobile = mobile - mobile_center
    centered_reference = reference - reference_center
    covariance = centered_mobile.T @ centered_reference
    left, _singular_values, right_t = np.linalg.svd(covariance)
    if np.linalg.det(left @ right_t) < 0:
        left[:, -1] *= -1
    rotation = left @ right_t
    return centered_mobile @ rotation + reference_center
