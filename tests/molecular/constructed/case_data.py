from dataclasses import dataclass

import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    constructed_molecule: stk.ConstructedMolecule
    atom_ids: dict[int, list[int]]
    name: str
