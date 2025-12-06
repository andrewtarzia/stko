from dataclasses import dataclass

import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    building_block: stk.BuildingBlock
    name: str
