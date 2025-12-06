import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            constructed_molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock("NCCN", [stk.PrimaryAminoFactory()]),
                        stk.BuildingBlock("O=CCCC=O", [stk.AldehydeFactory()]),
                    ),
                    repeating_unit="AB",
                    num_repeating_units=2,
                ),
            ),
            atom_ids={
                0: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                1: [10, 11, 12, 13, 14, 15, 16, 17],
                2: [18, 19, 20, 21, 22, 23, 24, 25, 26, 27],
                3: [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38],
            },
            name=name,
        ),
        lambda name: CaseData(
            constructed_molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock("NCCN", [stk.PrimaryAminoFactory()]),
                        stk.BuildingBlock("O=CCCC=O", [stk.AldehydeFactory()]),
                    ),
                    repeating_unit="AB",
                    num_repeating_units=2,
                    optimizer=stk.Collapser(scale_steps=False),
                ),
            ),
            atom_ids={
                0: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                1: [10, 11, 12, 13, 14, 15, 16, 17],
                2: [18, 19, 20, 21, 22, 23, 24, 25, 26, 27],
                3: [28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38],
            },
            name=name,
        ),
        lambda name: CaseData(
            constructed_molecule=stk.ConstructedMolecule(
                topology_graph=stk.cage.FourPlusSix(
                    (
                        stk.BuildingBlock(
                            smiles="NCCN",
                            functional_groups=[stk.PrimaryAminoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles="O=CC(C=O)C=O",
                            functional_groups=[stk.AldehydeFactory()],
                        ),
                    )
                ),
            ),
            atom_ids={
                0: [0, 1, 2, 3, 4, 5, 6, 7],
                1: [8, 9, 10, 11, 12, 13, 14, 15],
                2: [16, 17, 18, 19, 20, 21, 22, 23],
                3: [24, 25, 26, 27, 28, 29, 30, 31],
                4: [32, 33, 34, 35, 36, 37, 38, 39],
                5: [40, 41, 42, 43, 44, 45, 46, 47],
                6: [48, 49, 50, 51, 52, 53, 54, 55],
                7: [56, 57, 58, 59, 60, 61, 62, 63],
                8: [64, 65, 66, 67, 68, 69, 70, 71],
                9: [72, 73, 74, 75, 76, 77, 78, 79],
            },
            name=name,
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
