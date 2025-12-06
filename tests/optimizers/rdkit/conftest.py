from dataclasses import dataclass

import pytest
import stk


@dataclass(frozen=True, slots=True)
class CaseData:
    molecule: stk.Molecule
    unoptimised_energy: float


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            unoptimised_energy=18.487222241888883,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)cc3c2c1"
            ),
            unoptimised_energy=138.73518984157926,
        ),
        CaseData(
            molecule=stk.BuildingBlock("CCCCCC"),
            unoptimised_energy=16.949470583172072,
        ),
        CaseData(
            molecule=stk.BuildingBlock("c1ccccc1"),
            unoptimised_energy=14.760522274728105,
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrCCBr",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles="BrCNCCBr",
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit="AB",
                    num_repeating_units=2,
                    optimizer=stk.MCHammer(),
                ),
            ),
            unoptimised_energy=4890.375436145811,
        ),
    ],
)
def case_uff_molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            unoptimised_energy=26.422861818517934,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)cc3c2c1"
            ),
            unoptimised_energy=119.52858996294171,
        ),
        CaseData(
            molecule=stk.BuildingBlock("CCCCCC"),
            unoptimised_energy=4.036151622820159,
        ),
        CaseData(
            molecule=stk.BuildingBlock("c1ccccc1"),
            unoptimised_energy=18.64619045056115,
        ),
        CaseData(
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles="BrCCBr",
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles="BrCNCCBr",
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit="AB",
                    num_repeating_units=2,
                    optimizer=stk.MCHammer(),
                ),
            ),
            unoptimised_energy=1006.4204295019977,
        ),
    ],
)
def case_mmff_molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param


@pytest.fixture(
    scope="session",
    params=[
        stk.BuildingBlock("NCCN"),
        stk.BuildingBlock(
            "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)cc3c2c1"
        ),
        stk.BuildingBlock("CCCCCC"),
        stk.BuildingBlock("c1ccccc1"),
    ],
)
def case_etkdg_molecule(request: pytest.FixtureRequest) -> stk.BuildingBlock:
    return request.param
