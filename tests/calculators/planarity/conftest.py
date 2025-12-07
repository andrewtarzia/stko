from dataclasses import dataclass

import numpy as np
import pytest
import stk

import stko

_macrocycle = stk.ConstructedMolecule(
    topology_graph=stk.macrocycle.Macrocycle(
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
    ),
)

_square_planar = stk.ConstructedMolecule(
    topology_graph=stk.metal_complex.SquarePlanar(
        metals=stk.BuildingBlock(
            smiles="[Pd+2]",
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2)) for _ in range(4)
            ),
            position_matrix=np.array([[0, 0, 0]]),
        ),
        ligands=stk.BuildingBlock(
            smiles="NBr",
            functional_groups=(stk.PrimaryAminoFactory(),),
        ),
        optimizer=stk.MCHammer(),
    )
)
uff = stko.UFF()
_square_planar_uff = uff.optimize(_square_planar)
_octahedral = stk.ConstructedMolecule(
    topology_graph=stk.metal_complex.Octahedral(
        metals=stk.BuildingBlock(
            smiles="[Fe+2]",
            functional_groups=(
                stk.SingleAtom(stk.Fe(0, charge=2)) for _ in range(6)
            ),
            position_matrix=np.array([[0, 0, 0]]),
        ),
        ligands=stk.BuildingBlock(
            smiles="NBr",
            functional_groups=(stk.PrimaryAminoFactory(),),
        ),
        optimizer=stk.MCHammer(),
    )
)


@dataclass(frozen=True, slots=True)
class CaseData:
    molecule: stk.Molecule
    plane_ids: tuple[int, ...] | None
    deviation_ids: tuple[int, ...] | None
    plane_deviation: float
    planarity_parameter: float
    plane_span: float


@pytest.fixture(
    scope="session",
    params=[
        CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=7.0468262994917374,
            plane_span=2.8619224646264096,
            planarity_parameter=0.7376285949589063,
        ),
        CaseData(
            molecule=stk.BuildingBlock(
                "C(#Cc1cccc2ccncc21)c1ccc2[nH]c3ccc(C#Cc4cccc5cnccc54)cc3c2c1"
            ),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=30.72874008769579,
            plane_span=4.190622373823865,
            planarity_parameter=0.7774786041256974,
        ),
        CaseData(
            molecule=stk.BuildingBlock("CCCCCC"),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=10.798240107278682,
            plane_span=2.833250157796366,
            planarity_parameter=0.7469082286580176,
        ),
        CaseData(
            molecule=stk.BuildingBlock("c1ccccc1"),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            plane_span=0.0,
            planarity_parameter=0.0,
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
                ),
            ),
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=16.368844071942494,
            plane_span=3.108677868237704,
            planarity_parameter=0.6376957812606647,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=20.467847459970862,
            plane_span=2.667602762686369,
            planarity_parameter=0.7279960670882043,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=tuple(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6  # noqa: PLR2004
            ),
            deviation_ids=tuple(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6  # noqa: PLR2004
            ),
            plane_deviation=1.6024957200580106,
            plane_span=0.5007799125181283,
            planarity_parameter=0.20031196500725132,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=tuple(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6  # noqa: PLR2004
            ),
            deviation_ids=None,
            plane_deviation=20.33209530254546,
            plane_span=2.667602762686369,
            planarity_parameter=0.7283475685355212,
        ),
        CaseData(
            molecule=_macrocycle,
            plane_ids=None,
            deviation_ids=tuple(
                i.get_id()
                for i in _macrocycle.get_atoms()
                if i.get_atomic_number() == 6  # noqa: PLR2004
            ),
            plane_deviation=1.7382478774834076,
            plane_span=0.5007799125181284,
            planarity_parameter=0.20158568952378236,
        ),
        CaseData(
            molecule=_square_planar,
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            plane_span=0.0,
            planarity_parameter=0.0,
        ),
        CaseData(
            molecule=_square_planar_uff,
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=0.0,
            plane_span=0.0,
            planarity_parameter=0.0,
        ),
        CaseData(
            molecule=_octahedral,
            plane_ids=None,
            deviation_ids=None,
            plane_deviation=15.350980683338454,
            plane_span=5.734353982049391,
            planarity_parameter=1.4635207224067897,
        ),
    ],
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param
