import pytest
import stk

import stko

from .case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda name: CaseData(
            building_block=stk.BuildingBlock(
                smiles="C1=CN=CC=N1",
                functional_groups=stko.functional_groups.ThreeSiteFactory(
                    smarts="[#6]~[#7X2]~[#6]"
                ),
            ),
            name=name,
        ),
        lambda name: CaseData(
            building_block=stk.BuildingBlock(
                smiles="C1=CC(=CN=C1)C2=CC=C(C=C2)C3=CN=CC=C3",
                functional_groups=stko.functional_groups.ThreeSiteFactory(
                    smarts="[#6]~[#7X2]~[#6]"
                ),
            ),
            name=name,
        ),
        lambda name: CaseData(
            building_block=stk.BuildingBlock(
                smiles="C1=CC(=CC(=C1)C2=CC=NC=C2)C3=CC=NC=C3",
                functional_groups=stko.functional_groups.ThreeSiteFactory(
                    smarts="[#6]~[#7X2]~[#6]"
                ),
            ),
            name=name,
        ),
        lambda name: CaseData(
            building_block=stk.BuildingBlock(
                smiles="C1C=C(C2=CC3=C(OC4=C3C=C(C3C=CN=CC=3)C=C4)C=C2)C=CN=1",
                functional_groups=stko.functional_groups.ThreeSiteFactory(
                    smarts="[#6]~[#7X2]~[#6]"
                ),
            ),
            name=name,
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
