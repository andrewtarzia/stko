import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            molecule=stk.BuildingBlock("NCCN"),
            zmatrix=(
                "N\nC 1 1.44\nC 2 1.521 54.36\nN 3 1.42 111.18 1 -60.0\n"
                "H 4 3.343 114.09 2 49.63\nH 5 1.714 52.7 3 -138.83\n"
                "H 6 2.485 49.09 4 65.67\nH 7 1.846 88.39 5 -7.67\n"
                "H 8 2.367 85.15 6 -57.42\nH 9 1.818 92.23 7 -59.74\n"
                "H 10 2.919 69.69 8 120.27\nH 11 1.7610 51.39 9 -160.62"
            ),
        ),
        CaseData(
            molecule=stk.BuildingBlock("BrC#CBr"),
            zmatrix=("Br\nC 1 1.89\nC 2 1.221 0.38\nBr 3 1.92 179.61 1 60.23"),
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param
