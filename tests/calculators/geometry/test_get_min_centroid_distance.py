import numpy as np

import stko

from .case_data import CaseData


def test_get_min_centroid_distance(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.get_min_centroid_distance(case_data.molecule)
    assert np.isclose(
        result, case_data.min_centoid_distance, atol=1e-3, rtol=0
    )
