import json
import pathlib

import numpy as np

import stko
from tests.molecular.constructed.case_data import CaseData


def test_get_building_block_centroids(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.ConstructedAnalyser()
    result = {
        i: list(j)
        for i, j in analyser.get_building_block_centroids(
            case_data.constructed_molecule
        ).items()
    }

    parent = pathlib.Path(__file__).resolve().parent
    saved_dir = parent / "centroid_jsons"
    saved_json = saved_dir / f"{case_data.name}.json"
    if not saved_json.exists():
        with saved_json.open("w") as f:
            json.dump(result, f, indent=4)
        raise AssertionError

    with saved_json.open("r") as f:
        known = json.load(f)

    assert len(result) == len(known)

    for i, test in result.items():
        assert np.allclose(test, known[str(i)], rtol=0, atol=1e-3)
