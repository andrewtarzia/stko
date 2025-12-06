import json
import pathlib

import numpy as np

import stko
from tests.molecular.threesite.case_data import CaseData


def test_threesite_properties(case_data: CaseData) -> None:
    threesite_analysis = stko.molecule_analysis.DitopicThreeSiteAnalyser()

    properties = {
        "binder_adjacent_torsion": (
            threesite_analysis.get_binder_adjacent_torsion(
                case_data.building_block
            )
        ),
        "binder_distance": threesite_analysis.get_binder_distance(
            case_data.building_block
        ),
        "binder_centroid_angle": threesite_analysis.get_binder_centroid_angle(
            case_data.building_block
        ),
        "binder_binder_angle": threesite_analysis.get_binder_binder_angle(
            case_data.building_block
        ),
        "bite_angles": threesite_analysis.get_halfbite_angles(
            case_data.building_block
        ),
        "binder_angles": threesite_analysis.get_binder_angles(
            case_data.building_block
        ),
        "centroids": [
            tuple(i)
            for i in threesite_analysis.get_adjacent_centroids(
                case_data.building_block
            )
        ],
    }

    parent = pathlib.Path(__file__).resolve().parent
    saved_dir = parent / "property_jsons"
    saved_json = saved_dir / f"{case_data.name}_properties.json"
    if not saved_json.exists():
        with saved_json.open("w") as f:
            json.dump(properties, f, indent=4)
        raise AssertionError

    with saved_json.open("r") as f:
        known = json.load(f)

    for prop_name, prop_value in properties.items():
        if prop_name == "centroids":
            test = np.array([list(i) for i in prop_value])
            known_value = np.array(known[prop_name])
            assert np.linalg.norm(test[0] - test[1]) > 0.1  # noqa: PLR2004

        if isinstance(prop_value, float):
            assert np.isclose(prop_value, known[prop_name], atol=1e-3, rtol=0)
        elif isinstance(prop_value, tuple):
            for test, knownval in zip(
                prop_value, known[prop_name], strict=True
            ):
                assert np.isclose(test, knownval, rtol=0, atol=1e-2)
        elif isinstance(prop_value, np.ndarray):
            assert np.isclose(test, known_value, rtol=0, atol=1e-3)
