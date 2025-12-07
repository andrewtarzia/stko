import json
import pathlib

import numpy as np

import stko

from .case_data import CaseData


def test_geometry_properties(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.GeometryAnalyser()

    properties = {
        "radius_gyration": analyser.get_radius_gyration(case_data.molecule),
        "avg_centoid_distance": analyser.get_avg_centroid_distance(
            case_data.molecule
        ),
        "max_diameter": analyser.get_max_diameter(case_data.molecule),
        "metal_centroid_angles": {
            f"{i[0]}_{i[1]}": j
            for i, j in analyser.get_metal_centroid_metal_angle(
                case_data.molecule,
                metal_atom_nos=(26, 46),
            ).items()
        },
        "metal_atom_distances": {
            f"{i[0]}_{i[1]}": j
            for i, j in analyser.get_metal_distances(
                case_data.molecule,
                metal_atom_nos=(26, 46),
            ).items()
        },
        "min_atom_atom_distance": analyser.get_min_atom_atom_distance(
            case_data.molecule
        ),
        "min_centoid_distance": analyser.get_min_centroid_distance(
            case_data.molecule
        ),
    }

    parent = pathlib.Path(__file__).resolve().parent
    saved_dir = parent / "geom_jsons"
    saved_json = saved_dir / f"{case_data.name}_properties.json"
    if not saved_json.exists():
        with saved_json.open("w") as f:
            json.dump(properties, f, indent=4)
        raise AssertionError

    with saved_json.open("r") as f:
        known = json.load(f)

    for prop_name, prop_value in properties.items():
        if isinstance(prop_value, float):
            assert np.isclose(prop_value, known[prop_name], atol=1e-3, rtol=0)
        elif isinstance(prop_value, tuple):
            for test, knownval in zip(
                prop_value, known[prop_name], strict=True
            ):
                assert np.isclose(test, knownval, rtol=0, atol=1e-2)
        elif isinstance(prop_value, dict):
            for i in prop_value:
                assert np.isclose(
                    prop_value[i], known[prop_name][i], rtol=0, atol=1e-2
                )
