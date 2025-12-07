import json
import pathlib

import numpy as np

import stko

from .case_data import CaseData


def test_calculate_torsions(case_data: CaseData) -> None:
    analyser = stko.molecule_analysis.GeometryAnalyser()

    result = analyser.calculate_torsions(case_data.molecule)

    parent = pathlib.Path(__file__).resolve().parent
    saved_dir = parent / "geom_jsons"
    saved_json = saved_dir / f"{case_data.name}_torsions.json"
    if not saved_json.exists():
        with saved_json.open("w") as f:
            json.dump({"_".join(i): j for i, j in result.items()}, f, indent=4)
        raise AssertionError

    with saved_json.open("r") as f:
        known = {tuple(i.split("_")): j for i, j in json.load(f).items()}

    for four in result:
        assert len(result[four]) == len(known[four])
        assert np.allclose(result[four], known[four], rtol=0, atol=1e-3)
