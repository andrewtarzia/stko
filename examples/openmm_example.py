# ruff: noqa: T201

from pathlib import Path

import stk
import stko
from openff.toolkit import ForceField


def main() -> None:
    """Run the example."""
    bb1 = stk.BuildingBlock("C1CCC(C(C1)N)N", [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock(
        "C1=C(C=C(C=C1C=O)C=O)C=O", [stk.AldehydeFactory()]
    )
    output = Path("openmm_example_directory")
    output.mkdir(exist_ok=True, parents=True)
    cage = stk.ConstructedMolecule(stk.cage.FourPlusSix([bb1, bb2]))
    cage.write(output / "unopt_cage.mol")
    # Load the openff-2.1.0 force field appropriate for
    # vacuum calculations (without constraints)
    force_field = ForceField("openff_unconstrained-2.1.0.offxml")
    optimizer = stko.OpenMMForceField(
        force_field, partial_charges_method="mmff94"
    )
    cage = optimizer.optimize(cage)
    cage.write(output / "opt_cage.mol")


if __name__ == "__main__":
    main()