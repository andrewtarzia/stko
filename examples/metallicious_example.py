# ruff: noqa: T201

import logging
from pathlib import Path

import stk
from metallicious import (
    prepare_stk_initial_topology,
    stk_supramolecular_structure,
)
from openff.toolkit import ForceField

import stko


def main() -> None:
    """Run the example."""
    output = Path("metallicious_output")
    output.mkdir(exist_ok=True, parents=True)

    ligand_top = output / "ligand.top"
    ligand_coord = output / "ligand.pdb"

    # Build a MOC.
    if not ligand_top.exists():
        pd = stk.BuildingBlock(
            smiles="[Pd+2]",
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2)) for i in range(4)
            ),
            position_matrix=[[0.0, 0.0, 0.0]],
        )
        ditopic_bb = stk.BuildingBlock(
            smiles="C1=NC=CC(C2=CC=CC(C3=CC=NC=C3)=C2)=C1",
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts="[#6]~[#7X2]~[#6]", bonders=(1,), deleters=()
                )
            ],
        )

        prepare_stk_initial_topology(
            molecule=ditopic_bb,
            metal_names=[],
            metal_charge=0,
            output_dir=output,
            output_coord=ligand_coord,
            output_top=ligand_top,
            metal_vdw="uff",
            ligand_topol=None,
        )

        apdcage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M6L12Cube(
                building_blocks=(pd, ditopic_bb),
                reaction_factory=stk.DativeReactionFactory(
                    stk.GenericReactionFactory(
                        bond_orders={
                            frozenset(
                                {
                                    stk.GenericFunctionalGroup,
                                    stk.SingleAtom,
                                }
                            ): 9,
                        },
                    ),
                ),
                optimizer=stk.MCHammer(num_steps=1500),
            ),
        )
        apdcage.write(output / "cage_unopt.mol")
        apdcage.write(output / "cage_unopt.pdb")

    out_coord = output / "out.pdb"
    out_topol = output / "out.top"
    if not out_coord.exists() or not out_topol.exists():
        # Run metallicious.
        cage = stk_supramolecular_structure(
            filename=str(output / "p6_cage.pdb"),
            molecule=stk.BuildingBlock.init_from_file(
                output / "cage_unopt.mol"
            ),
            metal_charges={"Pd": 2},
            LJ_type="uff",
            donors=["N"],
        )

        cage.parametrize(
            molecule=stk.BuildingBlock.init_from_file(
                output / "cage_unopt.mol"
            ),
            out_coord=out_coord,
            out_topol=out_topol,
            prepare_initial_topology=True,
            homoleptic_ligand_topol=ligand_top,
            sub_dir=output / "init_top",
        )

    # Optimise with OpenMM.

    raise SystemExit

    ff_optimizer = stko.OpenMMForceField(
        # Load the openff-2.1.0 force field appropriate for
        # vacuum calculations (without constraints)
        force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
        partial_charges_method="espaloma-am1bcc",
        max_iterations=10,
    )
    ff_cage = ff_optimizer.optimize(cage)
    ff_cage.write(output / "ff_opt_cage.mol")
    print(
        "ff_opt_cage",
        stko.OpenMMEnergy(
            force_field=ForceField("openff_unconstrained-2.1.0.offxml"),
            partial_charges_method="espaloma-am1bcc",
        ).get_energy(ff_cage),
    )


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
