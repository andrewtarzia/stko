import logging
from itertools import combinations

import numpy as np
import rdkit.Chem.AllChem as rdkit  # noqa: N813
import stk

from stko._internal.internal_types import MoleculeT
from stko._internal.optimizers.optimizers import Optimizer
from stko._internal.optimizers.utilities import (
    get_metal_atoms,
    get_metal_bonds,
    to_rdkit_mol_without_metals,
)
from stko._internal.utilities.utilities import vector_angle

logger = logging.getLogger(__name__)


class MMFF(Optimizer):
    """Use the MMFF force field in :mod:`rdkit` to optimize molecules.

    See Also:
        * rdkit: https://www.rdkit.org/

    .. warning::
        this optimizer seems to be machine dependant, producing
        different energies after optimisation on Ubunut 18 vs. Ubuntu 20.

    Examples:
        .. testcode:: rdkit-mmff

            import stk
            import stko

            mol = stk.BuildingBlock('NCCNCCN')
            mmff = stko.MMFF()
            opt_mol = mmff.optimize(mol)

        .. testcode:: rdkit-mmff
            :hide:

            assert (
                stk.Smiles().get_key(mol) == stk.Smiles().get_key(opt_mol)
            )
            opt_e = stko.MMFFEnergy().get_energy(opt_mol)
            assert opt_e < stko.MMFFEnergy().get_energy(mol)

    """

    def __init__(self, ignore_inter_interactions: bool = True) -> None:
        self._ignore_inter_interactions = ignore_inter_interactions

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.MMFFOptimizeMolecule(
            rdkit_mol,
            ignoreInterfragInteractions=self._ignore_inter_interactions,
        )
        return mol.with_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )


class UFF(Optimizer):
    """Use the UFF force field in :mod:`rdkit` to optimize molecules.

    .. warning::
        this optimizer seems to be machine dependant, producing
        different energies after optimisation on Ubunut 18 vs. Ubuntu 20.

    See Also:
        * rdkit: https://www.rdkit.org/

    Examples:
        .. testcode:: rdkit-uff

            import stk
            import stko

            mol = stk.BuildingBlock('NCCNCCN')
            uff = stko.UFF()
            opt_mol = uff.optimize(mol)

        .. testcode:: rdkit-uff
            :hide:

            assert (
                stk.Smiles().get_key(mol) == stk.Smiles().get_key(opt_mol)
            )
            opt_e = stko.UFFEnergy().get_energy(opt_mol)
            assert opt_e < stko.UFFEnergy().get_energy(mol)

    """

    def __init__(self, ignore_inter_interactions: bool = True) -> None:
        self._ignore_inter_interactions = ignore_inter_interactions

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.UFFOptimizeMolecule(
            rdkit_mol,
            ignoreInterfragInteractions=self._ignore_inter_interactions,
        )
        return mol.with_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )


class ETKDG(Optimizer):
    """Uses ETKDG v2 algorithm in :mod:`rdkit` to optimize a structure.

    See Also:
        * ETKDG: http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00654
        * rdkit: https://www.rdkit.org/

    Parameters:
        random_seed:
            The random seed to use.

    Examples:
        .. testcode:: rdkit-etkdg

            import stk
            import stko

            mol = stk.BuildingBlock('NCCNCCN')
            etkdg = stko.ETKDG()
            opt_mol = etkdg.optimize(mol)

        .. testcode:: rdkit-etkdg
            :hide:

            assert (
                stk.Smiles().get_key(mol) == stk.Smiles().get_key(opt_mol)
            )

    """

    def __init__(self, random_seed: int = 12) -> None:
        self._random_seed = random_seed

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        params = rdkit.ETKDGv2()
        params.clearConfs = True
        params.random_seed = self._random_seed

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.EmbedMolecule(rdkit_mol, params)
        return mol.with_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )


class MetalOptimizer(Optimizer):
    """Applies forcefield optimizers in :mod:`rdkit` that can handle metals.

    Parameters:
        metal_binder_distance:
            Distance in Angstrom.

        metal_binder_forceconstant:
            Force constant to use for restricted metal-ligand bonds.

        max_iterations:
            Number of iteractions to run.


    See Also:
        * rdkit: https://www.rdkit.org/

    Notes:
        By default, :meth:`optimize` will run a restricted optimization
        using constraints and the UFF. To implement this, metal atoms are
        replaced by noninteracting H atoms, and constraints are applied
        to maintain the metal centre geometry.
        Restrictions are applied to the ligand with respect to its input
        structure. So if that is poorly optimised, then the output will be
        also.

        Warning: this optimizer seems to be machine dependant, producing
        different energies after optimisation on Ubunut 18 vs. Ubuntu 20.

    Examples:
        :class:`MetalOptimizer` allows for the restricted optimization of
        :class:`ConstructedMolecule` instances containing metals. Note that
        this optimizer algorithm is not very robust to large bonds and may
        fail.

        .. testcode:: rdkit-metal-opt

            import stk
            import stko

            # Produce a Pd+2 atom with 4 functional groups.
            palladium_atom = stk.BuildingBlock(
                smiles='[Pd+2]',
                functional_groups=(
                    stk.SingleAtom(stk.Pd(0, charge=2))
                    for i in range(4)
                ),
                position_matrix=[[0., 0., 0.]],
            )

            # Build a building block with two functional groups using
            # the SmartsFunctionalGroupFactory.
            bb1 = stk.BuildingBlock(
                smiles=(
                    'C1=NC=CC(C2=CC=CC(C3=C'
                    'C=NC=C3)=C2)=C1'
                ),
                functional_groups=[
                    stk.SmartsFunctionalGroupFactory(
                        smarts='[#6]~[#7X2]~[#6]',
                        bonders=(1, ),
                        deleters=(),
                    ),
                ],
            )

            cage1 = stk.ConstructedMolecule(
                stk.cage.M3L6(
                    building_blocks=(palladium_atom, bb1),
                    # Ensure that bonds between the GenericFunctionalGroups
                    # of the ligand and the SingleAtom functional groups
                    # of the metal are dative.
                    reaction_factory=stk.DativeReactionFactory(
                        stk.GenericReactionFactory(
                            bond_orders={
                                frozenset({
                                    stk.GenericFunctionalGroup,
                                    stk.SingleAtom
                                }): 9
                            }
                        )
                    ),
                    optimizer=stk.MCHammer(),
                )
            )

            # Define an optimizer.
            optimizer = stko.MetalOptimizer()

            # Optimize.
            opt_cage1 = optimizer.optimize(mol=cage1)

        .. testcode:: rdkit-metal-opt
            :hide:

            assert (
                stk.Smiles().get_key(cage1) == stk.Smiles().get_key(opt_cage1)
            )


    """

    def __init__(
        self,
        metal_binder_distance: float = 1.6,
        metal_binder_forceconstant: float = 1.0e2,
        max_iterations: int = 500,
    ) -> None:
        self._metal_binder_distance = metal_binder_distance
        self._metal_binder_forceconstant = metal_binder_forceconstant
        self._max_iterations = max_iterations

    def _apply_metal_centre_constraints(
        self,
        mol: stk.Molecule,
        ff: rdkit.ForceField,
        metal_bonds: list[stk.Bond],
    ) -> None:
        """Applies UFF metal centre constraints.

        Parameters:
            mol:
                The molecule to be optimized.

            ff:
                Forcefield to apply constraints to. Generally use UFF.

            metal_bonds:
                List of bonds including metal atoms.

        """
        # Add constraints to UFF to hold metal geometry in place.
        for bond in metal_bonds:
            idx1 = bond.get_atom1().get_id()
            idx2 = bond.get_atom2().get_id()
            # Add distance constraints in place of metal bonds.
            # Target distance set to a given metal_binder_distance.
            ff.UFFAddDistanceConstraint(
                idx1=idx1,
                idx2=idx2,
                relative=False,
                minLen=self._metal_binder_distance,
                maxLen=self._metal_binder_distance,
                forceConstant=self._metal_binder_forceconstant,
            )

        # Also implement angular constraints to all atoms in the
        # metal complex.
        for bonds in combinations(metal_bonds, r=2):
            bond1, bond2 = bonds
            bond1_atoms = [bond1.get_atom1(), bond1.get_atom2()]
            bond2_atoms = [bond2.get_atom1(), bond2.get_atom2()]
            pres_atoms = list(set(bond1_atoms + bond2_atoms))
            # If there are more than 3 atoms, implies two
            # independant bonds.
            if len(pres_atoms) > 3:  # noqa: PLR2004
                continue
            for atom in pres_atoms:
                if atom in bond1_atoms and atom in bond2_atoms:
                    idx2 = atom.get_id()
                elif atom in bond1_atoms:
                    idx1 = atom.get_id()
                elif atom in bond2_atoms:
                    idx3 = atom.get_id()
            pos1 = next(mol.get_atomic_positions(atom_ids=[idx1]))
            pos2 = next(mol.get_atomic_positions(atom_ids=[idx2]))
            pos3 = next(mol.get_atomic_positions(atom_ids=[idx3]))
            v1 = pos1 - pos2
            v2 = pos3 - pos2
            angle = vector_angle(v1, v2)
            ff.UFFAddAngleConstraint(
                idx1=idx1,
                idx2=idx2,
                idx3=idx3,
                relative=False,
                minAngleDeg=np.degrees(angle),
                maxAngleDeg=np.degrees(angle),
                forceConstant=1.0e5,
            )

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        # Find all metal atoms and atoms they are bonded to.
        metal_atoms = get_metal_atoms(mol)
        metal_bonds, ids_to_metals = get_metal_bonds(
            mol=mol, metal_atoms=metal_atoms
        )

        # Perform a forcefield optimisation that
        # only optimises non metal atoms that are not bonded to the
        # metal.

        # Write rdkit molecule with metal atoms and bonds deleted.
        edit_mol = to_rdkit_mol_without_metals(
            mol=mol, metal_atoms=metal_atoms, metal_bonds=metal_bonds
        )

        # Non-bonded interactions need to be explicitly turned on (if
        # desired) because at the point of initialisation, the
        # molecules are technically separate (not bonded) and are
        # treated as fragments.
        rdkit.SanitizeMol(edit_mol)
        ff = rdkit.UFFGetMoleculeForceField(
            edit_mol,
            ignoreInterfragInteractions=False,
        )

        # Constrain the metal centre.
        self._apply_metal_centre_constraints(
            mol=mol,
            ff=ff,
            metal_bonds=metal_bonds,
        )

        # Optimisation with UFF in RDKit. This method uses constraints
        # on the metal centre to attempt to enforce the metal geometry
        # described by the metal topology.
        ff.Minimize(maxIts=self._max_iterations)

        # Update stk molecule from optimized molecule. This should
        # only modify atom positions, which means metal atoms will be
        # reinstated.
        new_position_matrix = edit_mol.GetConformer().GetPositions()
        return mol.with_position_matrix(new_position_matrix)
