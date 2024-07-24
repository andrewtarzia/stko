import gzip
import re
import shutil
from collections import abc, deque
from itertools import chain
from pathlib import Path

import rdkit.Chem.AllChem as rdkit  # noqa: N813
import stk
from rdkit.Geometry import Point3D

# This dictionary gives easy access to the rdkit bond types.
bond_dict = {
    "1": rdkit.rdchem.BondType.SINGLE,
    "am": rdkit.rdchem.BondType.SINGLE,
    "2": rdkit.rdchem.BondType.DOUBLE,
    "3": rdkit.rdchem.BondType.TRIPLE,
    "ar": rdkit.rdchem.BondType.AROMATIC,
}


class MAEExtractor:
    """Extracts the lowest energy conformer from a .maegz file.

    Macromodel conformer searches produce -out.maegz files containing
    all of the conformers found during the search and their energies
    and other data.

    Initializing this class with a :class:`.ConstructedMolecule` finds
    the ``-out.maegz`` file of that :class:`.ConstructedMolecule` and
    converts it to a ``.mae`` file. It then creates and additional
    ``.mae`` file holding only the lowest energy conformer found.

    """

    maegz_path: Path
    """``-out.maegz`` file generated by the conformer search."""

    mae_path: Path
    """``.mae`` file holding the conformers from the conformer search."""

    content: str
    """The content of the ``.mae`` file holding all the conformers.

    This holds other data such as their energies too.

    """

    energies: list[tuple[float | None, int]]
    """Energies of the lowest energy confromers and their id."""

    min_energy: float | None
    """The minimum energy found in the ``.mae`` file."""

    path: Path
    """``.mae`` file holding the extracted lowest energy conformer."""

    def __init__(self, run_name: str, n: int = 1) -> None:
        self.maegz_path = Path(f"{run_name}-out.maegz")
        self.maegz_to_mae()
        self.extract_conformers(n)

    def extract_conformers(self, n: int) -> None:
        """Creates ``.mae`` files holding the lowest energy conformers.

        Parameters:
            n:
                The number of conformers to extract.

        """
        for i in range(n):
            # Get the id of the lowest energy conformer.
            num = self.lowest_energy_conformers(n)[i][1]
            # Get the structure block corresponding to the lowest
            # energy conformer.
            content = self.content.split("f_m_ct")
            new_mae = "f_m_ct".join([content[0], content[num]])

            # Write the structure block in its own .mae file, named
            # after conformer extracted.
            if n == 1:
                # Write the structure block in its own .mae file, named
                # after conformer extracted.
                new_name = self.mae_path.with_name(
                    f"{self.mae_path.stem}EXTRACTED_{num}.mae"
                )
            else:
                new_name = self.mae_path.with_name(
                    f"{self.mae_path.stem}EXTRACTED_{num}_conf_{i}.mae"
                )

            with new_name.open("w") as mae_file:
                mae_file.write(new_mae)

            if i == 0:
                # Save the path of the newly created file.
                self.path = new_name

    def extract_energy(self, block: str) -> float | None:
        """Extracts the energy value from a ``.mae`` energy data block.

        Parameters:
            block:
                An ``.mae`` energy data block.

        Returns:
            The energy value extracted from `block` or ``None`` if
            one is not found.

        """
        block_list = block.split(":::")
        for name, value in zip(
            block_list[0].split("\n"), block_list[1].split("\n"), strict=False
        ):
            if "r_mmod_Potential_Energy" in name:
                return float(value)
        return None

    def lowest_energy_conformers(
        self,
        n: int,
    ) -> list[tuple[float | None, int]]:
        """Returns the id and energy of the lowest energy conformers.

        Parameters:
            n:
                The number of lowest energy conformers to return.

        Returns:
            A :class:`list` of the form

            .. code-block:: python

                returned = [(123.3, 23), (143.89, 1), (150.6, 12), ...]

            Where each :class:`tuple` holds the energy of the
            `n` lowest energy conformers and the id, respectively.

        """
        # Open the .mae file holding all the conformers and load its
        # content.
        self.content = self.mae_path.read_text()
        # Split the content across curly braces. This divides the
        # various sections of the .mae file.
        content_split = re.split(r"[{}]", self.content)

        # Go through all the datablocks in the the .mae file. For each
        # energy block extract the energy and store it in the
        # `energies` list. Store the `index`  (conformer id) along with
        # each extracted energy.
        self.energies = []
        prev_block = deque([""], maxlen=1)
        index = 1
        for block in content_split:
            if (
                "f_m_ct" in prev_block[0]
                and "r_mmod_Potential_Energy" in block
            ):
                energy = self.extract_energy(block)
                self.energies.append((energy, index))
                index += 1

            prev_block.append(block)

        # Selecting the lowest energy n conformers
        confs = sorted(self.energies)[:n]
        # Define the energy of the lowest energy conformer
        self.min_energy = confs[0][0]
        # Return a list with id and energy of the lowest energy
        # conformers.
        return confs

    def maegz_to_mae(self) -> None:
        """Converts the .maegz file to a .mae file."""
        self.mae_path = self.maegz_path.with_suffix(".mae")
        with (
            gzip.open(self.maegz_path, "r") as maegz_file,
            self.mae_path.open("wb") as mae_file,
        ):
            mae_file.write(maegz_file.read())


def mol_from_mae_file(mae_path: Path | str) -> rdkit.Mol:  # noqa: PLR0912, C901
    """Creates a ``rdkit`` molecule from a ``.mae`` file.

    Parameters:
        mol2_file:
            The full path of the ``.mae`` file from which an rdkit molecule
            should be instantiated.

    Returns:
        An ``rdkit`` instance of the molecule held in `mae_file`.

    """
    mae_path = Path(mae_path)

    mol = rdkit.EditableMol(rdkit.Mol())
    conf = rdkit.Conformer()

    content = re.split(r"[{}]", mae_path.read_text())

    prev_block = deque([""], maxlen=1)
    for block in content:
        if "m_atom[" in prev_block[0]:
            atom_block = block
        if "m_bond[" in prev_block[0]:
            bond_block = block
        prev_block.append(block)

    labels, data_block, *_ = atom_block.split(":::")
    labels = [
        label
        for label in labels.split("\n")
        if not label.isspace() and label != ""
    ]

    data_block = [
        a.split()
        for a in data_block.split("\n")
        if not a.isspace() and a != ""
    ]

    for line in data_block:
        words = [word for word in line if word != '"']
        if len(labels) != len(words):
            msg = (
                "Number of labels does"
                " not match number of columns"
                " in .mae file."
            )
            raise RuntimeError(msg)

        for label, data in zip(labels, words, strict=False):
            if "x_coord" in label:
                x = float(data)
            if "y_coord" in label:
                y = float(data)
            if "z_coord" in label:
                z = float(data)
            if "atomic_number" in label:
                atom_num = int(data)

        atom_sym = rdkit.GetPeriodicTable().GetElementSymbol(atom_num)
        atom_coord = Point3D(x, y, z)
        atom_id = mol.AddAtom(rdkit.Atom(atom_sym))
        conf.SetAtomPosition(atom_id, atom_coord)

    labels, data_block, *_ = bond_block.split(":::")
    labels = [
        label
        for label in labels.split("\n")
        if not label.isspace() and label != ""
    ]
    data_block = [
        a.split()
        for a in data_block.split("\n")
        if not a.isspace() and a != ""
    ]

    for line in data_block:
        if len(labels) != len(line):
            msg = (
                "Number of labels does"
                " not match number of "
                "columns in .mae file."
            )
            raise RuntimeError(msg)

        for label, data in zip(labels, line, strict=False):
            if "from" in label:
                atom1 = int(data) - 1
            if "to" in label:
                atom2 = int(data) - 1
            if "order" in label:
                bond_order = str(int(data))
        mol.AddBond(atom1, atom2, bond_dict[bond_order])

    mol = mol.GetMol()
    mol.AddConformer(conf)
    return mol


def move_generated_macromodel_files(
    basename: str, output_dir: Path | str
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for filename in Path().glob(f"{basename}*"):
        # Do not move the output_dir.
        if filename == output_dir:
            continue
        shutil.move(filename, output_dir / filename)


def has_h_atom(bond: stk.Bond) -> bool:
    """Check if a bond has a H atom.

    Parameters:
        bond:
            Bond to test if it has a H atom.

    Returns:
        Returns `True` if bond has H atom.

    """
    if bond.get_atom1().get_atomic_number() == 1:
        return True
    return bond.get_atom2().get_atomic_number() == 1


def has_metal_atom(bond: stk.Bond, metal_atoms: list[stk.Atom]) -> bool:
    """Check if a bond has a metal atom.

    Parameters:
        bond:
            Bond to test if it has a metal atom.

        metal_atoms:
            List of metal atoms.

    Returns:
        Returns `True` if bond has metal atom.

    """
    if bond.get_atom1() in metal_atoms:
        return True
    return bond.get_atom2() in metal_atoms


def metal_atomic_numbers() -> abc.Iterable:
    return chain(range(21, 31), range(39, 49), range(72, 81))


def get_metal_atoms(mol: stk.Molecule) -> list[stk.Atom]:
    """Return a list of metal atoms in molecule."""
    metals = set(metal_atomic_numbers())
    return [
        atom for atom in mol.get_atoms() if atom.get_atomic_number() in metals
    ]


def get_metal_bonds(
    mol: stk.Molecule,
    metal_atoms: list[stk.Atom],
) -> tuple[list, list]:
    """Return a list of bonds in molecule that contain metal atoms."""
    metal_bonds = []
    ids_to_metals = []
    for bond in mol.get_bonds():
        if bond.get_atom1() in metal_atoms:
            metal_bonds.append(bond)
            ids_to_metals.append(bond.get_atom2().get_id())
        elif bond.get_atom2() in metal_atoms:
            metal_bonds.append(bond)
            ids_to_metals.append(bond.get_atom1().get_id())

    return metal_bonds, ids_to_metals


def to_rdkit_mol_without_metals(
    mol: stk.Molecule,
    metal_atoms: list[stk.Atom],
    metal_bonds: list[stk.Bond],
) -> rdkit.Mol:
    """Create :class:`rdkit.Mol` with metals replaced by H atoms.

    Parameters:
        mol:
            The molecule to be optimized.

        metal_atoms:
            List of metal atoms.

        metal_bonds:
            List of bonds including metal atoms.

    Returns:
        RDKit molecule with metal atoms replaced with H atoms.

    """
    edit_mol = rdkit.EditableMol(rdkit.Mol())
    for atom in mol.get_atoms():
        if atom in metal_atoms:
            # In place of metals, add H's that will be constrained.
            # This allows the atom ids to not be changed.
            rdkit_atom = rdkit.Atom(1)
            rdkit_atom.SetFormalCharge(0)
        else:
            rdkit_atom = rdkit.Atom(atom.get_atomic_number())
            rdkit_atom.SetFormalCharge(atom.get_charge())
        edit_mol.AddAtom(rdkit_atom)

    for bond in mol.get_bonds():
        if bond in metal_bonds:
            # Do not add bonds to metal atoms (replaced with H's).
            continue
        edit_mol.AddBond(
            beginAtomIdx=bond.get_atom1().get_id(),
            endAtomIdx=bond.get_atom2().get_id(),
            order=rdkit.BondType(bond.get_order()),
        )

    edit_mol = edit_mol.GetMol()
    rdkit_conf = rdkit.Conformer(mol.get_num_atoms())
    for atom_id, atom_coord in enumerate(mol.get_position_matrix()):
        rdkit_conf.SetAtomPosition(atom_id, atom_coord)
        edit_mol.GetAtomWithIdx(atom_id).SetNoImplicit(True)  # noqa:FBT003
    edit_mol.AddConformer(rdkit_conf)

    return edit_mol


def get_long_bond_ids(
    mol: stk.ConstructedMolecule,
    reorder: bool = False,
) -> tuple[tuple[int, int], ...]:
    """Return tuple of long bond ids in a ConstructedMolecule."""
    long_bond_ids = []
    for bond_infos in mol.get_bond_infos():
        if bond_infos.get_building_block() is None:
            if reorder:
                ba1 = bond_infos.get_bond().get_atom1().get_id()
                ba2 = bond_infos.get_bond().get_atom2().get_id()
                if ba1 < ba2:
                    ids = (
                        bond_infos.get_bond().get_atom1().get_id(),
                        bond_infos.get_bond().get_atom2().get_id(),
                    )
                else:
                    ids = (
                        bond_infos.get_bond().get_atom2().get_id(),
                        bond_infos.get_bond().get_atom1().get_id(),
                    )
            else:
                ids = (
                    bond_infos.get_bond().get_atom1().get_id(),
                    bond_infos.get_bond().get_atom2().get_id(),
                )
            long_bond_ids.append(ids)

    return tuple(long_bond_ids)
