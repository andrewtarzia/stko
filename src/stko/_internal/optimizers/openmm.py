from openff.interchange import Interchange
from openff.toolkit import ForceField, Molecule
from openmm import Integrator, LangevinIntegrator, State
from openmm.app import Simulation
from openmm.unit import Quantity, kelvin, picosecond, picoseconds

from stko._internal.optimizers.optimizers import Optimizer
from stko._internal.types import MoleculeT


class OpenMMForceField(Optimizer):
    def __init__(
        self,
        force_field: ForceField,
        box_vectors: Quantity | None = None,
        integrator: Integrator | None = None,
        num_steps: int = 10_000,
    ) -> None:
        if integrator is None:
            integrator = LangevinIntegrator(
                300 * kelvin, 1 / picosecond, 0.002 * picoseconds
            )
        self._integrator = integrator
        self._force_field = force_field
        self._box_vectors = box_vectors
        self._num_steps = num_steps

    def optimize(self, mol: MoleculeT) -> MoleculeT:
        molecule = Molecule.from_rdkit(mol)
        topology = molecule.to_topology()
        if self._box_vectors is not None:
            topology.box_vectors = self._box_vectors
        interchange = Interchange.from_smirnoff(self._force_field, topology)
        system = interchange.to_openmm()
        simulation = Simulation(topology, system, self._integrator)
        simulation.minimizeEnergy()
        simulation.step(self._num_steps)
        state = simulation.context.getState(
            getPositions=True,
            getEnergy=True,
        )
        return self._update_stk_molecule(mol, state)

    def _update_stk_molecule(
        self,
        molecule: MoleculeT,
        state: State,
    ) -> MoleculeT:
        positions = state.getPositions(asNumpy=True)
        return molecule.with_position_matrix(positions * 10)
