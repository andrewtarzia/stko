import stk
import stko


def main():
    bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
    bb2 = stk.BuildingBlock(
        smiles='O=CC(C=O)C=O',
        functional_groups=[stk.AldehydeFactory()],
    )
    cage1 = stk.ConstructedMolecule(
        topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
    )

    broken_bonds_by_id = []
    disconnectors = []
    for bi in cage1.get_bond_infos():
        if bi.get_building_block() is None:
            a1id = bi.get_bond().get_atom1().get_id()
            a2id = bi.get_bond().get_atom2().get_id()
            broken_bonds_by_id.append(sorted((a1id, a2id)))
            disconnectors.extend((a1id, a2id))

    print(broken_bonds_by_id)
    print(disconnectors)
    print('--')
    extracted_graph = stko.TopologyExtractor()
    tg_info = extracted_graph.extract_topology(
        molecule=cage1,
        broken_bonds_by_id=broken_bonds_by_id,
        disconnectors=set(disconnectors),
    )
    print(tg_info.get_vertex_centroids())
    print(tg_info.get_connectivities())
    print(tg_info.get_edge_pairs())
    cage1.write('output_directory/tg_cage.mol')
    tg_info.write('output_directory/tg_info.pdb')

    new_toplogy_graph = stko.GenericTopologyGraph(
        topology_info=tg_info,
        factory=stko.TopologyFactory(
            topology_graph=stk.cage.Cage,
            vertices=(stk.Vertex(), ),
            edges=(stk.Edge(), ),
        )
    )
    print(new_toplogy_graph.get_topology_graph())


if __name__ == "__main__":
    main()
