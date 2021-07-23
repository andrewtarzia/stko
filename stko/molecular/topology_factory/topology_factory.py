"""
Topology Factory
================

#. :class:`.TopologyFactory`

Class for defining a topology from a molecule and disconnections.

"""


class TopologyFactory:
    """
    Factory of stko.GenericTopologyGraph classes.

    Examples
    --------

    .. code-block:: python

    """

    def __init__(self, topology_graph, vertices, edges):
        """
        Initialize a :class:`TopologyFactory`.

        """

        self._topology_graph = topology_graph
        self._vertices = vertices
        self._edges = edges

    def get_topology_graph(self, topology_info):
        new_topology_graph = self._topology_graph.clone()

        vertex_prototypes = []
        edge_prototypes = []

        new_topology_graph._vertex_prototypes = vertex_prototypes
        new_topology_graph._edge_prototypes = edge_prototypes

        return new_topology_graph
