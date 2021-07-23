"""
Topology Graph
==============

#. :class:`.GenericTopologyGraph`

Class for defining a topology from a molecule and disconnections.

"""


class GenericTopologyGraph:
    """
    Class containing stk.TopologyGraph classes.

    Examples
    --------

    .. code-block:: python

    """

    def __init__(self, topology_info, factory):
        """
        Initialize a :class:`TopologyFactory`.

        """

        self._define_graph(topology_info, factory)

    def _define_graph(self, topology_info, factory):
        self._topology_graph = factory.get_topology_graph(topology_info)

    def get_topology_graph(self):
        return self._topology_graph
