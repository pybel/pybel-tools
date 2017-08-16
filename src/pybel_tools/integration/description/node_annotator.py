# -*- coding: utf-8 -*-

import abc

import six

from ...filters import get_nodes_by_namespace

__all__ = [
    'NodeAnnotator'
]


@six.add_metaclass(abc.ABCMeta)
class NodeAnnotator(object):
    """The base class for node annotators."""

    def __init__(self, namespace):
        """
        :param str or list[str] namespace: The name of the namespace or namespaces that this node annotator services
        """
        self.namespace = namespace
        self.download_successful = False

    @abc.abstractmethod
    def get_description(self, name):
        """Gets the description for the given name in this annotator's namespace."""

    def get_label(self, name):
        """Gets the label for the given name in. If not overridden, uses each node's name as its label."""

    def populate_by_graph(self, graph):
        """Optional hook for populating the annotator based on the nodes in a graph. Override this if
        your node annotator downloads data ahead of time, such as grouping requests to an external API.

        :param pybel.BELGraph graph: A BEL graph
        """

    def annotate(self, graph):
        """Annotates all nodes in this annotator's namespace

        :param pybel.BELGraph graph: A BEL graph
        """
        self.populate_by_graph(graph)

        for node in get_nodes_by_namespace(graph, self.namespace):
            name = graph.get_node_name(node)

            description = self.get_description(name)
            graph.set_node_description(node, description)

            label = self.get_label(name)
            if label:
                graph.set_node_label(node, label)
