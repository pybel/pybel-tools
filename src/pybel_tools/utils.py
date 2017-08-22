# -*- coding: utf-8 -*-

"""This module contains functions useful throughout PyBEL Tools"""

import datetime
import itertools as itt
import json
import logging
import os
import time
from collections import Counter, defaultdict
from itertools import zip_longest
from operator import itemgetter

import jinja2
import networkx as nx
import pandas as pd
from pkg_resources import get_distribution
from pybel.constants import (
    ANNOTATIONS,
    CITATION_TYPE,
    CITATION_NAME,
    CITATION_REFERENCE,
    CITATION_DATE,
    CITATION_AUTHORS,
    CITATION_COMMENTS,
    RELATION,
)

log = logging.getLogger(__name__)

CENTRALITY_SAMPLES = 200


def multidict_list(it):
    result = defaultdict(list)
    for k, v in it:
        result[k].append(v)
    return dict(result)


def multidict_set(it):
    result = defaultdict(set)
    for k, v in it:
        result[k].add(v)
    return dict(result)


def pairwise(iterable):
    """ Iterate over pairs in list s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = itt.tee(iterable)
    next(b, None)
    return zip(a, b)


def graph_edge_data_iter(graph):
    """Iterates over the edge data dictionaries

    :param graph: A BEL graph
    :type graph: pybel.BELGraph
    :return: An iterator over the edge dictionaries in the graph
    :rtype: iter
    """
    for _, _, data in graph.edges_iter(data=True):
        yield data


def count_defaultdict(dict_of_lists):
    """Takes a dictionary and applies a counter to each list

    :param dict_of_lists: A dictionary of lists
    :type dict_of_lists: dict or collections.defaultdict
    :return: A dictionary of {key: Counter(values)}
    :rtype: dict
    """
    return {
        k: Counter(v)
        for k, v in dict_of_lists.items()
    }


def get_value_sets(dict_of_iterables):
    """Takes a dictionary of lists/iterables/counters and gets the sets of the values

    :param dict_of_iterables: A dictionary of lists
    :type dict_of_iterables: dict or collections.defaultdict
    :return: A dictionary of {key: set of values}
    :rtype: dict
    """
    return {
        k: set(v)
        for k, v in dict_of_iterables.items()
    }


def count_dict_values(dict_of_counters):
    """Counts the number of elements in each value (can be list, Counter, etc)

    :param dict_of_counters: A dictionary of things whose lengths can be measured (lists, Counters, dicts)
    :type dict_of_counters: dict or collections.defaultdict
    :return: A Counter with the same keys as the input but the count of the length of the values list/tuple/set/Counter
    :rtype: collections.Counter
    """
    return Counter({
        k: len(v)
        for k, v in dict_of_counters.items()
    })


def _check_has_data(d, sd, key):
    return sd in d and key in d[sd]


def check_has_annotation(data, key):
    """Checks that ANNOTATION is included in the data dictionary and that the key is also present

    :param data: The data dictionary from a BELGraph's edge
    :type data: dict
    :param key: An annotation key
    :param key: str
    :return: If the annotation key is present in the current data dictionary
    :rtype: bool

    For example, it might be useful to print all edges that are annotated with 'Subgraph':

    >>> from pybel import BELGraph
    >>> graph = BELGraph()
    >>> ...
    >>> for u, v, data in graph.edges_iter(data=True):
    >>>     if not check_has_annotation(data, 'Subgraph')
    >>>         continue
    >>>     print(u, v, data)
    """
    return _check_has_data(data, ANNOTATIONS, key)


def set_percentage(x, y):
    """What percentage of x is contained within y?

    :param set x: A set
    :param set y: Another set
    :return: The percentage of x contained within y
    :rtype: float
    """
    a, b = set(x), set(y)

    if not a:
        return 0.0

    return len(a & b) / len(a)


def tanimoto_set_similarity(x, y):
    """Calculates the tanimoto set similarity

    :param set x: A set
    :param set y: Another set
    :return: The similarity between
    :rtype: float
    """
    a, b = set(x), set(y)
    union = a | b

    if not union:
        return 0.0

    return len(a & b) / len(union)


def min_tanimoto_set_similarity(x, y):
    """Calculates the tanimoto set similarity using the minimum size

    :param set x: A set
    :param set y: Another set
    :return: The similarity between
    :rtype: float
        """
    a, b = set(x), set(y)

    if not a or not b:
        return 0.0

    return len(a & b) / min(len(a), len(b))


def calculate_single_tanimoto_set_distances(target, dict_of_sets):
    """Returns a dictionary of distances keyed by the keys in the given dict. Distances are calculated
    based on pairwise tanimoto similarity of the sets contained

    :param set target: A set
    :param dict_of_sets: A dict of {x: set of y}
    :type dict_of_sets: dict
    :return: A similarity dicationary based on the set overlap (tanimoto) score between the target set and the sets in
            dos
    :rtype: dict
    """
    target_set = set(target)

    return {
        k: tanimoto_set_similarity(target_set, s)
        for k, s in dict_of_sets.items()
    }


def calculate_tanimoto_set_distances(dict_of_sets):
    """Returns a distance matrix keyed by the keys in the given dict. Distances are calculated
    based on pairwise tanimoto similarity of the sets contained

    :param dict_of_sets: A dict of {x: set of y}
    :type dict_of_sets: dict
    :return: A similarity matrix based on the set overlap (tanimoto) score between each x as a dict of dicts
    :rtype: dict
    """
    result = defaultdict(dict)

    for x, y in itt.combinations(dict_of_sets, 2):
        result[x][y] = result[y][x] = tanimoto_set_similarity(dict_of_sets[x], dict_of_sets[y])

    for x in dict_of_sets:
        result[x][x] = 1.0

    return dict(result)


def calculate_global_tanimoto_set_distances(dict_of_sets):
    """Calculates an alternative distance matrix based on the following equation:

    .. math:: distance(A, B)=1- \|A \cup B\| / \| \cup_{s \in S} s\|

    :param dict_of_sets: A dict of {x: set of y}
    :type dict_of_sets: dict
    :return: A similarity matrix based on the alternative tanimoto distance as a dict of dicts
    :rtype: dict
    """
    universe = set(itt.chain.from_iterable(dict_of_sets.values()))
    universe_size = len(universe)

    result = defaultdict(dict)

    for x, y in itt.combinations(dict_of_sets, 2):
        result[x][y] = result[y][x] = 1.0 - len(dict_of_sets[x] | dict_of_sets[y]) / universe_size

    for x in dict_of_sets:
        result[x][x] = 1.0 - len(x) / universe_size

    return dict(result)


def all_edges_iter(graph, u, v):
    """Lists all edges between the given nodes

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: A BEL node
    :param tuple v: A BEL node
    :return: A list of (node, node, key)
    :rtype: list[tuple]
    """
    if u not in graph.edge or v not in graph.edge[u]:
        raise ValueError('Graph has no edges')

    for k in graph.edge[u][v].keys():
        yield u, v, k


def barh(d, plt, title=None):
    """A convenience function for plotting a horizontal bar plot from a Counter"""
    labels = sorted(d, key=d.get)
    index = range(len(labels))

    plt.yticks(index, labels)
    plt.barh(index, [d[v] for v in labels])

    if title is not None:
        plt.title(title)


def barv(d, plt, title=None, rotation='vertical'):
    """A convenience function for plotting a vertical bar plot from a Counter"""
    labels = sorted(d, key=d.get, reverse=True)
    index = range(len(labels))
    plt.xticks(index, labels, rotation=rotation)
    plt.bar(index, [d[v] for v in labels])

    if title is not None:
        plt.title(title)


def citation_to_tuple(citation):
    """Converts a citation dictionary to a tuple. Can be useful for sorting and serialization purposes

    :param dict citation: A citation dictionary
    :return: A citation tuple
    :rtype: tuple
    """
    return tuple([
        citation.get(CITATION_TYPE),
        citation.get(CITATION_NAME),
        citation.get(CITATION_REFERENCE),
        citation.get(CITATION_DATE),
        citation.get(CITATION_AUTHORS),
        citation.get(CITATION_COMMENTS)
    ])


def is_edge_consistent(graph, u, v):
    """Checks if all edges between two nodes have the same relation

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: The source BEL node
    :param tuple v: The target BEL node
    :return: If all edges from the source to target node have the same relation
    :rtype: bool
    """
    if not graph.has_edge(u, v):
        raise ValueError('{} does not contain an edge ({}, {})'.format(graph, u, v))

    return 0 == len(set(d[RELATION] for d in graph.edge[u][v].values()))


def safe_add_edge(graph, u, v, key, attr_dict, **attr):
    """Adds an edge while preserving negative keys, and paying no respect to positive ones

    :param pybel.BELGraph graph: A BEL Graph
    :param tuple u: The source BEL node
    :param tuple v: The target BEL node
    :param int key: The edge key. If less than zero, corresponds to an unqualified edge, else is disregarded
    :param dict attr_dict: The edge data dictionary
    :param dict attr: Edge data to assign via keyword arguments
    """
    if key < 0:
        graph.add_edge(u, v, key=key, attr_dict=attr_dict, **attr)
    else:
        graph.add_edge(u, v, attr_dict=attr_dict, **attr)


def safe_add_edges(graph, edges):
    """Adds an iterable of edges to the graph

    :param pybel.BELGraph graph: A BEL Graph
    :param iter[tuple edges: An iterable of 4-tuples of (source, target, key, data)
    """
    for source, target, key, attr_dict in edges:
        safe_add_edge(graph, source, target, key=key, attr_dict=attr_dict)


def load_differential_gene_expression(data_path, gene_symbol_column='Gene.symbol', logfc_column='logFC'):
    """Quick and dirty loader for differential gene expression data

    :param str data_path:
    :param str gene_symbol_column:
    :param str logfc_colun:
    :return: A dictionary of {gene symbol: log fold change}
    :rtype: dict
    """
    df = pd.read_csv(data_path)
    df = df.loc[df[gene_symbol_column].notnull(), [gene_symbol_column, logfc_column]]

    return {
        k: v
        for _, k, v in df.itertuples()
    }


def prepare_c3(data, y_axis_label='y', x_axis_label='x'):
    """Prepares C3 JSON for making a bar chart from a Counter

    :param data: A dictionary of {str: int} to display as bar chart
    :type data: Counter or dict or collections.defaultdict
    :param str y_axis_label: The Y axis label
    :param str x_axis_label: X axis internal label. Should be left as default 'x')
    :return: A JSON dictionary for making a C3 bar chart
    :rtype: dict
    """
    labels, values = [], []
    for k, v in sorted(data.items(), key=itemgetter(1), reverse=True):
        labels.append(k)
        values.append(v)

    return json.dumps([
        [x_axis_label] + list(labels),
        [y_axis_label] + list(values)
    ])


def prepare_c3_time_series(data, y_axis_label='y', x_axis_label='x'):
    """Prepares C3 JSON for making a time series

    :param data: A list of tuples [(year, count)]
    :type data: list
    :param str y_axis_label: The Y axis label
    :param str x_axis_label: X axis internal label. Should be left as default 'x')
    :return: A JSON dictionary for making a C3 bar chart
    :rtype: dict
    """
    years, counter = zip(*data)

    years = [
        datetime.date(year, 1, 1).isoformat()
        for year in years
    ]

    return json.dumps([
        [x_axis_label] + list(years),
        [y_axis_label] + list(counter)
    ])


def get_version():
    """Gets the current PyBEL Tools version

    :return: The current PyBEL Tools version
    :rtype: str
    """
    return get_distribution('pybel_tools').version


def build_template_environment(here):
    """Builds a custom templating enviroment so Flask apps can get data from lots of different places
    
    :param str here: Give this the result of :code:`os.path.dirname(os.path.abspath(__file__))`
    :rtype: jinja2.Environment
    """
    template_environment = jinja2.Environment(
        autoescape=False,
        loader=jinja2.FileSystemLoader(os.path.join(here, 'templates')),
        trim_blocks=False
    )

    template_environment.globals['STATIC_PREFIX'] = here + '/static/'

    return template_environment


def build_template_renderer(file):
    """In your file, give this function the current file

    :param str file: The location of the current file. Pass it :code:`__file__` like in the example below.
    
    >>> render_template = build_template_renderer(__file__)
    """
    here = os.path.dirname(os.path.abspath(file))
    template_environment = build_template_environment(here)

    def render_template(template_filename, **context):
        """Renders a template as a unicode string

        :param str template_filename: The name of the file to render in the template directory
        :param dict context: The variables to template
        :rtype: str
        """
        return template_environment.get_template(template_filename).render(context)

    return render_template


def enable_cool_mode():
    log.info('enabled cool mode')
    logging.getLogger('pybel.parser').setLevel(50)


def calculate_betweenness_centality(graph, k=CENTRALITY_SAMPLES):
    """Calculates the betweenness centrality over nodes in the graph. Tries to do it with a certain number of samples,
    but then tries a complete approach if it fails.

    :param pybel.BELGraph graph: A BEL graph
    :param int k: The number of samples to use
    :rtype: collections.Counter[tuple,float]
    """
    try:
        res = Counter(nx.betweenness_centrality(graph, k=k))
        return res
    except:
        return Counter(nx.betweenness_centrality(graph))


def grouper(n, iterable, fillvalue=None):
    "grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def get_iso_8601_date():
    return time.strftime('%Y%m%d')


def hash_str_to_int(hash_str, length=16):
    """Hashes a tuple to the given number of digits

    :param str hash_str: Basically anything that can be pickled deterministically
    :param int length: The length of the hash to keep
    :rtype: int
    """
    return int(hash_str, 16) % (10 ** length)


def get_circulations(t):
    """Iterate over all possible circulations of an ordered collection (tuple or list)

    :param tuple or list t:
    :rtype: iter
    """
    for i in range(len(t)):
        yield t[i:] + t[:i]


def canonical_circulation(t, key=None):
    """Get get a canonical representation of the ordered collection by finding its minimum circulation with the
    given sort key

    :param tuple or list t:
    :param key: A function for sort
    :return: The
    """
    return min(get_circulations(t), key=key)
