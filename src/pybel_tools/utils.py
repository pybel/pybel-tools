# -*- coding: utf-8 -*-

"""This module contains functions useful throughout PyBEL Tools."""

import datetime
import itertools as itt
import json
import logging
import typing
from collections import Counter, defaultdict
from operator import itemgetter
from typing import Callable, Dict, Iterable, List, Mapping, Optional, Set, Sized, Tuple, TypeVar, Union

import networkx as nx

from pybel import BELGraph

log = logging.getLogger(__name__)

CENTRALITY_SAMPLES = 200
X = TypeVar('X')
Y = TypeVar('Y')


def pairwise(iterable: Iterable[X]) -> Iterable[Tuple[X, X]]:
    """Iterate over pairs in list.

    s -> (s0,s1), (s1,s2), (s2, s3), ...
    """
    a, b = itt.tee(iterable)
    next(b, None)
    return zip(a, b)


def group_as_sets(pairs: Iterable[Tuple[X, Y]]) -> Mapping[X, Set[Y]]:
    """Group elements in the iterable in a dictionary of sets."""
    rv = defaultdict(set)
    for x, y in pairs:
        rv[x].add(y)
    return dict(rv)


def group_as_lists(pairs: Iterable[Tuple[X, Y]]) -> Mapping[X, List[Y]]:
    """Group elements in the iterable in a dictionary of lists."""
    rv = defaultdict(list)
    for x, y in pairs:
        rv[x].append(y)
    return dict(rv)


def count_defaultdict(dict_of_lists: Mapping[X, List[Y]]) -> Mapping[X, typing.Counter[Y]]:
    """Count the number of elements in each value of the dictionary."""
    return {
        k: Counter(v)
        for k, v in dict_of_lists.items()
    }


def count_dict_values(dict_of_counters: Mapping[X, Sized]) -> typing.Counter[X]:
    """Count the number of elements in each value (can be list, Counter, etc).

    :param dict_of_counters: A dictionary of things whose lengths can be measured (lists, Counters, dicts)
    :return: A Counter with the same keys as the input but the count of the length of the values list/tuple/set/Counter
    """
    return Counter({
        k: len(v)
        for k, v in dict_of_counters.items()
    })


def set_percentage(x: Iterable[X], y: Iterable[X]) -> float:
    """Return what percentage of x is contained within y.

    :param x: A set
    :param y: Another set
    :return: The percentage of x contained within y
    """
    a, b = set(x), set(y)

    if not a:
        return 0.0

    return len(a & b) / len(a)


def tanimoto_set_similarity(x: Iterable[X], y: Iterable[X]) -> float:
    """Calculate the tanimoto set similarity."""
    a, b = set(x), set(y)
    union = a | b

    if not union:
        return 0.0

    return len(a & b) / len(union)


def min_tanimoto_set_similarity(x: Iterable[X], y: Iterable[X]) -> float:
    """Calculate the tanimoto set similarity using the minimum size.

    :param x: A set
    :param y: Another set
    :return: The similarity between
    """
    a, b = set(x), set(y)

    if not a or not b:
        return 0.0

    return len(a & b) / min(len(a), len(b))


def calculate_single_tanimoto_set_distances(
    target: Iterable[X],
    dict_of_sets: Mapping[Y, Set[X]],
) -> Mapping[Y, float]:
    """Return a dictionary of distances keyed by the keys in the given dict.

    Distances are calculated based on pairwise tanimoto similarity of the sets contained

    :param target: A set
    :param dict_of_sets: A dict of {x: set of y}
    :return: A similarity dictionary based on the set overlap (tanimoto) score between the target set and the sets in
            dos
    """
    target_set = set(target)

    return {
        k: tanimoto_set_similarity(target_set, s)
        for k, s in dict_of_sets.items()
    }


def calculate_tanimoto_set_distances(
    dict_of_sets: Mapping[X, Set],
) -> Mapping[X, Mapping[X, float]]:
    """Return a distance matrix keyed by the keys in the given dict.

    Distances are calculated based on pairwise tanimoto similarity of the sets contained.

    :param dict_of_sets: A dict of {x: set of y}
    :return: A similarity matrix based on the set overlap (tanimoto) score between each x as a dict of dicts
    """
    result: Dict[X, Dict[X, float]] = defaultdict(dict)

    for x, y in itt.combinations(dict_of_sets, 2):
        result[x][y] = result[y][x] = tanimoto_set_similarity(dict_of_sets[x], dict_of_sets[y])

    for x in dict_of_sets:
        result[x][x] = 1.0

    return dict(result)


def calculate_global_tanimoto_set_distances(dict_of_sets: Mapping[X, Set]) -> Mapping[X, Mapping[X, float]]:
    r"""Calculate an alternative distance matrix based on the following equation.

    .. math:: distance(A, B)=1- \|A \cup B\| / \| \cup_{s \in S} s\|

    :param dict_of_sets: A dict of {x: set of y}
    :return: A similarity matrix based on the alternative tanimoto distance as a dict of dicts
    """
    universe = set(itt.chain.from_iterable(dict_of_sets.values()))
    universe_size = len(universe)

    result: Dict[X, Dict[X, float]] = defaultdict(dict)

    for x, y in itt.combinations(dict_of_sets, 2):
        result[x][y] = result[y][x] = 1.0 - len(dict_of_sets[x] | dict_of_sets[y]) / universe_size

    for x in dict_of_sets:
        result[x][x] = 1.0 - len(x) / universe_size

    return dict(result)


def barh(d, plt, title=None):
    """Plot a horizontal bar plot from a Counter."""
    labels = sorted(d, key=d.get)
    index = range(len(labels))

    plt.yticks(index, labels)
    plt.barh(index, [d[v] for v in labels])

    if title is not None:
        plt.title(title)


def barv(d, plt, title=None, rotation='vertical'):
    """Plot a vertical bar plot from a Counter."""
    labels = sorted(d, key=d.get, reverse=True)
    index = range(len(labels))
    plt.xticks(index, labels, rotation=rotation)
    plt.bar(index, [d[v] for v in labels])

    if title is not None:
        plt.title(title)


def prepare_c3(
    data: Union[List[Tuple[str, int]], Mapping[str, int]],
    y_axis_label: str = 'y',
    x_axis_label: str = 'x',
) -> str:
    """Prepare a C3 JSON for making a bar chart from a Counter.

    :param data: A dictionary of {str: int} to display as bar chart
    :param y_axis_label: The Y axis label
    :param x_axis_label: X axis internal label. Should be left as default 'x')
    :return: A JSON dictionary for making a C3 bar chart
    """
    if not isinstance(data, list):
        data = sorted(data.items(), key=itemgetter(1), reverse=True)

    try:
        labels, values = zip(*data)
    except ValueError:
        log.info(f'no values found for {x_axis_label}, {y_axis_label}')
        labels, values = [], []

    return json.dumps([
        [x_axis_label] + list(labels),
        [y_axis_label] + list(values),
    ])


def prepare_c3_time_series(data: List[Tuple[int, int]], y_axis_label: str = 'y', x_axis_label: str = 'x') -> str:
    """Prepare a C3 JSON string dump for a time series.

    :param data: A list of tuples [(year, count)]
    :param y_axis_label: The Y axis label
    :param x_axis_label: X axis internal label. Should be left as default 'x')
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


def calculate_betweenness_centality(graph: BELGraph, number_samples: int = CENTRALITY_SAMPLES) -> Counter:
    """Calculate the betweenness centrality over nodes in the graph.

    Tries to do it with a certain number of samples, but then tries a complete approach if it fails.
    """
    try:
        res = nx.betweenness_centrality(graph, k=number_samples)
    except Exception:
        res = nx.betweenness_centrality(graph)
    return Counter(res)


T = TypeVar('T', List, Tuple)


def get_circulations(elements: T) -> Iterable[T]:
    """Iterate over all possible circulations of an ordered collection (tuple or list).

    Example:
    >>> list(get_circulations([1, 2, 3]))
    [[1, 2, 3], [2, 3, 1], [3, 1, 2]]

    """
    for i in range(len(elements)):
        yield elements[i:] + elements[:i]


def canonical_circulation(elements: T, key: Optional[Callable[[T], bool]] = None) -> T:
    """Get the canonical representation of the ordered collection.

    Does this by finding its minimum circulation with the given sort key
    """
    return min(get_circulations(elements), key=key)
