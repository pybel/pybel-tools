# -*- coding: utf-8 -*-

"""Generate HTML summaries of BEL graphs."""

from __future__ import annotations

import logging
from typing import Optional, TextIO, Tuple

from pybel import BELGraph
from pybel.dsl import BaseEntity
from ..jinja_utils import build_template_renderer
from ...summary import BELGraphSummary
from ...utils import prepare_c3, prepare_c3_time_series

__all__ = [
    'to_html',
    'to_html_file',
    'to_html_path',
    'prepare_c3',
]

log = logging.getLogger(__name__)

CONFIDENCES = 'None', 'Very Low', 'Low', 'Medium', 'High', 'Very High'

render_template = build_template_renderer(__file__)
render_template.environment.globals.update(
    prepare_c3=prepare_c3,
    prepare_c3_time_series=prepare_c3_time_series,
)


def to_html_path(graph: BELGraph, path: str) -> None:
    """Write the graph as HTML to a file at the given path."""
    with open(path, 'w') as file:
        to_html_file(graph, file)


def to_html_file(graph: BELGraph, file: Optional[TextIO] = None) -> None:
    """Write the graph as HTML to a file."""
    html = to_html(graph)
    print(html, file=file)


def to_html(graph: BELGraph) -> str:
    """Render the graph as an HTML string.

    Common usage may involve writing to a file like:

    >>> from pybel.examples import sialic_acid_graph
    >>> with open('html_output.html', 'w') as file:
    ...     print(to_html(sialic_acid_graph), file=file)
    """
    summary = BELGraphSummary.from_graph(graph)

    confidence_data = [
        (label, summary.confidence_count.get(label, 0))
        for label in CONFIDENCES
    ]

    return render_template(
        'index.html',
        graph=graph,
        summary=summary,
        confidence_data=confidence_data,
    )


# def get_network_summary_dict(graph: BELGraph) -> Mapping:
#     """Create a summary dictionary."""
#     summary = BELGraphSummary.from_graph(graph)
#
#     contradictory_triplets = list(itt.chain(
#         (get_triplet_tuple(a, b, c) + ('Separate',) for a, b, c in summary.separate_unstable_correlation_triples),
#         (get_triplet_tuple(a, b, c) + ('Mutual',) for a, b, c in summary.mutually_unstable_correlation_triples),
#         (get_triplet_tuple(a, b, c) + ('Jens',) for a, b, c in summary.jens_unstable),
#         (get_triplet_tuple(a, b, c) + ('Increase Mismatch',) for a, b, c in summary.increase_mismatch_triplets),
#         (get_triplet_tuple(a, b, c) + ('Decrease Mismatch',) for a, b, c in summary.decrease_mismatch_triplets),
#     )),
#
#     regulatory_pairs = [
#         get_pair_tuple(u, v)
#         for u, v in summary.regulatory_pairs
#     ]
#
#     unstable_pairs = list(itt.chain(
#         (get_pair_tuple(u, v) + ('Chaotic',) for u, v in summary.chaotic_pairs),
#         (get_pair_tuple(u, v) + ('Dampened',) for u, v in get_dampened_pairs(graph)),
#     ))
#
#     contradictory_pairs = [
#         get_pair_tuple(u, v) + (relation,)
#         for u, v, relation in summary.contradictory_pairs
#     ]
#
#     rv.update(summary.to_dict())
#
#     return rv


BELPairTuple = Tuple[str, str, str, str]
BELTripleTuple = Tuple[str, str, str, str, str, str]


def get_pair_tuple(a: BaseEntity, b: BaseEntity) -> BELPairTuple:
    """Get the pair as a tuple of BEL/hashes."""
    return a.as_bel(), a.sha512, b.as_bel(), b.sha512


def get_triplet_tuple(a: BaseEntity, b: BaseEntity, c: BaseEntity) -> BELTripleTuple:
    """Get the triple as a tuple of BEL/hashes."""
    return a.as_bel(), a.sha512, b.as_bel(), b.sha512, c.as_bel(), c.sha512
