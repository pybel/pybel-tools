# -*- coding: utf-8 -*-

"""Utilities for displaying graphs with inline HTML in Jupyter Notebooks"""

from random import sample

from pybel.io import to_jsons
from .utils import render_template, default_color_map
from ..mutation import add_canonical_names

__all__ = [
    'to_jupyter',
    'to_jupyter_str'
]

DEFAULT_WIDTH = 1000
DEFAULT_HEIGHT = 650


def generate_id():
    """Generates a random string of letters"""
    return "".join(sample('abcdefghjkmopqrstuvqxyz', 16))


def to_jupyter(graph, width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT, color_map=None, replace_cnames=False):
    """Displays the BEL graph inline in a Jupyter notebook.

    To use successfully, make run as the last statement in a cell inside a Jupyter notebook.

    :param pybel.BELGraph graph: A BEL graph
    :param int width: The width of the visualization window to render
    :param int height: The height of the visualization window to render
    :param dict color_map: A dictionary from PyBEL internal node functions to CSS color strings like #FFEE00. Defaults
                    to :data:`default_color_map`
    :return: An IPython notebook Javascript object
    :rtype: :class:`IPython.display.Javascript`
    """
    from IPython.display import Javascript

    return Javascript(to_jupyter_str(
        graph,
        width=width,
        height=height,
        color_map=color_map,
        replace_cnames=replace_cnames
    ))


def to_jupyter_str(graph, width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT, color_map=None, replace_cnames=False):
    """Returns the string to be javascript-ified by the Jupyter notebook function :class:`IPython.display.Javascript`

    :param pybel.BELGraph graph: A BEL graph
    :param int width: The width of the visualization window to render
    :param int height: The height of the visualization window to render
    :param dict color_map: A dictionary from PyBEL internal node functions to CSS color strings like #FFEE00. Defaults
                    to :data:`default_color_map`
    :return: The javascript string to turn into magic
    :rtype: str
    """
    add_canonical_names(graph, replace=replace_cnames)
    gjson = to_jsons(graph)

    d3_code = render_template('pybel_vis.js')
    chart_id = generate_id()

    color_map = default_color_map if color_map is None else color_map

    javascript_vars = """
        var chart = "{}";
        var width = {};
        var height = {};
        var graph = {};
        const color_map = {};
    """.format(chart_id, width, height, gjson, color_map)

    require_code = """
        require.config({
          paths: {
              d3: '//cdnjs.cloudflare.com/ajax/libs/d3/4.5.0/d3.min'
          }
        });

        var elementInnerHTML = "<div id='" + chart + "'></div>";

        element.append(elementInnerHTML);

        var chartQualified = "#" + chart;

        require(['d3'], function(d3) {
            return init_d3_force(d3, graph, chartQualified, width, height, color_map);
        });
    """

    result = d3_code + javascript_vars + require_code

    return result
