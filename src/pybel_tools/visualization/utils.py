# -*- coding: utf-8 -*-

from pybel.constants import (
    PROTEIN,
    PATHOLOGY,
    BIOPROCESS,
    MIRNA,
    COMPLEX,
    COMPOSITE,
    REACTION,
    GENE,
    ABUNDANCE,
    RNA
)
from ..utils import build_template_renderer

__all__ = [
    'render_template',
    'default_color_map'
]

#: Renders templates from pybel_tools.visualization.templates folder
render_template = build_template_renderer(__file__)

#: The color map defining the node colors in visualization
default_color_map = {
    PROTEIN: "#1F77B4",
    PATHOLOGY: "#FF7F0E",
    BIOPROCESS: "#2CA02C",
    MIRNA: "#D62728",
    COMPLEX: "#9467bd",
    COMPOSITE: "#9467bd",
    REACTION: "#8c564b",
    GENE: "#e377c2",
    ABUNDANCE: "#bcbd22",
    RNA: "#17becf"
}
