# -*- coding: utf-8 -*-

"""Constants for CausalR."""

import enum

from pybel.constants import (
    ANALOGOUS_TO, ASSOCIATION, BIOMARKER_FOR, CAUSES_NO_CHANGE, DECREASES, DIRECTLY_DECREASES, DIRECTLY_INCREASES,
    EQUIVALENT_TO, HAS_COMPONENT, HAS_MEMBER, HAS_PRODUCT, HAS_REACTANT, HAS_VARIANT, INCREASES, IS_A,
    NEGATIVE_CORRELATION, POSITIVE_CORRELATION, PROGONSTIC_BIOMARKER_FOR, RATE_LIMITING_STEP_OF, REGULATES,
    SUBPROCESS_OF, TRANSCRIBED_TO, TRANSLATED_TO,
)

__all__ = [
    'causal_effect_dict',
    'default_edge_ranking',
    'Effect',
]

causal_effect_dict = {
    INCREASES: 1,
    DIRECTLY_INCREASES: 1,
    DECREASES: -1,
    DIRECTLY_DECREASES: -1,
    NEGATIVE_CORRELATION: -1,
    POSITIVE_CORRELATION: 1,
}

default_edge_ranking = {
    INCREASES: 3,
    DIRECTLY_INCREASES: 4,
    DECREASES: 3,
    DIRECTLY_DECREASES: 4,
    RATE_LIMITING_STEP_OF: 0,
    CAUSES_NO_CHANGE: 0,
    REGULATES: 0,
    NEGATIVE_CORRELATION: 2,
    POSITIVE_CORRELATION: 2,
    ASSOCIATION: 1,
    HAS_MEMBER: 0,
    HAS_PRODUCT: 0,
    HAS_COMPONENT: 0,
    HAS_VARIANT: 0,
    HAS_REACTANT: 0,
    TRANSLATED_TO: 0,
    TRANSCRIBED_TO: 0,
    IS_A: 0,
    SUBPROCESS_OF: 0,
    ANALOGOUS_TO: 0,
    BIOMARKER_FOR: 0,
    PROGONSTIC_BIOMARKER_FOR: 0,
    EQUIVALENT_TO: 0
}


class Effect(enum.Enum):
    """Represents the possible effect of a root node on a target."""

    inhibition = -1
    no_effect = 0
    activation = 1
    ambiguous = None
