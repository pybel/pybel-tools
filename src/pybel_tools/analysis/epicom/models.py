# -*- coding: utf-8 -*-

from sqlalchemy import Column, Float, ForeignKey, Integer
from sqlalchemy.orm import relationship

from pybel.manager.models import (
    Base, NAMESPACE_TABLE_NAME, NAME_TABLE_NAME, NETWORK_TABLE_NAME, NamespaceEntry,
    Network,
)

TABLE_BASE_NAME = 'epicom'
SCORE_TABLE_NAME = 'epicom_score'


class Score(Base):
    """Represents the score for a sub-graph in a given network"""
    __tablename__ = SCORE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    network_id = Column(Integer, ForeignKey('{}.id'.format(NETWORK_TABLE_NAME)))
    network = relationship(Network)

    annotation_id = Column(Integer, ForeignKey('{}.id'.format(NAMESPACE_TABLE_NAME)))
    annotation = relationship(NamespaceEntry)

    drug_id = Column(Integer, ForeignKey('{}.id'.format(NAME_TABLE_NAME)))
    drug = relationship(NamespaceEntry)

    score = Column(Float, nullable=False, unique=False, index=True)
