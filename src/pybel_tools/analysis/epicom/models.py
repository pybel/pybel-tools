# -*- coding: utf-8 -*-

from sqlalchemy import Column, Float, ForeignKey, Integer
from sqlalchemy.orm import relationship

from pybel.manager.models import (
    ANNOTATION_ENTRY_TABLE_NAME, AnnotationEntry, Base, NAMESPACE_ENTRY_TABLE_NAME,
    NETWORK_TABLE_NAME, NamespaceEntry, Network,
)

TABLE_BASE_NAME = 'epicom'
SCORE_TABLE_NAME = 'epicom_score'


class Score(Base):
    """Represents the score for a subgraph in a given network"""
    __tablename__ = SCORE_TABLE_NAME

    id = Column(Integer, primary_key=True)

    network_id = Column(Integer, ForeignKey('{}.id'.format(NETWORK_TABLE_NAME)))
    network = relationship(Network)

    annotation_id = Column(Integer, ForeignKey('{}.id'.format(ANNOTATION_ENTRY_TABLE_NAME)))
    annotation = relationship(AnnotationEntry)

    drug_id = Column(Integer, ForeignKey('{}.id'.format(NAMESPACE_ENTRY_TABLE_NAME)))
    drug = relationship(NamespaceEntry)

    score = Column(Float, nullable=False, unique=False, index=True)
