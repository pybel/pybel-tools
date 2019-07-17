# -*- coding: utf-8 -*-

"""Database models for the results of EpiCom Reloaded."""

from sqlalchemy import Column, Float, ForeignKey, Integer
from sqlalchemy.orm import relationship

from pybel.manager.models import Base, NamespaceEntry, Network

TABLE_BASE_NAME = 'epicom'
SCORE_TABLE_NAME = 'epicom_score'


class Score(Base):
    """Represents the score for a sub-graph in a given network."""

    __tablename__ = SCORE_TABLE_NAME
    id = Column(Integer, primary_key=True)

    network_id = Column(Integer, ForeignKey(f'{Network.__tablename__}.id'))
    network = relationship(Network)

    annotation_id = Column(Integer, ForeignKey(f'{NamespaceEntry.__tablename__}.id'))
    annotation = relationship(NamespaceEntry)

    drug_id = Column(Integer, ForeignKey(f'{NamespaceEntry.__tablename__}.id'))
    drug = relationship(NamespaceEntry)

    score = Column(Float, nullable=False, unique=False, index=True)
