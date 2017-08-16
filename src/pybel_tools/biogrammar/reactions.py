# -*- coding: utf-8 -*-

"""

This module checks that the substrates and products of reactions are feasible/correct

1. Convert all abundances from CHEBI to modern identifier, like SMILES or inchi-key
2. Download reaction data from (source??)

"""

#: ID, EC Number, Reccomended Name, Reaction, Rxn ID Brenda, Rxn ID KEGG,
#: Rxn ID MetaCyc, Rxn ID SABIO-RK, Brenda Pathway Name, KEGG Pathway ID
#: Kegg pathway name, MetaCyc pathway ID, MetaCyc Pathway Name, Stoichiometry Check
#: Missing Substrate, Missing Produt, Commentary KEGG, Commentary MetaCyc, Remark
BKMS_REACT_DATA = 'http://bkm-react.tu-bs.de/download/Reactions_BKMS.tar.gz'
