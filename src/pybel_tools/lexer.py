# -*- coding: utf-8 -*-

from pygments.lexer import RegexLexer, words
from pygments.token import Name, Comment, Text, String, Keyword, Operator, Punctuation

__all__ = ['BELLexer']


class BELLexer(RegexLexer):
    """
    For BEL markup. Supported keywords from open BEL specification 1.0 and 2.0

    :Version: 1
    :Author: Nanditha Mallesh

    """

    name = 'BEL'
    aliases = ['bel']
    filenames = ['*.bel']

    # flags = re.DOTALL

    _keywords_values = (
        'a', 'abundance', 'act', 'analogous', 'association', 'biologicalProcess', 'biomarkerFor', 'bp', 'cat',
        'catalyticActivity', 'causesNoChange', 'cellSecretion', 'cellSurfaceExpression', 'chap', 'chaperoneActivity',
        'complex', 'complexAbundance', 'composite', 'compositeAbundance', 'decreases', 'deg', 'degradation',
        'directlyDecreases', 'directlyIncreases', 'fus', 'fusion', 'g', 'geneAbundance', 'g', 'p', 'gtpBoundActivity',
        'hasComponent', 'hasComponents', 'hasMember', 'hasMembers', 'increases', 'isA', 'kin', 'kinaseActivity', 'list',
        'm', 'microRNAAbundance', 'molecularActivity', 'negativeCorrelation', 'orthologous', 'p', 'path', 'pathology',
        'pep', 'peptidaseActivity', 'phos', 'phosphataseActivity', 'pmod', 'positiveCorrelation', 'products',
        'prognosticBiomarkerFor', 'proteinAbundance', 'proteinModification', 'r', 'rateLimitingStepOf',
        'reactants', 'reaction', 'ribo', 'ribosylationActivity', 'rnaAbundance', 'rxn', 'sec', 'sub', 'subProcessOf',
        'substitution', 'surf', 'tloc', 'tprot', 'transcriedTo', 'transcriptionalActivity',
        'translatedTo', 'translocation', 'transportActivity', 'trunc', 'truncation', 'tscript', 'variant', 'var',
        'fragment', 'frag' 'location', 'loc', 'molecularActivity', 'ma', 'regulates', 'activity', 'cnc', 'neg', 'pos'
    )
    tokens = {
        'root': [
            (words(('SET', 'AS URL', 'STATEMENT_GROUP', 'DOCUMENT', 'DEFINE', 'NAMSPACE', 'UNSET'), suffix=r'\b'),
             Keyword.Constant),
            (words(_keywords_values, suffix=r'\b'), Name.Builtin),
            (r'\w+', Name),
            (r'--', Name.Builtin),
            (r'->', Name.Builtin),
            (r'-\|', Name.Builtin),
            (r'=>', Name.Builtin),
            (r'=\|', Name.Builtin),
            (r'=', Operator),
            (r'\A#!.+$', Comment.Hashbang),
            (r'#.*$', Comment.Single),
            (r'\s+', Text),
            ('"[\S\s]*?\s*"', String),
            ("'[\S\s]*?'\s*", String),
            ('[(){},:]', Punctuation),

        ],

    }
