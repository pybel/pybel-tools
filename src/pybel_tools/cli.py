# -*- coding: utf-8 -*-

"""A command line interface for PyBEL-Tools.

Why does this file exist, and why not put this in ``__main__``?
You might be tempted to import things from ``__main__`` later, but that will cause
problems--the code will get executed twice:
 - When you run `python3 -m pybel_tools` python will execute
   ``__main__.py`` as a script. That means there won't be any
   ``pybel_tools.__main__`` in ``sys.modules``.
 - When you import __main__ it will get executed again (as a module) because
   there's no ``pybel_tools.__main__`` in ``sys.modules``.
Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""

import os
import sys
from getpass import getuser
from typing import TextIO

import click

from bel_resources import parse_bel_resource, write_annotation, write_namespace
from pybel import BELGraph, Manager, from_lines
from pybel.cli import connection_option, graph_pickle_argument
from pybel.constants import NAMESPACE_DOMAIN_OTHER
from pybel.struct import get_pubmed_identifiers
from pybel.utils import get_version as pybel_version
from .utils import get_version


@click.group(help=f"PyBEL-Tools v{get_version()} Command Line Interface on {sys.executable}\n"
                  f"with PyBEL v{pybel_version()}")
@click.version_option()
def main():
    """Command Line Interface for PyBEL Tools."""


@main.group()
@connection_option
@click.pass_context
def io(ctx, connection: str):
    """Upload and conversion utilities."""
    ctx.obj = Manager(connection=connection)


@main.group()
def namespace():
    """Namespace file utilities."""


@namespace.command()
@click.argument('name')
@click.argument('keyword')
@click.argument('domain')
@click.argument('citation')
@click.option('--author', default=getuser())
@click.option('--description')
@click.option('--species')
@click.option('--version')
@click.option('--contact')
@click.option('--license')
@click.option('--values', default=sys.stdin, help="A file containing the list of names")
@click.option('--output', type=click.File('w'), default=sys.stdout)
def write(name, keyword, domain, citation, author, description, species, version, contact, license, values, output):
    """Build a namespace from items."""
    write_namespace(
        name, keyword, domain, author, citation, values,
        namespace_description=description,
        namespace_species=species,
        namespace_version=version,
        author_contact=contact,
        author_copyright=license,
        file=output,
    )


@namespace.command()
@click.option('-f', '--file', type=click.File('r'), default=sys.stdin, help="Path to input BEL Namespace file")
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout,
              help="Path to output converted BEL Annotation file")
def convert_to_annotation(file, output):
    """Convert a namespace file to an annotation file."""
    resource = parse_bel_resource(file)

    write_annotation(
        keyword=resource['Namespace']['Keyword'],
        values={k: '' for k in resource['Values']},
        citation_name=resource['Citation']['NameString'],
        description=resource['Namespace']['DescriptionString'],
        file=output,
    )


@main.group()
def annotation():
    """Annotation file utilities."""


@annotation.command()
@click.option('-f', '--file', type=click.File('r'), default=sys.stdin, help="Path to input BEL Namespace file")
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout,
              help="Path to output converted BEL Namespace file")
@click.option('--keyword', help="Set custom keyword. useful if the annotation keyword is too long")
def convert_to_namespace(file, output, keyword):
    """Convert an annotation file to a namespace file."""
    resource = parse_bel_resource(file)

    write_namespace(
        namespace_keyword=(keyword or resource['AnnotationDefinition']['Keyword']),
        namespace_name=resource['AnnotationDefinition']['Keyword'],
        namespace_description=resource['AnnotationDefinition']['DescriptionString'],
        author_name='Charles Tapley Hoyt',
        namespace_domain=NAMESPACE_DOMAIN_OTHER,
        values=resource['Values'],
        citation_name=resource['Citation']['NameString'],
        file=output,
    )


@main.group()
def document():
    """BEL document utilities."""


@document.command()
@click.argument('name')
@click.argument('contact')
@click.argument('description')
@click.argument('pmids', nargs=-1)
@click.option('--version')
@click.option('--copyright')
@click.option('--authors')
@click.option('--licenses')
@click.option('--disclaimer')
@click.option('--output', type=click.File('wb'), default=sys.stdout)
def boilerplate(name, contact, description, pmids, version, copyright, authors, licenses, disclaimer, output):
    """Build a template BEL document with the given PubMed identifiers."""
    from .document_utils import write_boilerplate

    write_boilerplate(
        name=name,
        version=version,
        description=description,
        authors=authors,
        contact=contact,
        copyright=copyright,
        licenses=licenses,
        disclaimer=disclaimer,
        pmids=pmids,
        file=output,
    )


@document.command()
@click.argument('namespaces', nargs=-1)
@connection_option
@click.option('-p', '--path', type=click.File('r'), default=sys.stdin, help='Input BEL file path. Defaults to stdin.')
@click.option('-d', '--directory', help='Output folder. Defaults to current working directory {})'.format(os.getcwd()))
def serialize_namespaces(namespaces, connection: str, path, directory):
    """Parse a BEL document then serializes the given namespaces (errors and all) to the given directory."""
    from .definition_utils import export_namespaces

    graph = from_lines(path, manager=connection)
    export_namespaces(namespaces, graph, directory)


@io.command()
@graph_pickle_argument
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
def get_pmids(graph: BELGraph, output: TextIO):
    """Output PubMed identifiers from a graph to a stream."""
    for pmid in get_pubmed_identifiers(graph):
        click.echo(pmid, file=output)


if __name__ == '__main__':
    main()
