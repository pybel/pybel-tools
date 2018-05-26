# -*- coding: utf-8 -*-

"""
Module that contains the command line app

Why does this file exist, and why not put this in __main__?
You might be tempted to import things from __main__ later, but that will cause
problems--the code will get executed twice:
 - When you run `python3 -m pybel_tools` python will execute
   ``__main__.py`` as a script. That means there won't be any
   ``pybel_tools.__main__`` in ``sys.modules``.
 - When you import __main__ it will get executed again (as a module) because
   there's no ``pybel_tools.__main__`` in ``sys.modules``.
Also see (1) from http://click.pocoo.org/5/setuptools/#setuptools-integration
"""

from __future__ import print_function

import logging
import os
import sys
from getpass import getuser

import click

from pybel.constants import (
    LARGE_CORPUS_URL, NAMESPACE_DOMAIN_OTHER, SMALL_CORPUS_URL, get_cache_connection,
)
from pybel.resources.definitions import (
    parse_bel_resource, write_namespace,
)
from pybel.struct import get_pubmed_identifiers, union
from pybel.utils import get_version as pybel_version
from .constants import DEFAULT_SERVICE_URL, NAMED_COMPLEXES_URL
from .ioutils import convert_paths, get_paths_recursive, upload_recursive
from .utils import enable_cool_mode, get_version

log = logging.getLogger(__name__)

datefmt = '%H:%M:%S'
fmt = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

user_dump_path = os.path.join(os.path.expanduser('~'), '.pybel', 'data', 'users.tsv')


def set_debug(level):
    logging.basicConfig(level=level, format=fmt, datefmt=datefmt)
    pybel_log = logging.getLogger('pybel')
    pybel_log.setLevel(level)
    pbt_log = logging.getLogger('pybel_tools')
    pbt_log.setLevel(level)
    log.setLevel(level)


def set_debug_param(debug):
    if debug == 1:
        set_debug(20)
    elif debug == 2:
        set_debug(10)


@click.group(help="PyBEL-Tools v{} Command Line Interface on {}\n with PyBEL v{}".format(get_version(), sys.executable,
                                                                                         pybel_version()))
@click.version_option()
def main():
    """PyBEL Tools Command Line Interface"""


@main.group()
@click.option('-c', '--connection', help='Cache connection. Defaults to {}'.format(get_cache_connection()))
@click.pass_context
def ensure(ctx, connection):
    """Utilities for ensuring data"""
    from pybel.manager.cache_manager import Manager
    ctx.obj = Manager(connection=connection)


def help_ensure(debug, url, enrich_authors, manager):
    set_debug_param(debug)
    from pybel import from_url
    graph = from_url(url, manager=manager, citation_clearing=False, allow_nested=True)
    if enrich_authors:
        from .mutation.metadata import enrich_pubmed_citations
        enrich_pubmed_citations(graph, manager=manager)
    manager.insert_graph(graph, store_parts=True)


@ensure.command()
@click.option('--enrich-authors', is_flag=True)
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.pass_obj
def small_corpus(manager, enrich_authors, debug):
    """Caches the Selventa Small Corpus"""
    help_ensure(debug, SMALL_CORPUS_URL, enrich_authors, manager)


@ensure.command()
@click.option('--enrich-authors', is_flag=True)
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.pass_obj
def large_corpus(manager, enrich_authors, debug):
    """Caches the Selventa Large Corpus"""
    help_ensure(debug, LARGE_CORPUS_URL, enrich_authors, manager)


@ensure.command()
@click.option('--enrich-authors', is_flag=True)
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.pass_obj
def named_complexes(manager, enrich_authors, debug):
    """Caches GO Named Protein Complexes memberships"""
    help_ensure(debug, NAMED_COMPLEXES_URL, enrich_authors, manager)


@main.group()
@click.option('-c', '--connection', help='Cache connection. Defaults to {}'.format(get_cache_connection()))
@click.pass_context
def io(ctx, connection):
    """Upload and conversion utilities"""
    from pybel.manager.cache_manager import Manager
    ctx.obj = Manager(connection=connection)


@io.command()
@click.option('-p', '--path', help='Path or directory. Defaults to {}'.format(os.getcwd()), default=os.getcwd())
@click.option('-s', '--skip-check-version', is_flag=True, help='Skip checking the PyBEL version of the gpickle')
@click.option('--to-service', is_flag=True, help='Sends to PyBEL web service')
@click.option('--service-url', help='Service location. Defaults to {}'.format(DEFAULT_SERVICE_URL))
@click.option('--exclude-directory-pattern', help="Pattern to match for bad directories")
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.pass_obj
def upload(manager, path, skip_check_version, to_service, service_url, exclude_directory_pattern, debug):
    """Upload gpickles"""
    set_debug_param(debug)

    if os.path.isdir(path):
        log.info('uploading recursively from: %s', path)
        upload_recursive(path, connection=manager, exclude_directory_pattern=exclude_directory_pattern)

    elif os.path.isfile(path):
        from pybel import from_pickle
        graph = from_pickle(path, check_version=(not skip_check_version))
        if to_service:
            from pybel import to_web
            to_web(graph, service_url)
        else:
            from pybel import to_database
            to_database(graph, connection=manager, store_parts=True)


@io.command()
@click.argument('path')
@click.option('-u', '--url', help='Service location. Defaults to {}'.format(DEFAULT_SERVICE_URL))
@click.option('-s', '--skip-check-version', is_flag=True, help='Skip checking the PyBEL version of the gpickle')
def post(path, url, skip_check_version):
    """Posts the given graph to the PyBEL Web Service via JSON"""
    from pybel import from_pickle, to_web
    graph = from_pickle(path, check_version=(not skip_check_version))
    to_web(graph, url)


@io.command()
@click.option('-u', '--enable-upload', is_flag=True, help='Enable automatic database uploading')
@click.option('--enrich-citations', is_flag=True, help="Enrich citations and authors. Makes slower.")
@click.option('--no-citation-clearing', is_flag=True, help='Turn off citation clearing')
@click.option('-n', '--allow-nested', is_flag=True, help="Enable lenient parsing for nested statements")
@click.option('-d', '--directory', default=os.getcwd(),
              help='The directory to search. Defaults to cwd: {}'.format(os.getcwd()))
@click.option('-i', '--use-stdin', is_flag=True, help='Use stdin for paths')
@click.option('-w', '--send-pybel-web', is_flag=True, help='Send to PyBEL Web')
@click.option('--exclude-directory-pattern', help="Pattern to match for bad directories")
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.option('--uncool', is_flag=True, help='disable cool mode')
@click.option('-c', '--collate', is_flag=True, help='Make a single graph')
@click.pass_obj
def convert(manager, enable_upload, enrich_citations, no_citation_clearing, allow_nested, directory, use_stdin,
            send_pybel_web, exclude_directory_pattern, debug, uncool, collate):
    """Recursively walks the file tree and converts BEL scripts to gpickles. Optional uploader"""
    set_debug_param(debug)

    if not uncool:
        enable_cool_mode()

    if use_stdin:
        paths = (path for path in sys.stdin if path.endswith('.bel'))
    else:
        paths = get_paths_recursive(directory, exclude_directory_pattern=exclude_directory_pattern)

    successes, failures = convert_paths(
        paths=paths,
        connection=manager,
        upload=enable_upload,
        enrich_citations=enrich_citations,
        citation_clearing=(not no_citation_clearing),
        allow_nested=allow_nested,
        send=send_pybel_web,
        use_tqdm=True
    )

    for path, e in failures:
        click.echo('FAILED {} {}'.format(path, e))

    if collate:
        from pybel import to_pickle

        collated_graph = union(successes.values())
        collated_path = os.path.join(directory, '{}-collated.gpickle'.format(directory))
        to_pickle(collated_graph, file=collated_path)


@io.command()
@click.option('-d', '--directory', default=os.getcwd(),
              help='The directory to search. Defaults to cwd: {}'.format(os.getcwd()))
@click.option('-n', '--name', help='The name of the file. Defaults to directory name')
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.pass_obj
def merge_directory(manager, directory, name, debug):
    """Parses all BEL files in a directory and outputs it"""
    set_debug_param(debug)

    name = name or '{}-merged.gpickle'.format(directory)
    path = os.path.join(directory, name)
    if os.path.exists(path):
        click.echo('Path already exists. Quitting. [{}]'.format(path))

    from . import from_directory
    from pybel import to_pickle

    enable_cool_mode()

    graph = from_directory(directory, connection=manager)
    to_pickle(graph, file=path)


@main.group()
def namespace():
    """Namespace file utilities"""


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
@click.option('--functions')
@click.option('--output', type=click.File('w'), default=sys.stdout)
@click.option('--value-prefix', default='')
def write(name, keyword, domain, citation, author, description, species, version, contact, license, values,
          functions, output, value_prefix):
    """Builds a namespace from items"""
    write_namespace(
        name, keyword, domain, author, citation, values,
        namespace_description=description,
        namespace_species=species,
        namespace_version=version,
        author_contact=contact,
        author_copyright=license,
        functions=functions,
        file=output,
        value_prefix=value_prefix
    )


@namespace.command()
@click.option('-f', '--file', type=click.File('r'), default=sys.stdin, help="Path to input BEL Namespace file")
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout,
              help="Path to output converted BEL Annotation file")
def convert_to_annotation(file, output):
    """Convert a namespace file to an annotation file"""
    from pybel.resources.definitions import write_annotation

    resource = parse_bel_resource(file)

    write_annotation(
        keyword=resource['Namespace']['Keyword'],
        values={k: '' for k in resource['Values']},
        citation_name=resource['Citation']['NameString'],
        description=resource['Namespace']['DescriptionString'],
        file=output
    )


@main.group()
def annotation():
    """Annotation file utilities"""


@annotation.command()
@click.option('-f', '--file', type=click.File('r'), default=sys.stdin, help="Path to input BEL Namespace file")
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout,
              help="Path to output converted BEL Namespace file")
@click.option('--keyword', help="Set custom keyword. useful if the annotion keyword is too long")
def convert_to_namespace(file, output, keyword):
    """Convert an annotation file to a namespace file"""
    resource = parse_bel_resource(file)
    write_namespace(
        namespace_keyword=(keyword or resource['AnnotationDefinition']['Keyword']),
        namespace_name=resource['AnnotationDefinition']['Keyword'],
        namespace_description=resource['AnnotationDefinition']['DescriptionString'],
        author_name='Charles Tapley Hoyt',
        namespace_domain=NAMESPACE_DOMAIN_OTHER,
        values=resource['Values'],
        citation_name=resource['Citation']['NameString'],
        file=output
    )


@main.group()
def document():
    """BEL document utilities"""


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
    """Builds a template BEL document with the given PubMed identifiers"""
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
        file=output
    )


@document.command()
@click.argument('namespaces', nargs=-1)
@click.option('-c', '--connection', help='Cache connection. Defaults to {}'.format(get_cache_connection()))
@click.option('-p', '--path', type=click.File('r'), default=sys.stdin, help='Input BEL file path. Defaults to stdin.')
@click.option('-d', '--directory', help='Output folder. Defaults to current working directory {})'.format(os.getcwd()))
def serialize_namespaces(namespaces, connection, path, directory):
    """Parses a BEL document then serializes the given namespaces (errors and all) to the given directory"""
    from pybel import from_lines
    from .definition_utils import export_namespaces

    graph = from_lines(path, manager=connection)
    export_namespaces(namespaces, graph, directory)


@io.command()
@click.argument('path')
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
def get_pmids(path, output):
    """Outputs PubMed identifiers from a graph to a stream"""
    from pybel import from_pickle
    graph = from_pickle(path)
    for pmid in get_pubmed_identifiers(graph):
        click.echo(pmid, file=output)


if __name__ == '__main__':
    main()
