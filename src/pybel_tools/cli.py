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

import json
import logging
import sys
from getpass import getuser

import click
import hashlib
import os

from pybel import from_pickle, to_database, from_lines, from_url
from pybel.constants import (
    PYBEL_CONNECTION,
    SMALL_CORPUS_URL,
    LARGE_CORPUS_URL,
    get_cache_connection
)
from pybel.manager import Manager
from pybel.utils import get_version as pybel_version
from pybel.utils import parse_bel_resource
from pybel_tools.resources import deploy_directory
from .constants import GENE_FAMILIES, NAMED_COMPLEXES, DEFAULT_SERVICE_URL
from .definition_utils import (
    write_namespace,
    export_namespaces,
    hash_names,
    get_bel_resource_hash,
    get_bel_knowledge_hash,
)
from .document_utils import write_boilerplate
from .ioutils import convert_directory, upload_recursive, to_pybel_web, convert_paths
from .mutation.metadata import enrich_pubmed_citations
from .resources import get_namespace_history, get_annotation_history, get_knowledge_history
from .summary import get_pubmed_identifiers
from .utils import enable_cool_mode

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


@click.group(help="PyBEL-Tools Command Line Interface on {}\n with PyBEL v{}".format(sys.executable, pybel_version()))
@click.version_option()
def main():
    """PyBEL Tools Command Line Interface"""


@main.group()
@click.option('-c', '--connection', help='Cache connection. Defaults to {}'.format(get_cache_connection()))
@click.pass_context
def ensure(ctx, connection):
    """Utilities for ensuring data"""
    ctx.obj = Manager(connection=connection)


@ensure.command()
@click.option('--enrich-authors', is_flag=True)
@click.option('--use-edge-store', is_flag=True)
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.pass_obj
def small_corpus(manager, enrich_authors, use_edge_store, debug):
    """Caches the Selventa Small Corpus"""
    set_debug_param(debug)
    graph = from_url(SMALL_CORPUS_URL, manager=manager, citation_clearing=False, allow_nested=True)
    if enrich_authors:
        enrich_pubmed_citations(graph, manager=manager)
    manager.insert_graph(graph, store_parts=use_edge_store)


@ensure.command()
@click.option('--enrich-authors', is_flag=True)
@click.option('--use-edge-store', is_flag=True)
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.pass_obj
def large_corpus(manager, enrich_authors, use_edge_store, debug):
    """Caches the Selventa Large Corpus"""
    set_debug_param(debug)
    graph = from_url(LARGE_CORPUS_URL, manager=manager, citation_clearing=False, allow_nested=True)
    if enrich_authors:
        enrich_pubmed_citations(graph, manager=manager)
    manager.insert_graph(graph, store_parts=use_edge_store)


@ensure.command()
@click.option('--enrich-authors', is_flag=True)
@click.option('--use-edge-store', is_flag=True)
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.pass_obj
def gene_families(manager, enrich_authors, use_edge_store, debug):
    """Caches the HGNC Gene Family memberships"""
    set_debug_param(debug)
    graph = from_url(GENE_FAMILIES, manager=manager)
    if enrich_authors:
        enrich_pubmed_citations(graph, manager=manager)
    manager.insert_graph(graph, store_parts=use_edge_store)


@ensure.command()
@click.option('--enrich-authors', is_flag=True)
@click.option('--use-edge-store', is_flag=True)
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.pass_obj
def named_complexes(manager, enrich_authors, use_edge_store, debug):
    """Caches GO Named Protein Complexes memberships"""
    set_debug_param(debug)
    graph = from_url(NAMED_COMPLEXES, manager=manager)
    if enrich_authors:
        enrich_pubmed_citations(graph, manager=manager)
    manager.insert_graph(graph, store_parts=use_edge_store)


@main.group()
@click.option('-c', '--connection', help='Cache connection. Defaults to {}'.format(get_cache_connection()))
@click.option('--config', type=click.File('r'), help='Specify configuration JSON file')
@click.pass_context
def io(ctx, connection, config):
    """Upload and conversion utilities"""
    if config:
        file = json.load(config)
        ctx.obj = Manager(connection=file.get(PYBEL_CONNECTION, get_cache_connection()))
    else:
        ctx.obj = Manager(connection=connection)


@io.command()
@click.option('-p', '--path', default=os.getcwd())
@click.option('-r', '--recursive', is_flag=True,
              help='Recursively upload all gpickles in the directory given as the path')
@click.option('-s', '--skip-check-version', is_flag=True, help='Skip checking the PyBEL version of the gpickle')
@click.option('--to-service', is_flag=True, help='Sends to PyBEL web service')
@click.option('--service-url', help='Service location. Defaults to {}'.format(DEFAULT_SERVICE_URL))
@click.option('--exclude-directory-pattern', help="Pattern to match for bad directories")
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.pass_obj
def upload(manager, path, recursive, skip_check_version, to_service, service_url, exclude_directory_pattern, debug):
    """Quick uploader"""
    set_debug_param(debug)
    if recursive:
        log.info('uploading recursively from: %s', path)
        upload_recursive(path, connection=manager, exclude_directory_pattern=exclude_directory_pattern)
    else:
        graph = from_pickle(path, check_version=(not skip_check_version))
        if to_service:
            to_pybel_web(graph, service_url)
        else:
            to_database(graph, connection=manager)


@io.command()
@click.argument('path')
@click.option('-u', '--url', help='Service location. Defaults to {}'.format(DEFAULT_SERVICE_URL))
@click.option('-s', '--skip-check-version', is_flag=True, help='Skip checking the PyBEL version of the gpickle')
def post(path, url, skip_check_version):
    """Posts the given graph to the PyBEL Web Service via JSON"""
    graph = from_pickle(path, check_version=(not skip_check_version))
    to_pybel_web(graph, url)


@io.command()
@click.option('-u', '--enable-upload', is_flag=True, help='Enable automatic database uploading')
@click.option('--store-parts', is_flag=True, help='Automatically upload to database and edge store')
@click.option('--no-enrich-authors', is_flag=True, help="Don't enrich authors. Makes faster.")
@click.option('--no-enrich-genes', is_flag=True, help="Don't enrich HGNC genes")
@click.option('--no-enrich-go', is_flag=True, help="Don't enrich GO entries")
@click.option('--no-citation-clearing', is_flag=True, help='Turn off citation clearing')
@click.option('--allow-nested', is_flag=True, help="Enable lenient parsing for nested statements")
@click.option('-d', '--directory', default=os.getcwd(),
              help='The directory to search. Defaults to current working directory')
@click.option('-i', '--use-stdin', is_flag=True, help='Use stdin for paths')
@click.option('-w', '--send-pybel-web', is_flag=True, help='Send to PyBEL Web')
@click.option('--exclude-directory-pattern', help="Pattern to match for bad directories")
@click.option('--version-in-path', is_flag=True, help="Adds version to end of path string")
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
@click.option('-x', '--cool', is_flag=True, help='enable cool mode')
@click.pass_obj
def convert(manager, enable_upload, store_parts, no_enrich_authors, no_enrich_genes, no_enrich_go,
            no_citation_clearing, allow_nested, directory, use_stdin, send_pybel_web, exclude_directory_pattern,
            version_in_path, debug, cool):
    """Recursively walks the file tree and converts BEL scripts to gpickles. Optional uploader"""
    set_debug_param(debug)

    if cool:
        enable_cool_mode()

    kwargs = dict(
        connection=manager,
        upload=(enable_upload or store_parts),
        pickle=True,
        store_parts=store_parts,
        enrich_citations=(not no_enrich_authors),
        enrich_genes=(not no_enrich_genes),
        enrich_go=(not no_enrich_go),
        citation_clearing=(not no_citation_clearing),
        allow_nested=allow_nested,
        send=send_pybel_web,
        version_in_path=version_in_path,
    )

    if use_stdin:
        convert_paths(
            paths=(path for path in sys.stdin if path.endswith('.bel')),
            **kwargs
        )
    else:
        results = convert_directory(directory, exclude_directory_pattern=exclude_directory_pattern, **kwargs)

        for path, e in results:
            click.echo('FAILED {} {}'.format(path, e))


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


def _hash_helper(file):
    resource = parse_bel_resource(file)

    result = hash_names(
        resource['Values'],
        hash_function=hashlib.md5
    )

    click.echo(result)


@namespace.command()
@click.option('-f', '--file', type=click.File('r'), default=sys.stdin)
def semhash(file):
    """Semantic hash a namespace file"""
    _hash_helper(file)


@namespace.command()
@click.argument('namespace')
def history(namespace):
    """Hash all versions on Arty"""
    for path in get_namespace_history(namespace):
        h = get_bel_resource_hash(path.as_posix())
        click.echo('{}\t{}'.format(path, h))


@main.group()
def annotation():
    """Annotation file utilities"""


@annotation.command()
@click.argument('annotation')
def history(annotation):
    """Outputs the hashes for the annotation's versions"""
    for path in get_annotation_history(annotation):
        h = get_bel_resource_hash(path.as_posix())
        click.echo('{}\t{}'.format(path, h))


@annotation.command()
@click.option('-f', '--file', type=click.File('r'), default=sys.stdin)
def semhash(file):
    """Semantic hash a BEL annotation"""
    _hash_helper(file)


@main.group()
def document():
    """BEL document utilities"""


@document.command()
@click.argument('name')
def history(name):
    """Outputs the hashes for the BEL scripts' versions"""
    for path in get_knowledge_history(name):
        h = get_bel_knowledge_hash(path.as_posix())
        click.echo('{}\t{}'.format(path, h))


@document.command()
@click.argument('document-name')
@click.argument('contact')
@click.argument('description')
@click.argument('pmids', nargs=-1)
@click.option('--version')
@click.option('--copyright')
@click.option('--authors')
@click.option('--licenses')
@click.option('--output', type=click.File('wb'), default=sys.stdout)
def boilerplate(document_name, contact, description, pmids, version, copyright, authors, licenses, output):
    """Builds a template BEL document with the given PubMed identifiers"""
    write_boilerplate(
        document_name,
        contact,
        description,
        version,
        copyright,
        authors,
        licenses,
        pmids=pmids,
        file=output
    )


@document.command()
@click.argument('namespaces', nargs=-1)
@click.option('-c', '--connection', help='Cache connection. Defaults to {}'.format(get_cache_connection()))
@click.option('-p', '--path', type=click.File('r'), default=sys.stdin, help='Input BEL file path. Defaults to stdin.')
@click.option('-d', '--directory', help='Output directory. Defaults to current working directory')
def serialize_namespaces(namespaces, connection, path, directory):
    """Parses a BEL document then serializes the given namespaces (errors and all) to the given directory"""
    graph = from_lines(path, manager=connection)
    export_namespaces(namespaces, graph, directory)


@main.command()
@click.argument('directory')
@click.option('-v', '--debug', count=True, help="Turn on debugging. More v's, more debugging")
def upload_resources(directory, debug):
    """Uploads the resources in a directory to arty"""
    set_debug_param(debug)
    deploy_directory(directory)


@io.command()
@click.argument('path')
@click.option('-o', '--output', type=click.File('w'), default=sys.stdout)
def get_pmids(path, output):
    """Outputs PubMed identifiers from a graph to a stream"""
    graph = from_pickle(path)
    for pmid in get_pubmed_identifiers(graph):
        click.echo(pmid, file=output)


if __name__ == '__main__':
    main()
