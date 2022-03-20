import argparse
import collections
import csv
import logging
import os
from os.path import join
import gzip
import time

from .algo import *

log = logging.getLogger()


class TaxIdError(ValueError):
    '''Taxonomy ID couldn't be determined.'''


def compressed_open(fname, *opts):
    return fname.endswith('.gz') and gzip.open(fname, *opts) or open(fname, *opts)


def maybe_compressed(fn):
    fn_gz = fn + '.gz'
    if os.path.exists(fn):
        return fn
    elif os.path.exists(fn_gz):
        return fn_gz
    else:
        raise FileNotFoundError(fn)


def add_taxonomy_arguments(parser=None):
    if not parser:
        parser = argparse.ArgumentParser()
    parser.add_argument('--tax-dir', help='Directory of taxonomy db.')
    parser.add_argument('--tax-nodes', help='Taxonomy db nodes file.')
    parser.add_argument('--tax-names', help='Taxonomy db names file.')
    parser.add_argument('--tax-merged', help='Taxonomy db merged file.')
    parser.add_argument('--tax-gis', nargs='+', help='Taxonomy db gi to taxid file.')
    return parser


class TaxonomyDb(object):

    def __init__(
        self,
        tax_dir=None,
        gis=None,
        nodes=None,
        names=None,
        merged=None,
        gis_paths=None,
        nodes_path=None,
        names_path=None,
        merged_path=None,
        load_gis=False,
        load_nodes=False,
        load_names=False,
        load_merged=True,
        scientific_names_only=None,
    ):
        self.tax_dir = tax_dir
        self.gis_paths = gis_paths
        self.nodes_path = nodes_path
        self.names_path = names_path
        self.merged_path = merged_path
        self.missing_parents = set()
        if load_gis:
            if gis:
                self.gis = gis
            else:
                if gis_paths:
                    self.gis_paths = gis_paths
                elif tax_dir:
                    self.gis_paths = [maybe_compressed(join(tax_dir, 'gi_taxid_nucl.dmp')),
                                      maybe_compressed(join(tax_dir, 'gi_taxid_prot.dmp'))]
                else:
                    self.gis_paths = None
                if self.gis_paths:
                    self.gis = {}
                    for gi_path in self.gis_paths:
                        start_time = time.time()
                        log.info('Loading taxonomy gis: %s', gi_path)
                        self.gis.update(self.load_gi_single_dmp(gi_path))
                        log.info('Loaded gi mapping: %.2fs', time.time() - start_time)
        if load_nodes:
            if nodes:
                self.ranks, self.parents = nodes
            else:
                if nodes_path:
                    self.nodes_path = nodes_path
                elif tax_dir:
                    self.nodes_path = maybe_compressed(join(tax_dir, 'nodes.dmp'))
                else:
                    self.nodes_path = None
                if self.nodes_path:
                    start_time = time.time()
                    log.info('Loading taxonomy nodes: %s', self.nodes_path)
                    self.ranks, self.parents = self.load_nodes(self.nodes_path)
                    log.info('Loaded taxonomy nodes: %.2fs', time.time() - start_time)
        if load_names:
            if names:
                self.names = names
            else:
                if names_path:
                    self.names_path = names_path
                elif tax_dir:
                    self.names_path = maybe_compressed(join(tax_dir, 'names.dmp'))
                else:
                    self.names_path = None
                if self.names_path:
                    start_time = time.time()
                    log.info('Loading taxonomy names: %s', self.names_path)
                    self.names = self.load_names(self.names_path, scientific_only=scientific_names_only)
                    log.info('Loaded taxonomy names: %.2fs', time.time() - start_time)
        if load_merged:
            if merged:
                self.merged = merged
                self.add_merged_links()
            else:
                if merged_path:
                    self.merged_path = merged_path
                elif tax_dir:
                    self.merged_path = maybe_compressed(join(tax_dir, 'merged.dmp'))
                else:
                    self.merged_path = None

                if self.merged_path:
                    start_time = time.time()
                    log.info('Loading merged nodes: %s', self.merged_path)
                    self.merged = self.load_merged(self.merged_path)
                    self.add_merged_links()
                    log.info('Loaded merged nodes: %.2fs', time.time() - start_time)

        self._children = None

    def load_gi_single_dmp(self, dmp_path):
        '''Load a gi->taxid dmp file from NCBI taxonomy.'''
        gi_array = {}
        with compressed_open(dmp_path) as f:
            for i, line in enumerate(f):
                gi, taxid = line.strip().split('\t')
                gi = int(gi)
                taxid = int(taxid)
                gi_array[gi] = taxid
                if (i + 1) % 1000000 == 0:
                    log.info('Loaded %s gis', i)
        return gi_array

    def load_names(self, names_db, scientific_only=None):
        '''Load the names.dmp file from NCBI taxonomy.'''
        scientific_only = scientific_only if scientific_only is not None else True
        if scientific_only:
            names = {}
        else:
            names = collections.defaultdict(list)
        for line in compressed_open(names_db):
            parts = line.strip().split('|')
            taxid = int(parts[0])
            name = parts[1].strip()
            #unique_name = parts[2].strip()
            class_ = parts[3].strip()
            if scientific_only:
                if class_ == 'scientific name':
                    names[taxid] = name
            else:
                if class_ == 'scientific name':
                    names[taxid].insert(0, name)
                else:
                    names[taxid].append(name)
        return names

    def load_nodes(self, nodes_db):
        '''Load ranks and parents arrays from NCBI taxonomy.'''
        ranks = {}
        parents = {}
        with compressed_open(nodes_db) as f:
            for line in f:
                parts = line.strip().split('|')
                taxid = int(parts[0])
                parent_taxid = int(parts[1])
                rank = parts[2].strip()
                #embl_code = parts[3].strip()
                #division_id = parts[4].strip()
                parents[taxid] = parent_taxid
                ranks[taxid] = rank
        return ranks, parents

    def load_merged(self, merged_db):
        '''Load taxonomy nodes that have been merged with other nodes.'''
        merged = {}
        with compressed_open(merged_db) as f:
            for line in f:
                parts = line.strip().split('|')
                from_taxid = int(parts[0])
                to_taxid = int(parts[1])
                merged[from_taxid] = to_taxid
        return merged

    def add_merged_links(self):
        '''Add entries in names, parents, ranks for merged nodes.'''
        for from_node, to_node in self.merged.items():
            if self.names:
                self.names[from_node] = self.names[to_node]
            if self.ranks:
                self.ranks[from_node] = self.ranks[to_node]
            if self.parents:
                self.parents[from_node] = self.parents[to_node]

    @classmethod
    def from_args(cls, parser_args, *args, **kwargs):
        return cls(tax_dir=parser_args.tax_dir,
                   nodes_path=parser_args.tax_nodes,
                   names_path=parser_args.tax_names,
                   merged_path=parser_args.tax_merged,
                   *args, **kwargs)

    @property
    def children(self):
        if self._children:
            return self._children
        self._children = parents_to_children(self.parents)
        return self._children

    def coverage_lca(self, query_ids, lca_percent=None):
        return coverage_lca(query_ids, self.parents, lca_percent=lca_percent,
                            missing_parents=self.missing_parents)

    def parent_path(self, query_id, root_id=1, cache=None, warn_missing=True):
        return parent_path(query_id, self.parents, root_id=root_id, cache=cache, warn_missing=warn_missing,
                           missing_parents=self.missing_parents)

    def kraken_dfs_report(self, taxa_hits, total_reads=None):
        '''Return a kraken compatible DFS report of taxa hits.

        Args:
        db: (TaxonomyDb) Taxonomy db.
        taxa_hits: (collections.Counter) # of hits per tax id.
        total_reads: (int) Total reads (add to unclassified if taxa_hits is less)

        Return:
        []str lines of the report
        '''
        total_hits = sum(taxa_hits.values())
        lines = []
        if total_reads is not None:
            if total_reads >= total_hits:
                unclassified_hits = total_reads - total_hits
                total_hits = total_reads
            elif total_reads < total_hits:
                raise Exception('Total reads given is < number of hits in the report.')
        else:
            unclassified_hits = 0

        unclassified_hits += taxa_hits.get(0, 0)
        unclassified_hits += taxa_hits.get(-1, 0)
        if total_hits:
            self.kraken_dfs(lines, taxa_hits, total_hits, 1, 0)

        if unclassified_hits > 0:
            percent_covered = '%.2f' % (unclassified_hits / total_hits * 100)
            lines.append(
                '\t'.join([
                    str(percent_covered), str(unclassified_hits), str(unclassified_hits), 'U', '0', 'unclassified'
                ])
            )
        return reversed(lines)


    def kraken_dfs(self, lines, taxa_hits, total_hits, taxid, level):
        '''Recursively do DFS for number of hits per taxa.'''
        cum_hits = num_hits = taxa_hits.get(taxid, 0)
        for child_taxid in self.children[taxid]:
            cum_hits += self.kraken_dfs(lines, taxa_hits, total_hits, child_taxid, level + 1)
        percent_covered = '%.2f' % (cum_hits / total_hits * 100)
        rank = kraken_rank_code(self.ranks[taxid])
        name = self.names[taxid]
        if cum_hits > 0:
            lines.append('\t'.join([percent_covered, str(cum_hits), str(num_hits), rank, str(taxid), '  ' * level + name]))
        return cum_hits


def kraken_rank_code(rank):
    '''Get the kraken-based short 1 letter rank code for named ranks.'''
    if rank == 'species':
        return 'S'
    elif rank == 'genus':
        return 'G'
    elif rank == 'family':
        return 'F'
    elif rank == 'order':
        return 'O'
    elif rank == 'class':
        return 'C'
    elif rank == 'phylum':
        return 'P'
    elif rank == 'kingdom':
        return 'K'
    elif rank == 'superkingdom':
        return 'D'
    else:
        return '-'


def taxa_hits_from_tsv(f, taxid_column=2):
    '''Return a counter of hits from tsv.'''
    c = collections.Counter()
    for row in csv.reader(f, delimiter='\t'):
        taxid = int(row[taxid_column - 1])
        c[taxid] += 1
    return c


def parents_to_children(parents):
    '''Convert an array of parents to lists of children for each parent.

    Returns:
      (dict[list]) Lists of children
    '''
    children = collections.defaultdict(list)
    for node, parent in parents.items():
        if node == 1:
            continue
        if parent != 0:
            children[parent].append(node)
    return children


def collect_children(children, original_taxids):
    '''Collect nodes with all children recursively.'''
    taxids = original_taxids
    while taxids:
        taxid = taxids.pop()
        yield taxid
        for child_taxid in children[taxid]:
            taxids.add(child_taxid)


def collect_parents(parents, taxids):
    '''Collect nodes with all parents recursively.'''
    # The root taxid node is 1
    yield 1
    taxids_with_parents = set([1])
    for taxid in taxids:
        while taxid not in taxids_with_parents:
            yield taxid
            taxids_with_parents.add(taxid)
            taxid = parents[taxid]
