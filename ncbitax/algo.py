import collections
import logging

log = logging.getLogger()

def parent_path(query_id, parents, root_id=1, cache=None, warn_missing=True,
                missing_parents=None):
    '''Calculate the path to the root taxon or a specific taxon.

    Args:
      missing_parents (set): Cached taxids with missing parents to suppress re-warning.
    Return:
      (int) List of tax ids, else None if no valid path found
    '''
    if cache is not None and query_id in cache:
        return cache[query_id]

    path = []
    while query_id != root_id:
        path.append(query_id)
        if not parents.get(query_id):
            if warn_missing:
                if missing_parents is not None:
                    if query_id not in missing_parents:
                        missing_parents.add(query_id)
                        log.warn('Parent for query id: {} missing'.format(query_id))
                else:
                    log.warn('Parent for query id: {} missing'.format(query_id))
            break
        query_id = parents[query_id]
    result = None
    if query_id == root_id:
        path.append(root_id)
        result = list(path[1:])
    if cache is not None:
        cache[query_id] = result
    return result

def coverage_lca(query_ids, parents, lca_percent=None, missing_parents=None):
    '''Calculate the lca that will cover at least this percent of queries.

    Args:
      query_ids: []int list of nodes.
      parents: []int array of parents.
      lca_percent: (float) Cover at least this percent of queries.
      missing_parents (set): Cached taxids with missing parents to suppress re-warning.

    Return:
      (int) LCA
    '''
    if len(query_ids) == 1:
        return query_ids[0]
    lca_percent = lca_percent or 100
    lca_needed = lca_percent / 100 * len(query_ids)
    paths = []
    for query_id in query_ids:
        path = parent_path(query_id, parents, missing_parents=missing_parents)
        if path:
            paths.append(list(reversed([query_id] + path)))
    if not paths:
        return

    last_common = 1
    max_path_length = max(len(path) for path in paths)
    for level in range(max_path_length):
        valid_paths = (path for path in paths if len(path) > level)
        max_query_id, hits_covered = collections.Counter(path[level] for path in valid_paths).most_common(1)[0]
        if hits_covered >= lca_needed:
            last_common = max_query_id
        else:
            break
    return last_common
