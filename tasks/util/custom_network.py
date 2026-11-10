def make_node_id_map(g):
    mapping = {}
    for node in range(g.num_vertices()):
        mapping[g.vertex_properties['internal_id'][node]] = node
    return mapping

def make_node_id_map_pp(g):
    mapping = {}
    for node in range(g.num_vertices()):
        
        mapping[g.vertex_properties['internal_id'][node]] = node
    return mapping

def expand_compact_edges(edges):
    """
    Returns edges as [{"from": ..., "to": ...}, ...].

    Besides that list, accepts the compact format {"nodes": [id, ...], "edges": [[i, j], ...]},
    where i and j are positions in "nodes". Large networks are about a third of the size in it.
    Raises ValueError for any other structure.
    """
    if isinstance(edges, list):
        return edges
    if not isinstance(edges, dict) or set(edges) != {"nodes", "edges"}:
        raise ValueError('Edges must be a list of {"from": ..., "to": ...} or '
                         '{"nodes": [...], "edges": [[i, j], ...]}.')
    nodes, pairs = edges["nodes"], edges["edges"]
    if not isinstance(nodes, list) or not isinstance(pairs, list):
        raise ValueError('"nodes" and "edges" must be lists.')
    expanded = []
    for pair in pairs:
        if (not isinstance(pair, list) or len(pair) != 2
                or not all(type(i) is int and 0 <= i < len(nodes) for i in pair)):
            raise ValueError(f'Each edge must be a pair of positions in "nodes", got {pair!r}.')
        expanded.append({"from": nodes[pair[0]], "to": nodes[pair[1]]})
    return expanded


def add_edges(g, edge_list):
    """
    edge list is [{"fom":..., "to":...}, ...]
    """
    mapping = make_node_id_map(g)
    edge_id_list = []
    for edge in edge_list:
        a = mapping.get(edge['from'])
        b = mapping.get(edge['to'])
        if a is not None and b is not None:
            edge_id_list.append((a, b, 'protein-protein'))
    e_type = g.edge_properties["type"]
    g.add_edge_list(edge_id_list, eprops=[e_type])
    return g

def remove_ppi_edges(g):
    """
    removes all edges from g with type 'protein-protein'
    """
    edges_to_remove = []
    for edge in g.edges():
        if g.edge_properties["type"][edge] == 'protein-protein':
            edges_to_remove.append(edge)
    g.set_fast_edge_removal(fast=True)
    for edge in edges_to_remove:
        g.remove_edge(edge)
    g.set_fast_edge_removal(fast=False)
    return g

def filter_proteins(g, nodes_to_keep, drug_ids, seeds):
    """
    removes all nodes from g with internal_id but the ones in nodes_to_keep
    """
    mapping = make_node_id_map(g)
    nodes_to_keep = set(nodes_to_keep)
    seeds = set(seeds)
    nodes_to_remove = [node_id for node, node_id in mapping.items() if node not in nodes_to_keep and node_id not in drug_ids]
    g.remove_vertex(reversed(sorted(nodes_to_remove)), fast=True)
    
    # it is necessary to re-index
    new_seed_ids = []
    new_drug_ids = []
    for node in range(g.num_vertices()):
        if g.vertex_properties["type"][node] == 'drug':
            new_drug_ids.append(node)
        elif g.vertex_properties["type"][node] == 'protein':
            if g.vertex_properties["internal_id"][node] in seeds:
                new_seed_ids.append(node)
    return g, new_seed_ids, new_drug_ids