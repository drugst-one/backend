from tasks.util.custom_edges import add_edges
from tasks.task_hook import TaskHook
import graph_tool as gt
from drugstone.models import *
from drugstone.serializers import *
import os
import networkx as nx
import community as community_louvain
from collections import Counter
import matplotlib.colors as mcolors
import graph_tool.util as gtu

import distinctipy

def generate_color_palette(num_colors):
    palette = distinctipy.get_colors(num_colors)

    hex_colors = [mcolors.rgb2hex(color) for color in palette]

    return hex_colors


def add_cluster_groups_to_config(config, partition):
    unique_clusters = set(partition.values())
    color_palette = generate_color_palette(len(unique_clusters))
    
    for cluster_id, color in zip(unique_clusters, color_palette):
        group_name = f"Cluster {cluster_id}"
        group_id = f"cluster{cluster_id}"
        
        if not config["node_groups"].get(group_id):
            config["node_groups"][group_id] = {
                "group_name": group_name,
                "color": {
                    "border": color,
                    "background": color,
                    "highlight": {
                        "border": color,
                        "background": color
                    }
                },
                "shape": "circle",
                "type": "gene",
                "border_width": 0,
                "border_width_selected": 0,
                "font": {
                    "color": "#FFFFFF" if is_dark_color(color) else "#000000",  # Setze Schriftfarbe auf Weiß, wenn Hintergrund dunkel ist
                    "size": 14,
                    "face": "arial",
                    "stroke_width": 0,
                    "stroke_color": "#ffffff",
                    "align": "center",
                    "bold": False,
                    "ital": False,
                    "boldital": False,
                    "mono": False
                },
                "shadow": True,
                "group_id": group_id
            }
    return config

def is_dark_color(color):
    r, g, b, _ = mcolors.to_rgba(color)
    brightness = (r * 299 + g * 587 + b * 114) / 1000
    return brightness < 0.5


def louvain_clustering(task_hook: TaskHook):
    r"""
    Perform Louvain clustering.

    Parameters
    ----------
    

    num_threads : int, optional (default: 1)
      Number of threads. Requires that graph_tool is compiled with OpenMP support.
      Should not be exposed in the frontend.
      
    ignore_isolated : bool, optional (default: True)

    Returns
    -------
    results : {
        "algorithm": "pathway_enrichment",
        "geneset": map_genesets[geneset_lowest_pvalue],
        "pathway": pathway_lowest_pvalue,
        "filteredDf": filtered_df.to_json(orient='records'),
        "backgroundMapping": background_mapping,
        "backgroundMappingReverse": background_mapping_reverse,
        "mapGenesets": map_genesets,
        "mapGenesetsReverse": map_genesets_reverse,
        "geneSetsDict": gene_sets_dict,
        "network": result,
        "table_view": table_view_results,
        "gene_interaction_dataset": ppi_dataset,
        "drug_interaction_dataset": pdi_dataset,
        "parameters": task_hook.parameters,
        "geneSets": genesets,
        "geneSetPathways": gene_set_terms_dict,
        "config": add_group_to_config(task_hook.parameters["config"]),
        "node_attributes":
            {
                "is_seed": isSeed,
            },
    }
    
    "algorithm": "pathway_enrichment"
    "geneset": The name of the geneset of the calculated pathway (in the beginning: pathway with lowest p-value).
    "pathway": The name of the calculated pathway (in the beginning: pathway with lowest p-value).
    "filteredDf": A JSON string containing the filtered DataFrame with the pathway enrichment results.
    "backgroundMapping": A dictionary that maps the internal node IDs to the original node IDs.
    "backgroundMappingReverse": A dictionary that maps the original node IDs to the internal node IDs.
    "mapGenesets": A dictionary that maps the internal geneset IDs to the original geneset IDs.
    "mapGenesetsReverse": A dictionary that maps the original geneset IDs to the internal geneset IDs.
    "geneSetsDict": A dictionary that contains the genesets and their pathways.
    "network": The calculated result for the chosen pathway (in the beginning: pathway with lowest p-value).
    "table_view": A list of dictionaries containing the enriched pathways.
    "gene_interaction_dataset": The gene interaction dataset.
    "drug_interaction_dataset": The drug interaction dataset.
    "parameters": The parameters of the task.
    "geneSets": A list of the genesets.
    "geneSetPathways": A dictionary that contains the genesets and their pathways.
    "config": The configuration of the task.
    "node_attributes": A one-element list containing a dictionary with the following 
        attributes for all nodes in the returned network:
        "is_seed": A flag that specifies whether the node is a seed.
    

    Notes
    -----
    This implementation is based on graph-tool, a very efficient Python package for network
    analysis with C++ backend and multi-threading support. Installation instructions for graph-tool
    can be found at https://git.skewed.de/count0/graph-tool/-/wikis/installation-instructions.

    References
    ----------
    ..  [1] “Community Detection for NetworkX’s Documentation — Community Detection for NetworkX 2 Documentation.” Accessed April 8, 2024. https://python-louvain.readthedocs.io/en/latest/.
        [2] “NetworkX — NetworkX Documentation.” Accessed April 8, 2024. https://networkx.org/.
        [3] Roberts, Jack, Jon Crall, Kian-Meng Ang, and Yannick Brandt. “Alan-Turing-Institute/Distinctipy: V1.3.4.” Zenodo, January 10, 2024. https://doi.org/10.5281/zenodo.10480933.

    """
    

    # Type: list of str
    # Semantics: Names of the seed proteins. Use UNIPROT IDs for host proteins, and
    #            names of the for SARS_CoV2_<IDENTIFIER> (e.g., SARS_CoV2_ORF6) for
    #            virus proteins.
    # Reasonable default: None, has to be selected by user via "select for analysis"
    #            utility in frontend.
    # Acceptable values: UNIPROT IDs, identifiers of viral proteins.
    seeds = task_hook.parameters["seeds"]


    # Type: int.
    # Semantics: Number of threads used for running the analysis.
    # Example: 1.
    # Reasonable default: 1.
    # Note: We probably do not want to expose this parameter to the user.
    num_threads = task_hook.parameters.get("num_threads", 1)

    ppi_dataset = task_hook.parameters.get("ppi_dataset")

    pdi_dataset = task_hook.parameters.get("pdi_dataset")

    id_space = task_hook.parameters["config"].get("identifier", "symbol")

    custom_edges = task_hook.parameters.get("custom_edges", False)
    
    ignore_isolated = task_hook.parameters.get("ignore_isolated", True)
    
    # Parsing input file.
    task_hook.set_progress(1 / 4.0, "Parsing input.")
    
    filename = f"{id_space}_{ppi_dataset['name']}-{pdi_dataset['name']}"
    if ppi_dataset['licenced'] or pdi_dataset['licenced']:
        filename += "_licenced"
    filename = os.path.join(task_hook.data_directory, filename + ".gt")
    g = gt.load_graph(filename)
    if custom_edges:
        edges = task_hook.parameters.get("input_network")['edges']
        g = add_edges(g, edges)
        
    node_name_attribute = "internal_id"
    node_mapping = {}
    node_mapping_reverse = {}
    for seed in seeds:
        found_node = int(gtu.find_vertex(g, prop=g.vertex_properties[node_name_attribute], match=seed)[0])
        node_mapping[seed] = found_node
        node_mapping_reverse[found_node] = seed
        
    all_nodes_int = set([int(node_mapping[gene]) for gene in seeds if gene in node_mapping])
    edges_unique = set()
    for node in seeds:
        for neighbor in g.get_all_neighbors(node_mapping[node]):
            if int(neighbor) > int(node_mapping[node]) and int(neighbor) in all_nodes_int:
                first_key = next(iter(node_mapping_reverse))
                if isinstance(first_key, int):
                    neighbor_key = int(neighbor)
                else:
                    neighbor_key = str(int(neighbor))
                edges_unique.add((node, node_mapping_reverse[neighbor_key]))

    
    # Set number of threads if OpenMP support is enabled.
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(num_threads)
    
    edges = [{"from": source, "to":target} for source, target in edges_unique]
    nodes = task_hook.parameters.get("input_network")['nodes']
    
    G = nx.Graph()

    isSeed = {}
    seedSet = set(seeds)
    for node in nodes:
        if node["id"] in seedSet:
            isSeed[node["id"]] = True
            G.add_node(node_for_adding=node["id"])
    # find all edges in the network that are connected to a seed node

    for edge in edges:
        if edge["from"] in seedSet and edge["to"] in seedSet:
            G.add_edge(edge["from"], edge["to"])
            
    if ignore_isolated:
        # delete isolated nodes from networkx graph
        isolated_nodes = list(nx.isolates(G))
        G.remove_nodes_from(isolated_nodes)

    task_hook.set_progress(3 / 4.0, "Perform Louvain clustering.")

    partition = community_louvain.best_partition(G)
    
    task_hook.set_progress(3 / 4.0, "Parse clustering results.")

    config = add_cluster_groups_to_config(task_hook.parameters["config"], partition)
    task_hook.parameters["config"] = config
   
    table_view_results = []
    total_nodes = len(partition)
    cluster_counts = Counter(partition.values())
    for cluster_id, count in cluster_counts.items():
        percentage = round((count / total_nodes) * 100, 2)
        table_view_results.append({
            "cluster_id": cluster_id,
            "count": str(count) + "/"+ str(total_nodes),
            "percentage": percentage,
        })
    
    filtered_nodes = []
    for node in nodes:
        if node["id"] in seedSet:
            if node["id"] in partition:
                cluster = partition[node["id"]]
                group_id = f"cluster{cluster}"
                node["group"] = group_id
                node["cluster"] = str(cluster)
                filtered_nodes.append(node)
            else:
                # node in seeds but was isolated and those are ignored -> keep old group
                node["cluster"] = "none"
                filtered_nodes.append(node)
                


    # return the results.
    task_hook.set_progress(4 / 4.0, "Returning results.")
    
    result = {
        "nodes": filtered_nodes,
        "edges": [{"from": str(u), "to": str(v)} for u, v in G.edges()],
    }
    task_hook.parameters["algorithm"] = "louvain-clustering"
    task_hook.set_results({
        "algorithm": "louvain_clustering",
        "network":result,
        "table_view": table_view_results,
        "parameters": task_hook.parameters,
        "gene_interaction_dataset": ppi_dataset,
        "drug_interaction_dataset": pdi_dataset,
        "node_attributes":
            {
                "is_seed": isSeed,
            },
    })
