import copy
from collections import defaultdict
import json
from typing import List, Tuple, Set, OrderedDict
from django.db.models import Q
from django.db.models.functions import Upper
from drugstone import models
from drugstone.models import Protein, EnsemblGene, Task
from drugstone.serializers import ProteinProteinInteractionSerializer, ProteinSerializer


MAP_ID_SPACE_COMPACT_TO_DRUGSTONE = {
    'symbol:': 'symbol',
    'uniprot:': 'uniprot',
    'ensg:': 'ensg',
    'ncbigene:': 'entrez',
    'ensembl:': 'ensg',
    'entrez:': 'entrez'
}

def get_ppi_ds(source, licenced):
    ds = models.PPIDataset.objects.filter(name__iexact=source, licenced=licenced).last()
    if ds is None and licenced:
        return get_ppi_ds(source, False)
    return ds


PROTEIN_ATTRIBUTES = {"symbol": "symbol", "uniprot": "uniprot", "ensg": "ensg", "ensembl": "ensg",
                      "entrez": "entrez", "ncbigene": "entrez"}


def proteins_by_identifier(node_ids, protein_attribute: str):
    """Returns the proteins matching the IDs in one query, case-insensitive except for Entrez IDs.

    UPPER() values are compared with IN, which uses the functional indexes; one OR condition per ID
    cannot use an index and gets slow for thousands of IDs.
    """
    upper_ids = {str(n).upper() for n in node_ids}
    if protein_attribute == "symbol":
        return Protein.objects.annotate(lookup=Upper("gene")).filter(lookup__in=upper_ids)
    if protein_attribute == "uniprot":
        return Protein.objects.annotate(lookup=Upper("uniprot_code")).filter(lookup__in=upper_ids)
    if protein_attribute == "ensg":
        protein_ids = EnsemblGene.objects.annotate(lookup=Upper("name")).filter(
            lookup__in=upper_ids).values_list("protein_id", flat=True)
        return Protein.objects.filter(id__in=list(protein_ids))
    return Protein.objects.filter(entrez__in=[str(n) for n in node_ids])


def query_proteins_by_identifier(node_ids: Set[str], identifier: str, reviewed: bool) -> Tuple[List[dict], str]:
    """Queries the django database Protein table given a list of identifiers (node_ids) and a identifier name
    (identifier).
    The identifier name represents any protein attribute, e.g. uniprot or symbol.
    The identifier names vary from the Protein table names since they are the strings which are set by the user
    in the frontend, for readability they were changes from the original backend attributes.

    Args:
        node_ids (list): List of protein or gene identifiers. Note: Do not mix identifiers.
        identifier (str): Can be one of "symbol", "ensg", "uniprot"

    Returns:
        Tuple[List[dict], str]:
            Returns list of serialized protein entries for all matched IDs
            Returns name of backend attribute of Protein table
    """
    # query protein table
    if (len(node_ids) == 0):
        return list(), identifier
    protein_attribute = PROTEIN_ATTRIBUTES.get(identifier)
    if protein_attribute is None:
        return [], identifier
    proteins = proteins_by_identifier(node_ids, protein_attribute)
    if reviewed:
        proteins = proteins.filter(isReviewed=True)
    # load related rows in bulk instead of one query per protein (cellular components here, ENSG in the serializer)
    node_objects = proteins.prefetch_related("cellular_components", "ensg")

    
    cc_to_node_ids = {}
    for node in node_objects:
        cellular_components = node.cellular_components.all()
        components = []
        for cc in cellular_components:
            cc_string = cc.go_code + ":" + cc.display_name + ":" + cc.layer
            components.append(cc_string)
        
        node_id = ''
        if protein_attribute == 'symbol':
            node_id = node.gene
        elif protein_attribute == 'uniprot' or protein_attribute == 'ensg':
            node_id = node.uniprot_code
        elif protein_attribute == 'entrez':
            node_id = node.entrez
        cc_to_node_ids[node_id] = components
        
    nodes = list()
    node_map = defaultdict(list)
    if protein_attribute == 'ensg':
        for node in ProteinSerializer(many=True).to_representation(node_objects):
            for ensembl_id in node.get(protein_attribute):
                if ensembl_id.upper() in node_ids:
                    node = copy.copy(node)
                    node[identifier] = ensembl_id
                    id_node = node.get("uniprot")
                    if id_node in cc_to_node_ids:
                        node["cellular_component"] = cc_to_node_ids[id_node]
                    else:
                        node["cellular_component"] = []
                    node_map[ensembl_id].append(node)
    else:
        for node in ProteinSerializer(many=True).to_representation(node_objects):
            id_node = node.get(protein_attribute)
            if id_node in cc_to_node_ids:
                node["cellular_component"] = cc_to_node_ids[id_node]
            else:
                node["cellular_component"] = []
            node_map[node.get(protein_attribute)].append(node)
    for node_id, entries in node_map.items():
        nodes.append(aggregate_nodes(entries))
        
    layer_ids = {'GO:0005737': "Cytoplasm", 'GO:0005634': "Nucleus", 'GO:0005576': "Extracellular", 'GO:0009986': "Cell surface", 'GO:0005886': "Plasma membrane"}
    for node in nodes:
        ccs = node.get("cellular_component", [])
        if len(ccs) > 0:
            layers = set()
            for cc in ccs:
                splitted = cc.split(":")
                if len(splitted) == 4:
                    # go could be mapped
                    layer = "GO:" + cc.split(":")[3]
                    layers.add(layer)
            if len(layers) == 1:
                node["layer"] = layer_ids[list(layers)[0]]
            elif len(layers) == 0:
                node["layer"] = "Other"
            else:
                layer_names = [layer_ids[layer] for layer in layers]
                node["layer"] = f"Multiple ({', '.join(layer_names)})"
        else:
            node["layer"] = "Unknown"
            
    return nodes, protein_attribute

def get_protein_ids(id_space, proteins):
    attribute = PROTEIN_ATTRIBUTES.get(id_space)
    if attribute is None:
        return set()
    ids = set()
    for p in proteins:
        # ensg is a list, a protein can have several Ensembl genes
        value = p[attribute]
        ids |= set(value) if isinstance(value, list) else {value}
    return ids


def clean_proteins_from_compact_notation(node_ids: Set[str], identifier: str) -> List[str]:
    """Queries the django database Protein table given a list of identifiers (node_ids) and a identifier name
    (identifier).
    The identifier name represents any protein attribute, e.g. uniprot or symbol.
    The identifier names vary from the Protein table names since they are the strings which are set by the user
    in the frontend, for readability they were changes from the original backend attributes.

    Args:
        node_ids (list): List of protein or gene identifiers. Note: Do not mix identifiers.
        identifier (str): Can be one of "symbol", "ensg", "uniprot"

    Returns:
        Tuple[List[dict], str]:
            Returns list of serialized protein entries for all matched IDs
            Returns name of backend attribute of Protein table
    """
    # query protein table
    if len(node_ids) == 0:
        return list()

    symbol_set, ensg_set, uniprot_set, entrez_set = set(), set(), set(), set()

    id_map = {
        'symbol:': symbol_set,
        'uniprot:': uniprot_set,
        'ensg:': ensg_set,
        'ncbigene:': entrez_set,
        'ensembl:': ensg_set,
        'entrez:': entrez_set
    }
    clean_ids = set()
    for node_id in node_ids:
        added = False
        for id_space in id_map.keys():
            if node_id.startswith(id_space):
                id_map[id_space].add(node_id[len(id_space):].upper())
                added = True
                break
        if not added:
            clean_ids.add(node_id)

    for id_space, ids in id_map.items():
        if len(ids) == 0:
            continue
        protein_attribute = PROTEIN_ATTRIBUTES.get(id_space[:-1])
        if protein_attribute is None:
            continue
        proteins = ProteinSerializer(many=True).to_representation(
            proteins_by_identifier(ids, protein_attribute).prefetch_related("ensg"))
        # if protein could not be mapped
        clean_ids_temp = get_protein_ids(identifier, proteins)
        clean_ids_temp.discard('')
        target_attribute = PROTEIN_ATTRIBUTES.get(identifier)
        missing = [p for p in proteins if target_attribute and p[target_attribute] in ('', [])]
        if missing:
            # at least one protein could not be found in id space, use original id as placeholder
            if protein_attribute == 'ensg':
                # use the requested Ensembl ID, not all Ensembl genes of the protein
                ids_placeholder = {e for p in missing for e in p['ensg'] if e.upper() in ids}
            else:
                ids_placeholder = get_protein_ids(protein_attribute, missing)
            clean_ids_temp |= ids_placeholder
        clean_ids |= clean_ids_temp

    return list(clean_ids)


def aggregate_nodes(nodes: List[OrderedDict]):
    node = defaultdict(set)
    for n in nodes:
        for key, value in n.items():
            if isinstance(value, list):
                for e in value:
                    if e is not None and len(e) > 0:
                        node[key].add(e)
            elif (isinstance(value, bool)) or (value is not None and len(value) > 0):
                node[key].add(value)
    return {k: list(v) for k, v in node.items()}

def update_result(result, token: str):
    task = Task.objects.get(token=token)
    task.result = json.dumps(result)
    task.save()
    

def fetch_node_information(nodes, identifier, reviewed):
    id_map = {}
    nodes_clean = []
    for node in nodes:
        if not node["id"]:
            # skip empty node id ''
            continue
        upper = node["id"].upper()
        id_map[upper] = node["id"]
        node["id"] = upper
        nodes_clean.append(node)
    nodes = nodes_clean

    # extract ids for filtering
    node_ids = set([node["id"] for node in nodes])

    # query protein table
    nodes_mapped, id_key = query_proteins_by_identifier(node_ids, identifier, reviewed)

    # change data structure to dict in order to be quicker when merging
    nodes_mapped_dict = {}
    for node in nodes_mapped:
        if id_key in node:
            for id in node[id_key]:
                nodes_mapped_dict[id.upper()] = node
        # TODO find solution if target id space is empty
        # else:
        #     nodes_mapped_dict[node['id'].upper()] = node

    # merge fetched data with given data to avoid data loss
    for node in nodes:
        node["drugstoneType"] = "other"
        if node["id"] in nodes_mapped_dict:
            node["cellular_component"] = []
            node.update(nodes_mapped_dict[node["id"]])
            node["drugstoneType"] = "protein"
        node["id"] = id_map[node["id"]]
    return nodes

def fetch_edges_from_input(dataset: str, licenced: bool, edges: list) -> list:
    dataset_object = get_ppi_ds(dataset, licenced)
    edge_keys = set()
    for edge in edges:
        from_node = edge.get("from")
        to_node = edge.get("to")
        if from_node and to_node:
            if to_node.startswith('d') or from_node.startswith('d'):
                continue
            from_node = from_node[1:] if from_node.startswith('p') else from_node
            to_node = to_node[1:] if to_node.startswith('p') else to_node
            edge_keys.add((from_node, to_node))

    protein_ids = {int(node) for edge in edge_keys for node in edge}
    interaction_objects = models.ProteinProteinInteraction.objects.filter(
        Q(ppi_dataset=dataset_object) &
        Q(from_protein_id__in=protein_ids) &
        Q(to_protein_id__in=protein_ids)
    )
    serialized_data = ProteinProteinInteractionSerializer(many=True).to_representation(interaction_objects)
    serialized_data = [
        {
            'from': entry['protein_a'], 
            'to': entry['protein_b'], 
            **{k: v for k, v in entry.items() if k not in ['protein_a', 'protein_b']}
        }
        for entry in serialized_data
    ]
    found_edge_map = {
        (edge['from'], edge['to']): edge
        for edge in serialized_data
    }

    edges_to_return = []
    is_omnipath = dataset == "OmniPath"
    for edge in edges:
        from_node = edge.get("from")
        to_node = edge.get("to")

        if from_node and to_node:
            found_exact = (from_node, to_node) in found_edge_map
            found_reverse = (to_node, from_node) in found_edge_map
            
            if found_exact:
                # Exact match found - use DB edge with all properties
                edges_to_return.append(found_edge_map[(from_node, to_node)])
            
            if found_reverse:
                if is_omnipath:
                    # For OmniPath, both directions can exist (e.g., A->B stimulation, B->A inhibition)
                    # If reverse direction exists, also add it (even if exact match was found)
                    edges_to_return.append(found_edge_map[(to_node, from_node)])
                else:
                    # For non-OmniPath datasets (undirected), only use reverse if exact not found
                    if not found_exact:
                        edges_to_return.append(found_edge_map[(to_node, from_node)])
            
            if not found_exact and not found_reverse:
                # Edge not found in DB - return original edge
                edges_to_return.append(edge)

    return edges_to_return

def fetch_edges_for_proteins(ppi_dataset_name: str, licenced: bool, uniprot_codes: set, require_both_nodes: bool = False) -> list:
    """
    Fetch all edges from the database involving the given UniProt codes.
    
    Args:
        ppi_dataset_name: Name of the PPI dataset
        licenced: Whether to use licensed version
        uniprot_codes: Set of UniProt codes to search for
        require_both_nodes: If True, both from_protein AND to_protein must be in uniprot_codes.
                          If False, either from_protein OR to_protein can be in uniprot_codes.
    
    Returns:
        List of ProteinProteinInteraction objects
    """
    dataset_object = get_ppi_ds(ppi_dataset_name, licenced)
    if not dataset_object:
        return []
    
    if require_both_nodes:
        interaction_objects = models.ProteinProteinInteraction.objects.filter(
            Q(ppi_dataset=dataset_object) &
            Q(from_protein__uniprot_code__in=uniprot_codes) &
            Q(to_protein__uniprot_code__in=uniprot_codes)
        )
    else:
        interaction_objects = models.ProteinProteinInteraction.objects.filter(
            Q(ppi_dataset=dataset_object) &
            (Q(from_protein__uniprot_code__in=uniprot_codes) | Q(to_protein__uniprot_code__in=uniprot_codes))
        )
    
    return list(interaction_objects)

def map_edges(ppi_dataset, edges, nodes_mapped_dict, drugstone_mapping, drugstone_identifier="drugstone_id"):
    protein_to_node = {}
    for key, node in nodes_mapped_dict.items():
        for protein in node.get(drugstone_identifier) or []:
            protein_to_node.setdefault(protein, key)

    mapped_edges, other_edges, pairs = [], [], []
    for edge in edges:
        if edge["from"] in nodes_mapped_dict and edge["to"] in nodes_mapped_dict:
            mapped_edges.append(edge)
            pairs.extend(
                {"from": fr, "to": to}
                for fr in nodes_mapped_dict[edge["from"]][drugstone_identifier]
                for to in nodes_mapped_dict[edge["to"]][drugstone_identifier]
            )
        else:
            other_edges.append(edge)

    found = {}
    for edge in fetch_edges_from_input(ppi_dataset['name'], ppi_dataset['licenced'], pairs) if pairs else []:
        if "dataset" in edge:
            found.setdefault((edge["from"], edge["to"]), edge)

    is_omnipath = ppi_dataset['name'] == "OmniPath"
    result = []
    for edge in mapped_edges:
        annotated = {}
        for fr in nodes_mapped_dict[edge["from"]][drugstone_identifier]:
            for to in nodes_mapped_dict[edge["to"]][drugstone_identifier]:
                matches = [found.get((fr, to))]
                if is_omnipath or matches[0] is None:
                    matches.append(found.get((to, fr)))
                for match in matches:
                    if match is None:
                        continue
                    key = (protein_to_node[match["from"]], protein_to_node[match["to"]])
                    annotated.setdefault(key, {**match, "from": key[0], "to": key[1]})
        result.extend(annotated.values() or [edge])

    if other_edges:
        for edge in fetch_edges_from_input(ppi_dataset['name'], ppi_dataset['licenced'], other_edges):
            if drugstone_mapping and edge["from"] in drugstone_mapping and edge["to"] in drugstone_mapping:
                edge["from"] = drugstone_mapping[edge["from"]]
                edge["to"] = drugstone_mapping[edge["to"]]
            result.append(edge)
    return result
