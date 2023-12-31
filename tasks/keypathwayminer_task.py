import base64
import datetime
import itertools
import json
import random
import string
import time
from os.path import join

import requests

from tasks.task_hook import TaskHook

from drugstone.models import Protein, EnsemblGene

# Base URL
# url = 'http://172.25.0.1:9003/keypathwayminer/requests/'
url = 'https://exbio.wzw.tum.de/keypathwayminer/requests/'
attached_to_id = ''.join(random.choices(string.ascii_uppercase + string.digits, k=32))


def send_request(sub_url, data):
    """
    Send a POST request with form-data to a given sub-URL and retrieve the JSON response.
    Throws a RuntimeError if there was an error while submitting

    :param sub_url: Sub-URL to send the POST request to
    :param data: Data dictionary that is sent via the POST request
    :return: JSON-object from the server request
    """
    request_url = join(url, sub_url)
    try:
        response = requests.post(url=request_url, data=data)
    except Exception as e:
        print(e)

    # Check if submitting the job was successful
    if response.status_code != 200:
        raise RuntimeError(f'KPM server response code was "{response.status_code}", expected "200".')

    try:
        response_json = response.json()
    except json.decoder.JSONDecodeError:
        raise RuntimeError(f'The response could not be decoded as JSON, please check the URL:\n{request_url}')

    return response_json


def kpm_task(task_hook: TaskHook):
    """
    Run KeyPathwayMiner on given proteins and parameters remotely using the RESTful API of KPM-web
    Updates status of the TaskHook by polling the KPM-web server every second
    Writes results back to the TaskHook as 'networks'.

    :param task_hook: Needs to have 'k' set as a parameter (str or int) and a list of proteins set
    :return: None
    """
    # --- Fetch and generate the datasets
    dataset_name = 'indicatorMatrix'
    indicator_matrix_string = ''
    id_space = task_hook.parameters["config"].get("identifier", "symbol")
    proteins = []
    if id_space == 'symbol':
        proteins = Protein.objects.filter(gene__in=task_hook.seeds)
    elif id_space == 'entrez':
        proteins = Protein.objects.filter(entrez__in=task_hook.seeds)
    elif id_space == 'uniprot':
        proteins = Protein.objects.filter(uniprot_code__in=task_hook.seeds)
    elif id_space == 'ensg':
        protein_ids = {ensg.protein_id for ensg in EnsemblGene.objects.filter(name__in=task_hook.seeds)}
        proteins = Protein.objects.filter(id__in=protein_ids)
    protein_backend_ids = {p.id for p in proteins}
    for protein in proteins:
        indicator_matrix_string += f'{protein.uniprot_code}\t1\n'

    datasets = [
        {
            'name': dataset_name,
            'attachedToID': attached_to_id,
            'contentBase64': base64.b64encode(indicator_matrix_string.encode('UTF-8')).decode('ascii')
        }
    ]

    datasets_data = json.dumps(datasets)

    # --- Generate KPM settings
    k_val = str(task_hook.parameters['k'])
    kpm_settings = {
        'parameters': {
            'name': f'Drugstone run on {datetime.datetime.now()}',
            'algorithm': 'Greedy',
            'strategy': 'INES',
            'removeBENs': 'true',
            'unmapped_nodes': 'Add to negative list',
            'computed_pathways': 1,
            'graphID': 27,
            'l_samePercentage': 'false',
            'samePercentage_val': 0,
            'k_values': {
                'val': k_val,
                'val_step': '1',
                'val_max': k_val,
                'use_range': 'false',
                'isPercentage': 'false'
            },
            'l_values': {
                'val': '0',
                'val_step': '1',
                'val_max': '0',
                'use_range': 'false',
                'isPercentage': 'false',
                'datasetName': dataset_name
            }
        },
        'withPerturbation': 'false',
        'perturbation': [
            {
                'technique': 'Node-swap',
                'startPercent': '5',
                'stepPercent': '1',
                'maxPercent': '15',
                'graphsPerStep': '1'
            }
        ],
        'linkType': 'OR',
        'attachedToID': attached_to_id,
        'positiveNodes': '',
        'negativeNodes': ''
    }

    kpm_settings_data = json.dumps(kpm_settings)

    # --- Submit kpm job asynchronously
    kpm_job_data = {'kpmSettings': kpm_settings_data,
                    'datasets': datasets_data}

    submit_json = send_request('submitAsync', kpm_job_data)

    # Check if the submission was correct (add check whether parameters were correct)
    if not submit_json["success"]:
        print(f'Job submission failed. Server response:\n{submit_json}')
        raise RuntimeError(f'Job submission failed. Server response:\n{submit_json}')

    # Obtain questID for getting the result
    quest_id_data = {'questID': submit_json['questID']}
    # print(submit_json["resultUrl"])  # Remove in production

    # --- Retrieve status and update task_hook every 1s
    old_progress = -1
    while True:
        # Get status of job
        status_json = send_request('runStatus', quest_id_data)

        # Check if the questID exists (should)
        if not status_json['runExists']:
            raise RuntimeError(f'Job status retrieval failed. Run does not exist:\n{status_json}')

        # Set progress only when it changed
        progress = status_json['progress']
        if old_progress != progress:
            task_hook.set_progress(progress=progress, status='')
            old_progress = progress

        # Stop and go to results
        if status_json['completed'] or status_json['cancelled']:
            break

        time.sleep(1)

    # --- Retrieve results and write back
    results_json = send_request('results', quest_id_data)

    if not results_json['success']:
        raise RuntimeError(f'Job terminated but was unsuccessful:\n{results_json}')

    graphs_json = results_json['resultGraphs']
    # Build the networks
    network = None

    # Only build networks if the result is not empty
    if graphs_json:
        for graph in graphs_json:
            # Ignore the union set
            if graph['isUnionSet']:
                continue

            # Add nodes
            nodes = []
            for node in graph['nodes']:
                nodes.append(node['name'])

            # Add edges
            edges = []
            for edge in graph['edges']:
                edges.append({'from': edge['source'], 'to': edge['target']})

            # Add nodes and edges to network
            network = {'nodes': nodes, 'edges': edges}

    # Remapping everything from UniProt Accession numbers to internal IDs
    flat_map = lambda f, xs: (y for ys in xs for y in f(ys))
    uniprote_nodes = []
    uniprote_nodes.extend(network["nodes"])
    uniprote_nodes.extend(set(flat_map(lambda l: [l['from'], l['to']], network['edges'])))

    result_nodes = Protein.objects.filter(uniprot_code__in=uniprote_nodes)
    node_map = {}
    node_map_for_edges = {}


    for node in result_nodes:
        node_map_for_edges[node.uniprot_code] = node.id
        if id_space == 'symbol':
            node_map[node.uniprot_code] = [node.gene]
        if id_space == 'entrez':
            node_map[node.uniprot_code] = [node.entrez]
        if id_space == 'uniprot':
            node_map[node.uniprot_code] = [node.uniprot_code]
        if id_space == 'ensembl':
            node_map[node.uniprot_code] = [ensg.name for ensg in EnsemblGene.objects.filter(protein_id=node.id)]



    network["nodes"] = list(flat_map(lambda uniprot: node_map[uniprot], network["nodes"]))
    drugstone_edges = []
    for uniprot_edge in network['edges']:
        from_node = f'p{node_map_for_edges[uniprot_edge["from"]]}' if uniprot_edge['from'] in node_map_for_edges else uniprot_edge['from']
        to_node = f'p{node_map_for_edges[uniprot_edge["to"]]}' if uniprot_edge['to'] in node_map_for_edges else uniprot_edge['to']
        drugstone_edges.append({"from": from_node,"to": to_node})
    network['edges']=drugstone_edges

    node_types = {node: "protein" for node in network["nodes"]}
    is_seed = {node: node in set(map(lambda p: "p" + str(p), protein_backend_ids)) for node in network["nodes"]}
    result_dict = {
        "network": network,
        "target_nodes": [node for node in network["nodes"] if node not in task_hook.seeds],
        "node_attributes": {"node_types": node_types, "is_seed": is_seed}
    }
    task_hook.set_results(results=result_dict)
