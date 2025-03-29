import networkx as nx
from networkx.algorithms.flow import shortest_augmenting_path
import numpy as np
from tqdm import tqdm
import json

def get_metrics_network(G):

    degree_centrality = nx.degree_centrality(G)
    closeness_centrality = nx.closeness_centrality(G)

    dict_ = {
        'Average clustering coefficient' : nx.average_clustering(G),
        'Average degree centrality' : sum(degree_centrality.values())/len(degree_centrality),
        'Average closeness centrality' : sum(closeness_centrality.values())/len(closeness_centrality),
        'Node connectivity': nx.node_connectivity(G),
        # 'Average node connectivity': nx.average_node_connectivity(G), # Too slow!
    }
    return dict_

def get_parameters(n, sims, mode, **kwargs):

    all_dict_ = {
        'Average clustering coefficient' : [],
        'Average degree centrality': [],
        'Average closeness centrality' : [],
        'Node connectivity': [],
    }

    p = kwargs.get('p', 0.05)
    m = kwargs.get('m', 50)
    k = kwargs.get('k', 15)

    for _ in tqdm(range(sims)):

        if mode == 'random':
            net_name = f'Erdos-Renyi Network ER({n}, {p})'
            G = nx.gnp_random_graph(n, p)
        elif mode == 'ab':
            net_name = f'Albert-Barabasi AB({n},{m})'
            G = nx.barabasi_albert_graph(n, m)
        elif mode == 'ws':
            net_name = f'Watts-Strogatz WS({n},{m})'
            G = nx.watts_strogatz_graph(n=n, k=k, p=p)

        measures_dict = get_metrics_network(G)
        for key in measures_dict.keys():
            all_dict_[key].append(measures_dict[key])

    avg_dict_ = {}

    for key in all_dict_.keys():
        avg_dict_[key] = np.mean(all_dict_[key])

    print(f"for {net_name} these are the values of network metrics (in average):")
    print(avg_dict_)

    return avg_dict_


if __name__ == '__main__':

    n = 500
    sims = 20
    p = 0.05
    m = 20
    k = 15

    avg_dict_ = get_parameters(n, sims, mode='random', p=p) 
    with open('./Figures/result_ER.json', 'w') as fp:
        json.dump(avg_dict_, fp)
    avg_dict_ = get_parameters(n, sims, mode='ab', m=m)
    with open('./Figures/result_AB.json', 'w') as fp:
        json.dump(avg_dict_, fp)
    avg_dict_ = get_parameters(n, sims, mode='ws', p=p, k=k)
    with open('./Figures/result_WS.json', 'w') as fp:
        json.dump(avg_dict_, fp)