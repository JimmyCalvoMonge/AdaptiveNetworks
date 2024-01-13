import networkx as nx
import random
import numpy as np
import pandas as pd
import multiprocessing
import functools
import warnings
warnings.filterwarnings("ignore")
from SIR_adaptive_net_ import get_neighbors, vertxmaxer, infneighs

n_cores = min(multiprocessing.cpu_count() - 2, 30)
print(f'Using {n_cores} cores')


def infneighsLayer1(ntw, node, ns, ilist, IDtoSecondLayer):
    inf = [k for k in range(len(ilist)) if ilist[k] == 1]
    neighlist = get_neighbors(ntw, node, ns)
    infneigh = [neigh for neigh in neighlist if neigh in inf and neigh not in IDtoSecondLayer[node]]
    return [node, len(infneigh), infneigh]


def infneighsLayer2(ntw, node, ns, ilist, IDtoSecondLayer):
    inf = [k for k in range(len(ilist)) if ilist[k] == 1]
    neighlist = IDtoSecondLayer[node][1]
    infneigh = [neigh for neigh in neighlist if neigh in inf]
    return [node, len(infneigh), infneigh]


def sporadicRun(idx, n, T, v, T2, v2, pi, pr, ns, net, max_layer_2, it):

    i0 = random.randint(0, n-1)
    i = [1 if k == i0 else 0 for k in range(n)]
    s = [1 - k for k in i]
    r = [0]*n

    if not isinstance(v, list):
        v = [v]*n

    if not isinstance(T, list):
        T = [T]*n

    if not isinstance(v2, list):
        v2 = [v2]*n

    if not isinstance(T2, list):
        T2 = [T2]*n

    rednet = net.copy()
    # print(f'Edges initial: {rednet.number_of_edges()}')

    IDtoSecondLayer = [[node, 
                        random.sample(
                            [k for k in range(n) if k not in list(rednet.neighbors(node))],
                            int(max_layer_2/2))] 
                        for node in range(n)]
    
    for node in IDtoSecondLayer:
        if s[node[0]] == 1:
            for new_conn in node[1]:
                rednet.add_edge(node[0], new_conn)

    net_contact_dict = {
        node: len([neig for neig in rednet.neighbors(node)]) for node in range(n)
    }

    bepidemic_ = pd.DataFrame({})
    contacts_df_ = pd.DataFrame({})

    for t in range(it):

        # print(f'day = {t}')
        
        # Use edges variable
        A = nx.adjacency_matrix(rednet).todense()
        infn = np.matmul(A, i)
        infn = infn[0, :].tolist()[0]
        ninfp = [1-(1-pi)**nodes for nodes in infn]

        # Newly infected nodes
        x = [random.random() for _ in range(n)]
        newinf = [x[k] <= ninfp[k] for k in range(n)]
        upi = [1 if (s[k] == 1) and (newinf[k] == True) else 0 for k in range(n)]

        # Update S and I vectors: after infections 
        inew = [i[k] + upi[k] for k in range(n)]
        snew = [s[k] - upi[k] for k in range(n)]

        y = [random.random() for _ in range(n)]
        newrec = [y[k] <= pr for k in range(n)]

        # Update I and R vectors: after recovery 
        upr = [1 if (inew[k] == 1) and (newrec[k] == True) else 0 for k in range(n)]
        rnew = [r[k] + upr[k] for k in range(n)]
        inew = [inew[k] - upr[k] for k in range(n)]
        
        # Update variable vectors
        s = snew
        i = inew
        r = rnew

        # ===== The maximization process =====
        rednet = net.copy() # Define a dummy network: this is the network to play with

        # Get the info of all susceptible nodes at time t
        sus = [k for k in range(n) if s[k] == 1]

        # List of infected neighbors for each susceptible node
        # Filter only those that are not in the second layer
        infofsus = [infneighsLayer1(net, node, 1, inew, IDtoSecondLayer) for node in sus]

        # Get only the info of susceptible with  infected neighbours
        susinfs = [k for k in infofsus if k[1]!=0]

        # Layer 1: remove edges from contacts, as long as they are not in Layer 2
        # Create the list of all the edges to remove at time t: of the form {s, i-rem}
        redneighslay1 = [vertxmaxer(net, susinfs[k], [pi,pr], [v[susinfs[k][0]], T[susinfs[k][0]]], ns)[0] for k in range(len(susinfs))]

        list_ = [[susinfs[k][0], random.sample(susinfs[k][2], susinfs[k][1] - redneighslay1[k])] for k in range(len(redneighslay1))]
    
        for edges_list in list_:
            flat_edges = [(edges_list[0], edges_list[1][k]) for k in range(len(edges_list[1]))]
            try:
                rednet.remove_edges_from(flat_edges)
            except Exception:
                print('Error!')
                pass

        # Layer 2: we add sporadic contacts in this layer.
        # print(f'Edges before: {rednet.number_of_edges()}')
        # 2.1 Remove previous second layer connections:
        for node in IDtoSecondLayer:
            rednet.remove_edges_from([(node[0], new_conn) for new_conn in node[1]])
        # print(f'Edges after: {rednet.number_of_edges()}')
        
        # 2.2 List of infected layer 2 contacts for each susceptible node
        infofsus2 = [infneighsLayer2(net, node, 1, inew, IDtoSecondLayer) for node in sus]

        # 2.3 Filter only those that are in the second layer
        susinfs2 = [k for k in infofsus2 if k[1]!=0]

        # 2.4 Create the list of all new edges to add for sporadic layer at time t
        redneighslay2 = [vertxmaxer(net, [susinfs2[k][0], int(max_layer_2/2) - susinfs2[k][1]],
                                    [pi,pr], [v2[susinfs2[k][0]], T2[susinfs2[k][0]]], ns)[0] for k in range(len(susinfs2))]

        # For each susceptible node, if node is affected by adaptive scenario, then down resample Layer 2 contacts.
        # If node is not affected, continue to sample contacts in 2nd Layer as before.

        susnodes2 = [susinfs2[j][0] for j in range(len(susinfs2))]

        for node in IDtoSecondLayer:
            
            if s[node[0]] == 1:

                # Population to sample for sporadic contacts: all other nodes, which are not neighbors of this node.
                pop_to_sample = [k for k in range(n) if k not in list(rednet.neighbors(node[0]))]

                if node[0] in susnodes2:

                    # Sample as many nodes as the vertmaxer said, for this node.
                    j = susnodes2.index(node[0])
                    # print('For second layer maximizer said', redneighslay2[j], 'contacts')
                    list_2_ = random.sample(pop_to_sample, redneighslay2[j])
                    flat_edges = [(node[0], li) for li in list_2_]

                else:

                    # Just sample max_layer_2/2 nodes randomly from 2nd Layer.
                    list_2_ = random.sample(pop_to_sample, int(max_layer_2/2))
                    flat_edges = [(node[0], li) for li in list_2_]
                
                rednet.add_edges_from(flat_edges) # Add these edges, from the node to sporadic contacts
                IDtoSecondLayer[node[0]] = [node[0], list_2_] # Change IDtoSecondLayer

            else:

                if r[node[0]] == 1:

                    list_2_ = random.sample(pop_to_sample, int(max_layer_2/2))
                    flat_edges = [(node[0], li) for li in list_2_]
                    rednet.add_edges_from(flat_edges) # Add these edges, from the node to sporadic contacts
                    IDtoSecondLayer[node[0]] = [node[0], list_2_] # Change IDtoSecondLayer

                else:
                    IDtoSecondLayer[node[0]] = [node[0], []]
                # Other considerations: recovered nodes can have sporadic contacts. How about sporadic contacts for infected?

        # print(f'Edges final: {rednet.number_of_edges()}')
        # === Update done, compute day values === #

        sus_contacts = {}
        for node in sus:
            if net_contact_dict[node] != 0:
                sus_contacts[node] = len([neig for neig in rednet.neighbors(node) 
                                          if neig not in IDtoSecondLayer[node][1]])/net_contact_dict[node]
            else:
                sus_contacts[node] = 0 

        sus_contacts_ratios = np.nanmean(list(sus_contacts.values()))

        bepidemic_ = pd.concat([bepidemic_, pd.DataFrame({
            'day': [t+1],
            's': [sum(snew)],
            'i': [sum(inew)],
            'r': [sum(rnew)],
            'edgecount': [rednet.number_of_edges()],
            'suscedgecount': sus_contacts_ratios
        })])

        contacts_df_day = pd.DataFrame({
            'node': list(range(n)),
            'L2_contacts': [len(nd[1]) for nd in IDtoSecondLayer],
            'L1_contacts': [len([ng for ng in list(rednet.neighbors(node)) 
                                            if ng not in IDtoSecondLayer[node][1]]) 
                                            for node in range(n)]
        })
        contacts_df_day['day'] = t
        contacts_df_ = pd.concat([contacts_df_, contacts_df_day], ignore_index=False)

    bepidemic_['edgecount'] = bepidemic_['edgecount']/np.nanmax(bepidemic_['edgecount'])
    bepidemic_['suscedgecount'] = bepidemic_['suscedgecount']/np.nanmax(bepidemic_['suscedgecount'])

    bepidemic_['s'] = bepidemic_['s']/n
    bepidemic_['i'] = bepidemic_['i']/n
    bepidemic_['r'] = bepidemic_['r']/n
    
    # bepidemic_.to_csv('./Data/layerexperiment.csv', index=False)
    # contacts_df_.to_csv('./Data/contactinfo.csv', index=False)
    # print('End of experiment')

    bepidemic_['index'] = idx
    contacts_df_['index'] = idx

    return bepidemic_, contacts_df_


def sporadicExperiment():

    bepis = 50

    n = 500 # Total population
    it = 200 # Days
    
    T = 14 # Planning horizon Layer 1
    v = 0.05 # Risk perception Layer 1

    T2 = 14 # Planning horizon Layer 2
    v2 = 0.025 # Risk perception Layer 2

    pi = 0.02 # Probability of infection
    pr = 0.03 # Probability of recovery
    ns = 1
    net = nx.gnp_random_graph(n, 0.025) # Underlying network

    max_layer_2 = 20

    pool = multiprocessing.Pool(n_cores)
    bepidist0 = pool.map(functools.partial(sporadicRun, T=T, v=v, 
                                           T2=T2, v2=v2, 
                                           pi=pi, pr=pr, 
                                           ns=ns, net=net, 
                                           max_layer_2=max_layer_2, it=it), range(bepis))
    pool.close()
    pool.join()

    bepidist_ = pd.concat([bep[0] for bep in bepidist0], ignore_index=True)
    bepidist_ = bepidist_.groupby(['day'], as_index=False).agg({
        's': ['mean', 'max', 'min', 'std'],
        'i': ['mean', 'max', 'min', 'std'],
        'r': ['mean', 'max', 'min', 'std'],
        'edgecount': ['mean', 'max', 'min', 'std'],
        'suscedgecount': ['mean', 'max', 'min', 'std'],})
    bepidist_.columns = bepidist_.columns.droplevel()
    bepidist_.columns = ['day', 
                        's_mean', 's_max', 's_min', 's_std', 
                        'i_mean', 'i_max', 'i_min', 'i_std',
                        'r_mean', 'r_max', 'r_min', 'r_std', 
                        'edgecount_mean', 'edgecount_max', 'edgecount_min', 'edgecount_std', 
                        'suscedgecount_mean', 'suscedgecount_max', 'suscedgecount_min', 'suscedgecount_std']
    bepidist_.to_csv('./Data/sporadicexperiment.csv', index=False)

    contacts_df_ = pd.concat([bep[1] for bep in bepidist0], ignore_index=True)
    contacts_df_ = contacts_df_.groupby(['day'], as_index=False).agg({
        'L1_contacts': ['mean', 'max', 'min', 'std'],
        'L2_contacts': ['mean', 'max', 'min', 'std'],})
    contacts_df_.columns = contacts_df_.columns.droplevel()
    contacts_df_.columns = ['day', 
                    'L1_mean', 'L1_max', 'L1_min', 'L1_std', 
                    'L2_mean', 'L2_max', 'L2_min', 'L2_std']
    contacts_df_.to_csv('./Data/sporadicexperimentcontactinfo.csv', index=False)

    print("Sporadic contacts experiment done! :)")


def twoNetworksRun(idx, n, ns, v, T, v2, T2, pi, pr, it):

    net1 = nx.gnp_random_graph(n, 0.05) # Underlying network 1
    net2 = nx.gnp_random_graph(n, 0.02) # Underlying network 2

    i0 = random.randint(0, n-1)
    i = [1 if k == i0 else 0 for k in range(n)]
    s = [1 - k for k in i]
    r = [0]*n

    if not isinstance(v, list):
        v = [v]*n

    if not isinstance(T, list):
        T = [T]*n

    if not isinstance(v2, list):
        v2 = [v2]*n

    if not isinstance(T2, list):
        T2 = [T2]*n

    rednet1 = net1.copy()
    rednet2 = net2.copy()
    
    net_contact_dict1 = {
        node: len([neig for neig in rednet1.neighbors(node)]) for node in range(n)
    }
    net_contact_dict2 = {
        node: len([neig for neig in rednet2.neighbors(node)]) for node in range(n)
    }
    bepidemic_ = pd.DataFrame({})

    for t in range(it):

        # print(f'day = {t}')
        
        # Use edges variable
        A1 = nx.adjacency_matrix(rednet1).todense()
        A2 = nx.adjacency_matrix(rednet2).todense()
        A = np.logical_or(A1,A2).astype(int)

        infn = np.matmul(A, i)
        infn = infn[0, :].tolist()[0]
        ninfp = [1-(1-pi)**nodes for nodes in infn]

        # Newly infected nodes
        x = [random.random() for _ in range(n)]
        newinf = [x[k] <= ninfp[k] for k in range(n)]
        upi = [1 if (s[k] == 1) and (newinf[k] == True) else 0 for k in range(n)]

        # Update S and I vectors: after infections 
        inew = [i[k] + upi[k] for k in range(n)]
        snew = [s[k] - upi[k] for k in range(n)]

        y = [random.random() for _ in range(n)]
        newrec = [y[k] <= pr for k in range(n)]

        # Update I and R vectors: after recovery 
        upr = [1 if (inew[k] == 1) and (newrec[k] == True) else 0 for k in range(n)]
        rnew = [r[k] + upr[k] for k in range(n)]
        inew = [inew[k] - upr[k] for k in range(n)]
        
        # Update variable vectors
        s = snew
        i = inew
        r = rnew

        # ===== The maximization process =====
        rednet1 = net1.copy() # Define a dummy network: this is the network to play with
        rednet2 = net2.copy() # Define a dummy network: this is the network to play with

        # Get the info of all susceptible nodes at time t
        sus = [k for k in range(n) if s[k] == 1]

        # List of infected neighbors for each susceptible node
        # Filter only those that are not in the second layer
        infofsus1 = [infneighs(net1, node, 1, inew) for node in sus]

        # Get only the info of susceptible with  infected neighbours
        susinfs1 = [k for k in infofsus1 if k[1]!=0]

        # Layer 1: remove edges from contacts
        # Create the list of all the edges to remove at time t: of the form {s, i-rem}
        redneighs1 = [vertxmaxer(net1, susinfs1[k], [pi,pr], [v[susinfs1[k][0]], T[susinfs1[k][0]]], ns)[0] for k in range(len(susinfs1))]

        list1_ = [[susinfs1[k][0], random.sample(susinfs1[k][2], susinfs1[k][1] - redneighs1[k])] for k in range(len(redneighs1))]
    
        for edges_list in list1_:
            flat_edges = [(edges_list[0], edges_list[1][k]) for k in range(len(edges_list[1]))]
            try:
                rednet1.remove_edges_from(flat_edges)
            except Exception:
                print('Error!')
                pass

        # List of infected neighbors for each susceptible node
        # Filter only those that are not in the second layer
        infofsus2 = [infneighs(net2, node, 1, inew) for node in sus]

        # Get only the info of susceptible with  infected neighbours
        susinfs2 = [k for k in infofsus2 if k[1]!=0]

        # Layer 2: remove edges from contacts
        # Create the list of all the edges to remove at time t: of the form {s, i-rem}
        redneighs2 = [vertxmaxer(net2, susinfs2[k], [pi,pr], [v2[susinfs2[k][0]], T2[susinfs2[k][0]]], ns)[0] for k in range(len(susinfs2))]

        list2_ = [[susinfs2[k][0], random.sample(susinfs2[k][2], susinfs2[k][1] - redneighs2[k])] for k in range(len(redneighs2))]
    
        for edges_list in list2_:
            flat_edges = [(edges_list[0], edges_list[1][k]) for k in range(len(edges_list[1]))]
            try:
                rednet2.remove_edges_from(flat_edges)
            except Exception:
                print('Error!')
                pass

        sus_contacts1 = {}
        for node in sus:
            if net_contact_dict1[node] != 0:
                sus_contacts1[node] = len([neig for neig in rednet1.neighbors(node)])/net_contact_dict1[node]
            else:
                sus_contacts1[node] = 0
        sus_contacts_ratios1 = np.nanmean(list(sus_contacts1.values()))

        sus_contacts2 = {}
        for node in sus:
            if net_contact_dict2[node] != 0:
                sus_contacts2[node] = len([neig for neig in rednet2.neighbors(node)])/net_contact_dict2[node]
            else:
                sus_contacts2[node] = 0
        sus_contacts_ratios2 = np.nanmean(list(sus_contacts2.values()))

        bepidemic_ = pd.concat([bepidemic_, pd.DataFrame({
            'day': [t+1],
            's': [sum(snew)],
            'i': [sum(inew)],
            'r': [sum(rnew)],
            'edgecount1': [rednet1.number_of_edges()],
            'suscedgecount1': sus_contacts_ratios1,
            'edgecount2': [rednet2.number_of_edges()],
            'suscedgecount2': sus_contacts_ratios2,
        })])

    bepidemic_['s'] = bepidemic_['s']/n
    bepidemic_['i'] = bepidemic_['i']/n
    bepidemic_['r'] = bepidemic_['r']/n

    bepidemic_['edgecount1'] = bepidemic_['edgecount1']/np.nanmax(bepidemic_['edgecount1'])
    bepidemic_['suscedgecount1'] = bepidemic_['suscedgecount1']/np.nanmax(bepidemic_['suscedgecount1'])

    bepidemic_['edgecount2'] = bepidemic_['edgecount2']/np.nanmax(bepidemic_['edgecount2'])
    bepidemic_['suscedgecount2'] = bepidemic_['suscedgecount2']/np.nanmax(bepidemic_['suscedgecount2'])

    # bepidemic_.to_csv('./Data/2layerexperiment.csv', index=False)
    # print('End of experiment')

    return bepidemic_


def twoNetworksExperiment():

    bepis = 50

    n = 500 # Total population
    it = 200 # Days
    
    T = 14 # Planning horizon Layer 1
    v = 0.025 # Risk perception Layer 1

    T2 = 14 # Planning horizon Layer 2
    v2 = 0.0125 # Risk perception Layer 2

    pi = 0.009 # Probability of infection
    pr = 0.02 # Probability of recovery
    ns = 1

    pool = multiprocessing.Pool(n_cores)
    bepidist0 = pool.map(functools.partial(twoNetworksRun, n=n, ns=ns, 
                                           v=v, T=T, v2=v2, T2=T2, pi=pi, pr=pr, it=it), range(bepis))
    pool.close()
    pool.join()

    bepidist_ = pd.concat([bep for bep in bepidist0], ignore_index=True)

    bepidist_ = bepidist_.groupby(['day'], as_index=False).agg({
        's': ['mean', 'max', 'min', 'std'],
        'i': ['mean', 'max', 'min', 'std'],
        'r': ['mean', 'max', 'min', 'std'],
        'edgecount1': ['mean', 'max', 'min', 'std'],
        'suscedgecount1': ['mean', 'max', 'min', 'std'],
        'edgecount2': ['mean', 'max', 'min', 'std'],
        'suscedgecount2': ['mean', 'max', 'min', 'std']})
    bepidist_.columns = bepidist_.columns.droplevel()
    bepidist_.columns = ['day', 
                        's_mean', 's_max', 's_min', 's_std', 
                        'i_mean', 'i_max', 'i_min', 'i_std',
                        'r_mean', 'r_max', 'r_min', 'r_std', 
                        'edgecount1_mean', 'edgecount1_max', 'edgecount1_min', 'edgecount1_std', 
                        'suscedgecount1_mean', 'suscedgecount1_max', 'suscedgecount1_min', 'suscedgecount1_std',
                        'edgecount2_mean', 'edgecount2_max', 'edgecount2_min', 'edgecount2_std', 
                        'suscedgecount2_mean', 'suscedgecount2_max', 'suscedgecount2_min', 'suscedgecount2_std']

    bepidist_.to_csv('./Data/twonetworkexperiment.csv', index=False)

    print("Two network experiment done! :)")


if __name__ == '__main__':

    print("Start")

    sporadicExperiment()
    
    twoNetworksExperiment()

    print("End")
