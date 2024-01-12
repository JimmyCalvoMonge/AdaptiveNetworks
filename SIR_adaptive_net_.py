"""
Adaptive human behavior in epidemics on networks

Baltazar Espinoza: Mathematica Version
Jimmy Calvo: Python Version

November 2023
"""

import networkx as nx
import random
import numpy as np
import pandas as pd
import multiprocessing
import functools
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings("ignore")
import itertools
from scipy.stats import norm

n_cores = min(multiprocessing.cpu_count() - 2, 30)
print(f'Using {n_cores} cores')

# Get neighbors at levels (first, second, third level neighbors)
def get_neighbors(ntw, node, levels):
    if levels > 1:
        subgraph = nx.ego_graph(ntw, node, radius=levels)
        return list(subgraph.nodes())
    else:
        return list(ntw.neighbors(node))


def infneighs(ntw, node, ns, ilist):
    inf = [k for k in range(len(ilist)) if ilist[k] == 1]
    neighlist = get_neighbors(ntw, node, ns)
    infneigh = [neigh for neigh in neighlist if neigh in inf]
    return [node, len(infneigh), infneigh]

# Single peak function
def fs(x,y,z):
    return (x*z - z**2)**y

# The maxer code runs the optimization process for each vertex 
def vertxmaxer(ntwk, nodeid, dispars, u, neis):

    ntw = ntwk # network used 
    node = nodeid[0] # the node for which the maximizer is running and the information about inf neghbours
    maxed = nodeid[1]
    beta = dispars[0]
    gamma = dispars[1]
    nu = u[0] # utility parameters
    T = u[1] # time horizon
    ns = neis # Neighboorhood size to get the local prevalence

    b = 2*maxed # maximum possible # of edges for the utility function, the optimal utility is assumed at half of this
    delta = 0.99986 # discount factor

    # Vector storing the optimal s utility 
    EUS = [0]*(T+1)
    # The last step in s assumes to use all edges 
    EUS[-1] = fs(b, nu, maxed)

    # Vector storaging the optimal i utility 
    EUI = [0]*(T+1)

    # Vector storaging the optimal i utility 
    EUR = [0]*(T+1)

    # The last step in r assumes to use all edges 
    # EUR[-1] = fs(b, nu, maxed)

    # Vector storing the optimal # of edges for the susceptible node during the planning horizon 
    nedgest = [0]*(T+1)
    nedgest[-1] = maxed

    # Pt is the probability of disease transmission along an edge during a single time step
    Pt = beta
    Pir = gamma
    
    for t in range(T-1, -1, -1):

        # individual's utility starts at zero
        maxu = 0

        # Evaluating all the contact rates for susceptible individuals, all possible variable states
        for edgesd in range(maxed, 0, -1):

            Psi = 1 - (1 - Pt)**edgesd

            # Evaluating the expected utilities
            # Utility gained by being currently susceptible
            EUS[t] = fs(b, nu, edgesd) + delta*( (1-Psi)*EUS[t+1] + Psi*EUI[t+1] )

            EUI[t] = fs(b, nu, maxed) + delta*( (1-Pir)*EUI[t+1] + Pir*EUR[t+1] )

            EUR[t] = fs(b, nu, maxed) + delta*EUR[t+1]

            if EUS[t] > maxu:
                nedgest[t] = edgesd
                maxu = EUS[t]

    return nedgest


def epidemic(index, net, it, pi, pr, n):

    epidemic_ = pd.DataFrame({})

    i0 = random.randint(0, n-1)
    i = [1 if k == i0 else 0 for k in range(n)]
    s = [1 - k for k in i]
    r = [0]*n

    for t in range(it):

        A = nx.adjacency_matrix(net).todense()
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

        epidemic_ = pd.concat([epidemic_, pd.DataFrame({
            'day': [t+1],
            's': [sum(snew)],
            'i': [sum(inew)],
            'r': [sum(rnew)]
        })])
    
    return epidemic_


def episim(ntwk, epidemics, iterations, dispars, n):

    net = ntwk
    epis = epidemics
    it = iterations
    pi = dispars[0]
    pr = dispars[1] # disease parameters

    pool = multiprocessing.Pool(n_cores)
    epidist = pool.map(functools.partial(epidemic, net=net, it=it, pi=pi, pr=pr,n=n), range(epis))
    pool.close()
    pool.join()
    
    epidist = pd.concat(epidist, ignore_index=True)
    epidist = epidist.groupby(['day'], as_index=False).agg(
                      {'s':['mean', 'max', 'min', 'std'],
                       'i':['mean', 'max', 'min', 'std'],
                       'r':['mean', 'max', 'min', 'std']})
    epidist.columns = epidist.columns.droplevel()
    epidist.columns = ['day', 
                       's_mean', 's_min', 's_max', 's_std',
                       'i_mean', 'i_min', 'i_max', 'i_std',
                       'r_mean', 'r_min', 'r_max', 'r_std']

    return epidist


def bepidemic(index, net, it, pi, pr, v, T, ns, n):

    bepidemic_ = pd.DataFrame({})
    bepidemic_edge_count = pd.DataFrame({'index': [index]*n, 'node': list(range(n))})

    i0 = random.randint(0, n-1)
    i = [1 if k == i0 else 0 for k in range(n)]
    s = [1 - k for k in i]
    r = [0]*n

    rednet = net.copy()

    net_contact_dict = {
        node: len([neig for neig in rednet.neighbors(node)]) for node in range(n)
    }

    if not isinstance(v, list):
        v = [v]*n

    if not isinstance(T, list):
        T = [T]*n

    for t in range(it):

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

        # Get the info of all susceptible nodes at time t*)
        # sus=Flatten@Position[s,1];
        # infofsus=infneighs[net,#,ns,i]&/@sus;
        sus = [k for k in range(n) if s[k] == 1]
        infofsus = [infneighs(net, node, ns, inew) for node in sus]

        # (*Get only the info of susceptible with  infected neighbours*)
        # susinfs=Table[If[infofsus[[k,2]]!=0,infofsus[[k]]],{k,1,Length@sus}]/.Null->Sequence[];
        susinfs = [k for k in infofsus if k[1]!=0]

        # (* Use the maximizer to define how many edges to drop: it returns the number of edges to leave *)
        # redneighs=Table[vertxmaxer[net,susinfs[[k]],{pi,pr},{vuse[[k]],Tuse[[k]]},ns][[1]],{k,1,Length[susinfs]}];
        redneighs = [vertxmaxer(net, susinfs[k], [pi,pr], [v[susinfs[k][0]],T[susinfs[k][0]]], ns)[0]  for k in range(len(susinfs))]

        # (*Get the list of the form {susceptible  id,{list of infected neighbours to remove}}*)
        # list=Table[{susinfs[[k,1]],RandomSample[susinfs[[k,3]],susinfs[[k,2]]-redneighs[[k]]]},{k,1,Length[redneighs]}];
        list_ = [[susinfs[k][0], random.sample(susinfs[k][2], susinfs[k][1] - redneighs[k])] for k in range(len(redneighs))]

        # (*Create the list of all the edges to remove at time t: of the form {s,i-rem}*)
        # edgestorem=Flatten[
        # Table[
        # {list[[k,1]]\[UndirectedEdge]#}&/@list[[k,2]]
        # ,{k,1,Length@redneighs}]];

        for edges_list in list_:
            flat_edges = [(edges_list[0], edges_list[1][k]) for k in range(len(edges_list[1]))]
            try:
                rednet.remove_edges_from(flat_edges)
            except Exception:
                print('Error!')
                pass

        # (*Drop edges from the network*)
        # rednet = EdgeDelete[rednet,edgestorem];

        # {Count[s,1],Count[i,1],Count[r,1],EdgeCount[rednet],
        # N@Mean[Length[IncidenceList[rednet,#]]&/@sus]}

        sus_contacts = {}
        for node in sus:
            if net_contact_dict[node] != 0:
                sus_contacts[node] = len([neig for neig in rednet.neighbors(node)])/net_contact_dict[node]
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

        bepidemic_edge_count[f'day_{t}'] = [len(list(rednet.neighbors(node))) for node in range(n)]
        bepidemic_edge_count[f's_day_{t}'] = s
        bepidemic_edge_count[f'i_day_{t}'] = i
        bepidemic_edge_count[f'r_day_{t}'] = r

    bepidemic_['edgecount'] = bepidemic_['edgecount']/np.nanmax(bepidemic_['edgecount'])
    bepidemic_['suscedgecount'] = bepidemic_['suscedgecount']/np.nanmax(bepidemic_['suscedgecount'])

    return [bepidemic_ , bepidemic_edge_count]


def bepisim(ntwk, epidemics, iterations, dispars, u, neis, n, **kwargs):

    net = ntwk
    bepis = epidemics
    it = iterations
    pi = dispars[0]
    pr = dispars[1] # disease parameters
    v = u[0] # utility parameters
    T = u[1] # time horizon
    ns = neis # Neighboorhood size to get the local prevalence

    if bepis < n_cores:
        bepidist0 = [bepidemic(i, net=net, it=it,
                               pi=pi, pr=pr, v=v,
                               T=T, ns=ns, n=n) for i in range(bepis)]
    else:
        pool = multiprocessing.Pool(n_cores)
        bepidist0 = pool.map(functools.partial(bepidemic, net=net, it=it,
                                            pi=pi, pr=pr, v=v,
                                            T=T, ns=ns), range(bepis))
        pool.close()
        pool.join()

    bepidist = pd.concat([bep[0] for bep in bepidist0], ignore_index=True)

    bepidist = bepidist.groupby(['day'], as_index=False).agg(
                      {'s':['mean', 'max', 'min', 'std'],
                       'i':['mean', 'max', 'min', 'std'],
                       'r':['mean', 'max', 'min', 'std'],
                       'edgecount':['mean', 'max', 'min', 'std'],
                       'suscedgecount':['mean', 'max', 'min', 'std']})
    
    bepidist.columns = bepidist.columns.droplevel()
    bepidist.columns = ['day', 
                        's_mean', 's_min', 's_max', 's_std', 
                        'i_mean', 'i_min', 'i_max', 'i_std', 
                        'r_mean', 'r_min', 'r_max', 'r_std',
                        'edgecount_mean', 'edgecount_min', 'edgecount_max', 'edgecount_std', 
                        'suscedgecount_mean', 'suscedgecount_min', 'suscedgecount_max', 'suscedgecount_std']

    # Get contact information for each individual at each day. For all simulations.
    get_node_history = kwargs.get('get_node_history', False)
    if get_node_history:
        bepidist_contacts_ = pd.concat([bep[1] for bep in bepidist0], ignore_index=True)
        return bepidist, bepidist_contacts_

    return bepidist

# Figures 

def infected_comparison_fig(epidist, bepidist, fig, v, T, n, **kwargs):

    mean = (epidist['i_mean']/n).tolist()
    lower = (epidist['i_min']/n).tolist()
    upper = (epidist['i_max']/n).tolist()

    meanb = (bepidist['i_mean']/n).tolist()
    lowerb = (bepidist['i_min']/n).tolist()
    upperb = (bepidist['i_max']/n).tolist()

    axarr = fig.add_subplot(1,1,1)

    plt.plot(mean,'-b', label="Mean classic model")
    # fill the area with black color, opacity 0.15
    plt.fill_between(list(range(len(mean))), upper, lower, color="b", alpha=0.2)

    plt.plot(meanb,'-r', label="Mean adaptive model")
    # fill the area with black color, opacity 0.15
    plt.fill_between(list(range(len(meanb))), upperb, lowerb, color="r", alpha=0.2)
    
    plt.title(f"Infected proportion comparison (v={v}, T={T})")
    plt.legend(loc="upper right")

    if 'save' in kwargs:
        idx = kwargs.get('idx', 1)
        plt.savefig(f'./Figures/infected_comparison_{idx}.png')

    return fig


def edge_reduction_comparison_fig(bepidist, fig, v, T, **kwargs):

    mean_all = bepidist['edgecount_mean'].tolist()
    lower_all = bepidist['edgecount_min'].tolist()
    upper_all = bepidist['edgecount_max'].tolist()

    mean_sus = bepidist['suscedgecount_mean'].tolist()
    lower_sus = bepidist['suscedgecount_min'].tolist()
    upper_sus = bepidist['suscedgecount_max'].tolist()

    axarr = fig.add_subplot(1,1,1)

    plt.plot(mean_all[0:50],'-k', label="Network avg. edges reduction")
    plt.fill_between(list(range(len(mean_all[0:50]))), upper_all[0:50], lower_all[0:50], color="k", alpha=0.2)
    plt.plot(mean_sus[0:50], linestyle='--', color='b', label="Individuals avg. edges reduction")
    plt.fill_between(list(range(len(mean_sus[0:50]))), upper_sus[0:50], lower_sus[0:50], color="b", alpha=0.2)
    plt.title(f"Global and local behavioral responses (v={v}, T={T})")
    plt.legend(loc="lower right")

    if 'save' in kwargs:
        idx = kwargs.get('idx', 1)
        plt.savefig(f'./Figures/local_global_behavior_comparison_{idx}.png')

    return fig


def infected_edge_reduction_fig(bepidist, fig, v, T, n,**kwargs):

    meanb = (bepidist['i_mean']/n).tolist()
    lowerb = (bepidist['i_min']/n).tolist()
    upperb = (bepidist['i_max']/n).tolist()

    meanbs = (bepidist['s_mean']/n).tolist()
    lowerbs = (bepidist['s_min']/n).tolist()
    upperbs = (bepidist['s_max']/n).tolist()

    mean_sus = bepidist['suscedgecount_mean'].tolist()
    lower_sus = bepidist['suscedgecount_min'].tolist()
    upper_sus = bepidist['suscedgecount_max'].tolist()

    axarr = fig.add_subplot(1,1,1)

    plt.plot(mean_sus, linestyle='--', color='b', label="Individuals avg. edges reduction")
    plt.fill_between(list(range(len(mean_sus))), upper_sus, lower_sus, color="b", alpha=0.2)
    
    plt.plot(meanb,'-r', label="Mean adaptive model (Infected)")
    # fill the area with black color, opacity 0.15
    plt.fill_between(list(range(len(meanb))), upperb, lowerb, color="r", alpha=0.2)

    plt.plot(meanbs,'-g', label="Mean adaptive model")
    # fill the area with black color, opacity 0.15
    plt.fill_between(list(range(len(meanbs))), upperbs, lowerbs, color="g", alpha=0.2)

    plt.title(f"Susceptible/Infected proportions and individual effort (v={v}, T={T})")
    plt.legend(loc="lower right")

    if 'save' in kwargs:
        idx = kwargs.get('idx', 1)
        plt.savefig(f'./Figures/infected_and_effort_{idx}.png')

    return fig


def get_heatmap_data(n, ped):

    vs = [round(vv, 4) for vv in np.linspace(0.01, 0.1, 10)]
    Ts = list(range(7, 70, 7))
    combs = list(itertools.product(vs, Ts))

    epidist_all = pd.DataFrame({})
    bepidist_all = pd.DataFrame({})

    for comb in combs:
        
        v = comb[0]
        T = comb[1]

        print(f"""
        Simulation for ({v},{T}) started: ====================>
        """)

        try:

            print("Started no behavior net")
            net = nx.gnp_random_graph(n, ped)
            epidist = episim(ntwk=net, epidemics=10, iterations=200, dispars=[0.05, 0.04])
            epidist['V'] = v
            epidist['T'] = T

            print("Started behavior net")
            net = nx.gnp_random_graph(n, ped)
            bepidist = bepisim(ntwk=net, epidemics=10, iterations=200, dispars=[0.05, 0.04], u=[v, T], neis=1)
            bepidist['V'] = v
            bepidist['T'] = T

            epidist_all = pd.concat([epidist_all, epidist], ignore_index=True)
            bepidist_all = pd.concat([bepidist_all, bepidist], ignore_index=True)

        except Exception as e:
            print(f"Error with ({v},{T}): {e}")

    epidist_all.to_csv(f'./Data/epidist.csv', index=False)
    bepidist_all.to_csv(f'./Data/bepidist.csv', index=False)


def get_all_figures(n, ped):

    vs = [round(vv, 4) for vv in np.linspace(0.01, 0.1, 10)]
    # Ts = list(range(7, 70, 7))

    # Example with plots
    T = 7
    # v = 0.05
    nei = 1

    for idx, v in enumerate(vs):

        print(v, T, nei, "------------>")

        v_name = v
        T_name = T

        idx_use = f'{idx}_neigh_{nei}'

        print("Started no behavior net")
        net = nx.gnp_random_graph(n, ped)
        epidist = episim(ntwk=net, epidemics=10,
                         iterations=200, dispars=[0.05, 0.04], n=n)

        print("Started behavior net")
        net = nx.gnp_random_graph(n, ped)
        bepidist = bepisim(ntwk=net, epidemics=10,
                           iterations=200, dispars=[0.05, 0.04], u=[v, T], neis=nei, n=n)

        print("Figure 1 ----")

        fig = plt.figure()
        fig = infected_comparison_fig(epidist, bepidist,
                                      fig, v_name, T_name, n, save=True, idx=idx_use)

        print("Figure 2 ----")

        fig = plt.figure()
        fig = edge_reduction_comparison_fig(bepidist, fig,
                                            v_name, T_name, n, save=True, idx=idx_use)
        

        print('Figure 3 ----')

        fig = plt.figure()
        fig = infected_edge_reduction_fig(bepidist, fig, 
                                          v_name, T_name, n, save=True, idx=idx_use)

# Article example with distributions and subpopulations
        
def example_with_distributions(n, ped):

    print('Example with distributions')

    pop1_n = 200
    pop2_n = 300

    pop1_T = list(np.random.normal(7, 1, pop1_n))
    pop1_T = [int(t) for t in pop1_T]
    pop2_T = list(np.random.normal(16, 1, pop2_n))
    pop2_T = [int(t) for t in pop2_T]
    T = pop1_T + pop2_T

    pop1_v = list(np.random.normal(0.05, 0.0025, pop1_n))
    pop2_v = list(np.random.normal(0.08, 0.006, pop2_n))
    v = pop1_v + pop2_v

    # We have two populations:
    # Pop 1: n = 200, centered arround 7 days horizon and 0.05 utility parameter (more risk averse)
    # Pop 2: n = 500, centered arround 16 days horizon and 0.08 utility parameter (less risk averse)

    # Example semi-homogeneous:
    # T = [7]*pop1_n + [16]*pop2_n
    # v = [0.05]*pop1_n + [0.08]*pop2_n

    distribution_df = pd.DataFrame({
        'v': v,
        'T': T
    })
    distribution_df.to_csv('distribution_example.csv', index=False)

    # Plot the distributions:
    plt.figure()
    plt.hist(T)
    plt.savefig(f'./Figures/T_dist_example_use.png')

    plt.figure()
    plt.hist(v)
    plt.savefig(f'./Figures/v_dist_example_use.png')

    idx = 'DistributionExample'
    v_name = 'norm(0.05, 0.01)'
    T_name = 'norm(7,1) + norm(16,1)'

    print("Started no behavior net")
    net = nx.gnp_random_graph(n, ped)
    epidist = episim(ntwk=net, epidemics=50,
                        iterations=200, dispars=[0.05, 0.04], n=n)

    print("Started behavior net")
    net = nx.gnp_random_graph(n, ped)
    bepidist, bepidist_contacts_ = bepisim(ntwk=net, epidemics=50,
                       iterations=200,
                       dispars=[0.05, 0.04],
                       u=[v, T],
                       neis=1, n=n, get_node_history=True)

    epidist.to_csv('./Data/Distribution_epidist.csv', index=False)
    bepidist.to_csv('./Data/Distribution_bepidist.csv', index=False)
    bepidist_contacts_.to_csv('./Data/Distribution_bepidist_contacts.csv', index=False)

    print("Figure 1 ----")
    fig = plt.figure()
    fig = infected_comparison_fig(epidist, bepidist,
                                    fig, 
                                    v_name, T_name, n, save=True, idx=idx)

    print("Figure 2 ----")
    fig = plt.figure()
    fig = edge_reduction_comparison_fig(bepidist, fig,
                                    v_name, T_name, n, save=True, idx=idx)

    print('Figure 3 ----')
    fig = plt.figure()
    fig = infected_edge_reduction_fig(bepidist, fig, 
                                    v_name, T_name, n, save=True, idx=idx)

# Bifurcation example, not needed ...

def bifurcation_example(n,ped):

    pis = [round(pinf, 4) for pinf in np.linspace(0.00001, 0.1, 30)]
    prs = [0.04] # [round(pinf, 4) for pinf in np.linspace(0.0001, 0.1, 20)]
    combs = list(itertools.product(pis, prs))
    v = 0.05
    T = 7
    nei = 1

    epidist_all = pd.DataFrame({})
    bepidist_all = pd.DataFrame({})

    for idx, comb in enumerate(combs):

        print(f' ---- {comb} ----')

        print("Started no behavior net")
        net = nx.gnp_random_graph(n, ped)
        epidist = episim(ntwk=net, epidemics=10,
                         iterations=200, dispars=[comb[0], comb[1]], n=n)

        print("Started behavior net")
        net = nx.gnp_random_graph(n, ped)
        bepidist = bepisim(ntwk=net, epidemics=10,
                           iterations=200, dispars=[comb[0],
                                                    comb[1]], u=[v, T], neis=nei, n=n)
        epidist['pi'] = comb[0]
        bepidist['pi'] = comb[0]
        epidist['pr'] = comb[1]
        bepidist['pr'] = comb[1]
        epidist_all = pd.concat([epidist_all, epidist], ignore_index=True)
        bepidist_all = pd.concat([bepidist_all, bepidist], ignore_index=True)
        

    epidist_all.to_csv(f'./Data/epidist_bifurcation.csv')
    bepidist_all.to_csv(f'./Data/bepidist_bifurcation.csv')

# ----- Examples for appendix ------------------------------------------------ #

# ====== Adaptive SIR vs Network ============================================= #

def example_network_and_nodes(n,ped):

    print('Example with network and nodes')

    v = 0.05
    T = 7

    print("Started behavior net")
    net = nx.gnp_random_graph(n, ped)
    bepidist, bepidist_contacts_ = bepisim(ntwk=net, epidemics=1,
                       iterations=200,
                       dispars=[0.05, 0.04],
                       u=[v, T],
                       neis=1, n=n, get_node_history=True)

    bepidist.to_csv('./Data/ExampleNetworkNodes_bepidist.csv', index=False)
    bepidist_contacts_.to_csv('./Data/ExampleNetworkNodes_contacts.csv', index=False)

# ====== Vary the network type and network connectivity parameters =========== #

def barabasi_albert_network_experiment():
    
    print("Barabasi Albert experiment ===>")
    mms = range(5,51)

    nns = [500, 1000]
    mms = range(5,51)
    sample_mms = [mms[i] for i in range(0,len(mms),10)]

    for n in nns:

        print(f'Using n = {n} ======================================>')

        day_lags_df_all = pd.DataFrame({})

        for m in mms:

            print(f' ---- Using m == {m} ----')
            net = nx.barabasi_albert_graph(n, m)
            
            pool = multiprocessing.Pool(n_cores)
            bepidist0 = pool.map(functools.partial(bepidemic, net=net, it=200,
                                                pi=0.05, pr=0.04, v=0.05,
                                                T=7, ns=1, n=n), range(50))
            pool.close()
            pool.join()
            print('========================= ')

            day_lags_ = []
            for bep in bepidist0:
                try:
                    day_lags_.append(np.nanargmin(bep[0]['suscedgecount']) - np.nanargmin(bep[0]['edgecount']))
                except Exception as e:
                    print(f'Error: {e}')
                    pass

            day_lags_df_ = pd.DataFrame({
                'day_lags': day_lags_
            })
            day_lags_df_['m'] = m
            day_lags_df_all = pd.concat([day_lags_df_all, day_lags_df_], ignore_index=True)

            if m in sample_mms:
                bepidist_ = pd.concat([bep[0] for bep in bepidist0], ignore_index=True)
                bepidist_.to_csv(f'./Data/bepidist_contacts_ntwk_barabasi_albert_{n}_{m}.csv', index=False)

        day_lags_df_all.to_csv(f'./Data/bepidist_ntwk_barabasi_albert_{n}.csv', index=False)

    print("Barabasi Albert experiment done")


def small_world_network():
    
    print("Small World experiment ===>")

    mms = range(5, 51, 5)
    peds = [round(ped, 4) for ped in np.linspace(0.001, 0.1, 10)]
    combs = list(itertools.product(mms, peds))

    sample_mms = [5, 15, 25, 35, 45]
    sample_peds = [peds[i] for i in range(0, len(peds),5)]
    sample_combs = combs

    # mms = [35, 45]
    nns = [500, 1000]

    print(f"""
    m = {list(mms)}
    peds = {peds}
    sample_mms = {sample_mms}
    sample_peds = {sample_peds}
    """)

    for n in nns:

        print(f'Using n = {n} ======================================>')
        day_lags_df_all = pd.DataFrame({})

        for comb in combs:

            print(f' ---- Using {comb} ----')
            net = nx.watts_strogatz_graph(n=n,
                                          k=comb[0],
                                          p=comb[1])

            pool = multiprocessing.Pool(n_cores)
            bepidist0 = pool.map(functools.partial(bepidemic, net=net, it=200,
                                                pi=0.05, pr=0.04, v=0.05,
                                                T=7, ns=1,n=n), range(50))
            pool.close()
            pool.join()
            print('========================= ')

            day_lags_ = []
            for bep in bepidist0:
                try:
                    day_lags_.append(np.nanargmin(bep[0]['suscedgecount']) - np.nanargmin(bep[0]['edgecount']))
                except Exception as e:
                    print(f'Error: {e}')
                    pass

            day_lags_df_ = pd.DataFrame({
                'day_lags': day_lags_
            })

            day_lags_df_['m'] = comb[0]
            day_lags_df_['ped'] = comb[1]
            day_lags_df_all = pd.concat([day_lags_df_all, day_lags_df_], ignore_index=True)

            if comb in sample_combs:
                bepidist_ = pd.concat([bep[0] for bep in bepidist0], ignore_index=True)
                bepidist_.to_csv(f'./Data/bepidist_contacts_ntwk_small_world_{n}_{comb[0]}_{comb[1]}.csv',
                                 index=False)

        day_lags_df_all.to_csv(f'./Data/bepidist_ntwk_small_world_{n}.csv', index=False)
        
    print("Small World experiment done")


def network_connectivity_experiment():

    print('Network connectivity experiment ===>')

    nns = [500, 1000]
    conns = [round(conn, 4) for conn in np.linspace(0.001, 0.1, 100)]
    sample_conns = [conns[i] for i in range(0,len(conns),10)]

    for n in nns:

        print(f'Using n = {n} ======================================>')

        day_lags_df_all = pd.DataFrame({})

        for conn in conns:

            print(f' ---- Using ped == {conn} ----')
            net = nx.gnp_random_graph(n, conn)
            
            pool = multiprocessing.Pool(n_cores)
            bepidist0 = pool.map(functools.partial(bepidemic, net=net, it=200,
                                                pi=0.05, pr=0.04, v=0.05,
                                                T=7, ns=1, n=n), range(50))
            pool.close()
            pool.join()
            print('========================= ')

            day_lags_ = []
            for bep in bepidist0:
                try:
                    day_lags_.append(np.nanargmin(bep[0]['suscedgecount']) - np.nanargmin(bep[0]['edgecount']))
                except Exception as e:
                    print(f'Error: {e}')
                    pass

            day_lags_df_ = pd.DataFrame({
                'day_lags': day_lags_
            })
            day_lags_df_['ped'] = conn
            day_lags_df_all = pd.concat([day_lags_df_all, day_lags_df_], ignore_index=True)

            if conn in sample_conns:
                bepidist_ = pd.concat([bep[0] for bep in bepidist0], ignore_index=True)
                bepidist_.to_csv(f'./Data/bepidist_contacts_ntwk_conn_erdos_renyi_{n}_{conn}.csv', index=False)

        day_lags_df_all.to_csv(f'./Data/bepidist_ntwk_conn_erdos_renyi_{n}.csv', index=False)

    print("Network connectivity experiment done")

# ====== Simple Two Layers Example =========================================== #
# TODO

if __name__ == '__main__':

    print("Starting Experiments ====>")

    # get_all_figures(500, 0.05)

    # get_heatmap_data(500, 0.05)

    # example_with_distributions(500, 0.05)

    # bifurcation_example(500, 0.05)

    # example_network_and_nodes(500, 0.05)

    # network_connectivity_experiment()

    # barabasi_albert_network_experiment()

    small_world_network()