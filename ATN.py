import numpy as np
import random
import networkx as nx
import matplotlib.pyplot as plt
import scipy
import time
import pandas as pd
import csv
import os
import niche
from scipy import integrate

#Dynamics function, takes in current state (Bis) and model parameters, outputs deltas (dBis)
#Written to take in fixed parameters that are the same for all nodes (as in Romanuk 2009)
def atn(t,y,atn_args):
    
    G = atn_args[0]
    rs = atn_args[1]
    K = atn_args[2]
    x = atn_args[3]
    y_arg = atn_args[4]
    e = atn_args[5]
    q = atn_args[6]
    b0 = atn_args[7]
    #map from node id to index in state vector y (b curr)
    index_dict = atn_args[8]
    

    
    #Based on Williams and Martinez 2004 - https://link.springer.com/content/pdf/10.1140/epjb/e2004-00122-1.pdf
    b_curr = y
    dBs = np.zeros(len(b_curr))
    for i in G.nodes(): #by node id
                
        Gi = rs[i]*b_curr[index_dict[i]]*(1-b_curr[index_dict[i]]/K)       
        in_neigh_list = list(G.in_edges([i]))
        for inl in range(0,len(in_neigh_list)):
            in_neigh_list[inl] = in_neigh_list[inl][0]     
        gains_from_resources = 0
        for j in in_neigh_list:
            b_sum = 0 #sum over i's resources for Fij
            for k in in_neigh_list:
                b_sum += (b_curr[index_dict[k]]**(1+q))
            Fij = (b_curr[index_dict[j]]**(1+q))/(b_sum + (b0**(1+q)))
            curr_gain = x*y_arg*Fij*b_curr[index_dict[i]]
            gains_from_resources += curr_gain
                    
        out_neigh_list = list(G.out_edges([i]))
        for onl in range(0,len(out_neigh_list)):
            out_neigh_list[onl] = out_neigh_list[onl][1]
        
        loss_from_consumers = 0
        for j in out_neigh_list:

            b_sum = 0 #sum over j's resources for Fji           
            #get j's resources
            j_in_neigh_list = list(G.in_edges([j]))
            for inl in range(0,len(j_in_neigh_list)):
                j_in_neigh_list[inl] = j_in_neigh_list[inl][0]
            
            for k in j_in_neigh_list:
                b_sum += (b_curr[index_dict[k]]**(1+q))
            Fji = (b_curr[index_dict[i]]**(1+q))/(b_sum + (b0**(1+q)))

            curr_loss = x*y_arg*Fji*b_curr[index_dict[j]]/e
            loss_from_consumers += curr_loss

        dBi = Gi - x*b_curr[index_dict[i]] + gains_from_resources - loss_from_consumers
        dBs[index_dict[i]] = dBi
    return dBs

# Generate parameters for the dynamic model for a web based on the Romanuk 2009 defaults
# New and improved version with dictionaries and just returning scalars for fixed parameters for this version (Faster)
def get_params(G):

    ### FIXED PARAMETERS FOR ATN DYNAMIC MODEL PER ROMANUK 2009 ####
    
    # intrinsic growth rates, non-zero for basal species
    # Per Romanuck, ri=1 for basal species
    rs = {}
    for node in G.nodes():
        if G.in_degree(node) == 0: #basal species have no indegree
            rs[node] = 1
        else:
            rs[node] = 0
    # Ks - carrying capacities 
    K = 1
    # metabolic rates
    # Per Romanuk, xi = 0.5 for all species
    x = 0.5
    # yij is the maximum rate at which species i assimilates species j per unit metabolic rate of species i
    # Per Romanuk yij = 6 for all species combinations
    y_arg = 6 #name to distinguish from state vector y
    # eij is the conversion efficiency with which the biomass of species j lost due to consumption by species i is converted into the biomass of species i
    # Per Romanuk eij = 1 for all species combinations
    e = 1
    # q is the parameter to determine the functional response
    q = 0.2
    # b0ij is the half saturation density
    # per Romanuk b0ij = 0.5 for all i, j
    b0 = 0.5
    
    #create a map from node i's in G to an index in the state vector
    index_dict = {}
    ct = 0 
    for node in G.nodes():
        index_dict[node] = ct
        ct +=1
    
    return (G,rs,K,x,y_arg,e,q,b0,index_dict)

def plot_dynamic_solution(sol):
    x = sol.t
    for i in range(0,len(sol.y)):
        plt.plot(x, sol.y[i], label = "species " + str(i))
    #plt.legend()
    plt.show()
    
def plot_dynamic_solution_nodes(sol,ns):
    x = sol.t
    for i in range(0,len(sol.y)):
        if i in ns: 
            plt.plot(x, sol.y[i], label = "species " + str(i))
    #plt.legend()
    plt.show()
    
def plot_dynamic_solution_inv_es(sol,inv,intermediate,es):
    x = sol.t
    for i in range(0,len(sol.y)):
        if i in inv:
            plt.plot(x, sol.y[i], label = "invader", c="red")
        if i in intermediate:
            plt.plot(x, sol.y[i], label = "species " + str(i), c="purple")
        if i in es:
            plt.plot(x, sol.y[i], label = "es species", c="yellow")
    # TEMPPP
    plt.ylim([0,0.25])
    plt.legend()
    plt.show()
    
def prune_web(G, b_curr,threshold):
    #remove all species that drop below a certain biomass (extinction threshold)
    
    new_web = nx.DiGraph()
    new_web.add_nodes_from(G.nodes())
    new_web.add_edges_from(G.edges())
    
    # index in b array to id in graph
    id_dict = {}
    ct = 0 
    for i in G.nodes():
        id_dict[ct] = i
        ct +=1
    
    for i in range(0,len(b_curr)):
        if b_curr[i] < threshold:
            print("pruned node " + str(id_dict[i]) + " with biomass " + str(b_curr[i]))
            new_web.remove_node(id_dict[i])
    
    return new_web

# function used to run a specified web for a specified number of timesteps
def run_dynamics(dargs):
    
    web_file = dargs[0]
    timesteps = dargs[1]
    plot = dargs[2]
    
    if len(dargs) > 3:
        id_p1 = dargs[3]
        id_p2 = dargs[4]
        id_p3 = dargs[5]
    
    start = time.time()
    print("running dynamics for -- " + web_file + " " + str(timesteps) + " timesteps...")

    web,nis,ris,cis,b_init = niche.read_web_from_file(web_file)

    #run dynamics for 2000 timesteps and prune web
    atn_args = get_params(web)

    # b has to be a vector, so map the node ids to 0 indexed values
    index_dict = atn_args[8]
    if plot:
        print("index dict: ")
        print(index_dict)
    b_init_vector = np.zeros((len(web.nodes())))
    for key in b_init:
        b_init_vector[index_dict[key]] = b_init[key]

    sol = integrate.solve_ivp(lambda t, y: atn(t, y, atn_args),(0,timesteps),b_init_vector,t_eval=range(0,timesteps))
    
    if plot:
        print("overall:")
        plot_dynamic_solution(sol)
        if len(dargs) > 3:
            nns1 = []
            for ii in id_p1:
                if ii in index_dict:
                    nns1.append(index_dict[ii])
            print("group 1:")
            #plot_dynamic_solution_nodes(sol,nns)
            nns2 = []
            for ii in id_p2:
                if ii in index_dict:
                    nns2.append(index_dict[ii])
            print("group 2:")
            #plot_dynamic_solution_nodes(sol,nns)
            nns3 = []
            for ii in id_p3:
                if ii in index_dict:
                    nns3.append(index_dict[ii])
            print("group 3:")
            #plot_dynamic_solution_nodes(sol,nns)
            plot_dynamic_solution_inv_es(sol,nns1,nns2,nns3)
    b_final = []
    for j in range(0,len(sol.y)):
        b_final.append(sol.y[j][len(sol.t)-1])
    web_final = prune_web(web,b_final,1e-10)

    final_S = len(web_final.nodes())
    
    print(web_file + " final nodes: " + str(final_S))
    print("final number of species after " + str(timesteps) + " timesteps: "  + str(final_S))    
    
    #write results to file
    nis_final = {}
    ris_final = {}
    cis_final = {}
    b_final_dict = {}
    for node in web_final.nodes(): 
        nis_final[node] = nis[node]
        ris_final[node] = ris[node]
        cis_final[node] = cis[node]
        b_final_dict[node] = b_final[index_dict[node]]

    niche.write_web_to_file(web_final,nis_final,ris_final,cis_final,b_final_dict, web_file + "_2000")

    end = time.time()
    print("time for " + web_file + " " + str(end - start))
    
