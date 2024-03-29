import numpy as np
import random
import networkx as nx
import matplotlib.pyplot as plt
import scipy
import time
import pandas as pd
import csv
import os
from scipy.stats import beta

# Niche model generation
def niche_model(S,C):
        
    #generate parameters for a node
    def get_params():
        
        B = (1 - 2*C)/(2*C)
        
        #assign each species a niche value, ni
        ni = random.uniform(0,1)
      
        #assign each species a range, ri
        y = random.uniform(0,1)
        x = 1 - (1-y)**(1/B)
        #x = beta.rvs(1, B, size=1)[0]
        
        ri = x*ni
        
        #assign each species a center of the range, ci
        ci = random.uniform(ri/2,ni)
        
        return [ni,ri,ci]
            
    def assign_links(curr_web,nis,ris,cis):
        #figure out who each species eats
        for i in range(0,S):
            lower_bound = cis[i] - ris[i]/2
            upper_bound = cis[i] + ris[i]/2
            for j in range(0,S):
                if lower_bound < nis[j] and nis[j] < upper_bound:
                    curr_web.add_edge(j,i)
        return curr_web
            
    #helper function to find trophically identical species 
    def find_trophic_id(curr_web):    
                
        edge_check = {}
        trophic_id = []
        for node in curr_web.nodes(): 
                      
            #gather set of neighbors
            in_neighs = []
            in_neigh_list = list(curr_web.in_edges([node]))
            for in_neigh in in_neigh_list: 
                in_neighs.append(in_neigh[0])
            in_neighs = tuple(in_neighs)
            out_neighs = []
            out_neigh_list = list(curr_web.out_edges([node]))
            for out_neigh in out_neigh_list:
                out_neighs.append(out_neigh[1])
            out_neighs = tuple(out_neighs)

            #check to see if another species has the same neighbor set
            check_key = (in_neighs, out_neighs)
            if check_key not in edge_check:
                edge_check[check_key] = [node]
            else:
                edge_check[check_key].append(node)
                
        #save all but one of the identical species to be replaced
        for edge_k in edge_check:
            if (len(edge_check[edge_k]) > 1):
                r_ct = 0
                while (r_ct < len(edge_check[edge_k])):
                    to_replace = np.random.choice(edge_check[edge_k])
                    if (to_replace not in trophic_id):
                        trophic_id.append(to_replace)
                        r_ct += 1
            
        return trophic_id
   
    def assign_basal(nis,ris):
        #set ri to 0 for species with smallest ni (basal)
        min_ni = 1000000000
        min_ni_id = -1
        for key in nis.keys():
            if nis[key] < min_ni:
                min_ni = nis[key]
                min_ni_id = key

        ris[min_ni_id] = 0
        return ris
   
    def get_problem_species(web):
        isos = []
        for iso in nx.isolates(web):
            isos.append(iso)
        ps = isos + find_trophic_id(web)
        return ps
        
    nis = {}
    ris = {}
    cis = {}
        
    web = nx.DiGraph()
    for i in range(0,S):
        params = get_params()
        nis[i] = params[0]
        ris[i] = params[1]
        cis[i] = params[2]        
        #add a node for each species
        web.add_node(i)
    ris = assign_basal(nis,ris)
    web = assign_links(web,nis,ris,cis)

    #Check for and replace isolates and trophically identical species
    problem_species = get_problem_species(web)
    while (len(problem_species) > 0):
        
        # replace one at a time
        prob = problem_species[0]
        params = get_params()
        nis[prob] = params[0]
        ris[prob] = params[1]
        cis[prob] = params[2]
        ris = assign_basal(nis,ris)
        web2 = nx.DiGraph()
        web2.add_nodes_from(web)
        web2 = assign_links(web2,nis,ris,cis)
        web = web2
        problem_species = get_problem_species(web)
        
    return nis, ris, cis, web

# Returns whether the potential synthetic food web is "plausible" in terms of energy flow
# You don't want a loop with has no basal food source
def check_web(G):

    # detect loops
    loops = nx.algorithms.cycles.simple_cycles(G)
    
    # check if each loop has an "external energy source" 
    # this means that there is some path from one of the nodes in the loop back to a basal species
    # so, travelling back via in neighbors, one eventually reaches one that has no in neighbors itself
    lct = 0 
    for loop in loops:
        loop_okay = False
        for node in loop:
            if not loop_okay:
                #run DFS for every node in the loop (via in neighbors)
                visited = []
                stack = []
                stack.append(node)
                #keep_checking_node = True
                while len(stack) != 0: #and keep_checking_node:
                    top = stack.pop()
                    #print("top: " + str(top))
                    if top not in visited:
                        visited.append(top)
                    in_neigh_list = list(G.in_edges([top]))
                    for inl in range(0,len(in_neigh_list)):
                        in_neigh_list[inl] = in_neigh_list[inl][0]
                    for in_neigh in in_neigh_list:
                        if len(list(G.in_edges([in_neigh]))) == 0:
                            #reached a basal species, this loop is okay
                            loop_okay = True
                            #can immediately end this while loop
                            break
                        if in_neigh not in visited:
                            stack.append(in_neigh)
        lct += 1
        if not loop_okay:
            #if any loop isn't okay, the whole web isn't okay, return False
            return False
    
    return True # if all loops are okay, return that the web is plausible

# helper function to write the food web to a file in a standard format
# Three files: node list, edge list, niche model parameter + biomass details 
def write_web_to_file(weby,nis,ris,cis,b,stem):
    
    #node file
    node_file = open(stem + "_nodes.txt","w")
    
    for node in weby.nodes():
        
        node_file.write(str(node) + '\n')
    
    #directed edge file (u->v)
    edge_file = open(stem + "_edges.txt","w")
    for edge in weby.edges():
        edge_file.write(str(edge[0]) + ' ' + str(edge[1]) + '\n')
        
    columns = ["n", "r", "c", "b"]
    rows = weby.nodes()
    data_rows = []
    for i in rows:
        data_rows.append([nis[i],ris[i],cis[i],b[i]])
    data = np.array(data_rows)
    df = pd.DataFrame(data=data, index=rows, columns=columns)
    df.to_csv(stem + "_details.csv")
        
    node_file.close()
    edge_file.close()

# helper function to generate links based on niche model parameters (same as in niche function) 
def assign_links(webb,cc,rr,nn):
    # figure out who each species eats
    S = webb.nodes()
    for x in S:
        for y in S:
            # x eats those with ny in rx around cx
            lower_bound = cc[x] - rr[x]/2
            upper_bound = cc[x] + rr[x]/2
            if lower_bound < nn[y] and nn[y] < upper_bound:
                # x eats y, add edge from y to x 
                webb.add_edge(y,x)
    return webb

# helper function to read the food web from file using the standard format
def read_web_from_file(stem):
    
    web = nx.DiGraph()
    
    node_file = open(stem + "_nodes.txt","r")
    for node in node_file:
        web.add_node(int(node))
    
    edge_file = open(stem + "_edges.txt","r")
    for edge in edge_file:
        web.add_edge(int(edge.split()[0]),int(edge.split()[1]))
        
    nis = {}
    ris = {}
    cis = {}
    b_curr = {}
    deets = open(stem + "_details.csv")
    deets_reader = csv.reader(deets)
    next(deets_reader)
    for row in deets_reader:
        nis[int(row[0])] = float(row[1])
        ris[int(row[0])] = float(row[2])
        cis[int(row[0])] = float(row[3])
        b_curr[int(row[0])] = float(row[4])
    
    node_file.close()
    edge_file.close()
    deets.close()
    
    return web,nis,ris,cis,b_curr

def viz_web(G):
    nx.draw_random(G,with_labels=True)
    plt.show()
    
def viz_web_inv_es_int(G,inv_ids,es_ids,int_ids):
    colors = []
    posi = {}
    for i in G.nodes():
        if i in int_ids:
            ypos = 0.2
            xpos = np.random.random()*0.8 + 0.1
        elif i in es_ids:
            ypos = 0.1
            xpos = 0.5
        elif G.out_degree(i) > 0 and G.in_degree(i) == 0 :# has predators but no prey 
            #ypos = np.random.random()*0.2 + 0.1
            ypos = 0.1
            xpos = np.random.random()*0.8 + 0.1
        elif G.in_degree(i) > 0 and G.out_degree(i) == 0 : # has prey but no predators - none of these?
            #ypos = np.random.random()*0.2 + 0.7
            ypos = 0.7
            xpos = np.random.random()*0.8 + 0.1
        elif G.out_degree(i) > 0 and G.in_degree(i) > 0: # in the middle
            #ypos = np.random.random()*0.2 + 0.4
            ypos = 0.4
            xpos = np.random.random()*0.8 + 0.1
        else: #disconnected - shouldn't happen?
            ypos = np.random.random()*0.2 + 0.1
            xpos = np.random.random()*0.1
        posi[i] = [xpos,ypos]
        
        if i in inv_ids:
            colors.append('r')
        elif i in es_ids:
            colors.append('y')
        elif i in int_ids:
            colors.append('purple')
        else:
            colors.append('b')
    nx.draw_networkx(G,with_labels=True,pos=posi,node_color=colors)
    plt.show()
    
# get basal species (only consumers, no resources)
def get_basal_ids(net):
    basal_species = []
    
    for node in net.nodes():
        
        out_neighs = []
        out_neigh_list = list(net.out_edges([node]))
        for out_neigh in out_neigh_list:
            out_neighs.append(out_neigh[1])
            
        in_neighs = []
        in_neigh_list = list(net.in_edges([node]))
        for in_neigh in in_neigh_list: 
            in_neighs.append(in_neigh[0])
            
        # consumers but no resources 
        if len(out_neighs) > 0 and len(in_neighs) == 0:
            basal_species.append(node)
            
    return basal_species

# get intermediate species (both resources and consumers)
def get_intermediate_ids(net):
    int_species = []
    for node in net.nodes():
        
        #top predator would have no outgoing edges
        out_neighs = []
        out_neigh_list = list(net.out_edges([node]))
        for out_neigh in out_neigh_list:
            out_neighs.append(out_neigh[1])
            
        in_neighs = []
        in_neigh_list = list(net.in_edges([node]))
        for in_neigh in in_neigh_list: 
            in_neighs.append(in_neigh[0])
            
        if len(out_neighs) > 0 and len(in_neighs) > 0:
            int_species.append(node)
            
    return int_species
    