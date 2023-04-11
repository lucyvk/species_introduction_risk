import niche
import numpy as np

# get es metric assuming it's provided by one species instead 
# es_metric - change in es amount
# raw_final_es - the final es amount (not the change - used for looking at within-network variation)
def get_es_spec(b_init, b_final, es_spec):
    
    initial_es = es_linear(b_init[es_spec])
    
    if es_spec in b_final:
        final_es = es_linear(b_final[es_spec])
    else:
        final_es = 0    
    
    es_metric = (final_es - initial_es)
    
    return es_metric, final_es
    
#maps the value of species biomass to es value
def es_linear(b_val):
    return b_val*0.5

def get_rsd(vect):
    return abs(np.std(vect)/np.mean(vect))

# five metrics of interest returned:
# intro_metric1 - final biomass of introduced species 1
# intro_metric2 - final biomass of introduced species 2
# es_metric - % change in es amount
# bio_metric - fraction resident species remaining
# raw_final_es - the final es amount (not the change - used for looking at within-network variation)
def get_metrics(web_init, web_final, b_init,b_final,inv_id1,inv_id2,es_fn):
        
    #based off of aggregate of basal species biomasses
    init_int = niche.get_basal_ids(web_init)
    b_total_init = 0
    for node in init_int:
        b_total_init += b_init[node]
    
    initial_es = es_fn(b_total_init)
    initial_bio = len(b_init.keys())
    
    inv_adj = 0
    intro_metric1 = 0
    intro_metric2 = 0
    if inv_id1 and inv_id1 in b_final:
        intro_metric1 = b_final[inv_id1]
        inv_adj = inv_adj - 1
    if inv_id2 and inv_id2 in b_final:
        intro_metric2 = b_final[inv_id2]
        inv_adj = inv_adj - 1
    
    # total biomass of all the resident basal species
    final_int = niche.get_basal_ids(web_final)
    b_total_final = 0
    for node in final_int:
        if node != inv_id1 and node != inv_id2:
            b_total_final += b_final[node]
    
    final_es = es_fn(b_total_final)
    final_bio = len(b_final.keys()) + inv_adj # subtract out the invader counts
    
    es_metric = (final_es - initial_es)
    bio_metric = final_bio/initial_bio
    
    return intro_metric1, intro_metric2, es_metric, bio_metric, final_es

# get metrics assuming es provided by a single basal species - calculate for all
def get_metrics_single_es(b_init_web, b_final, webi):
    
    es_nodes = {}
    es_nodes_raw = {}
    for_rsd = []
    
    for node in niche.get_basal_ids(webi):
        ges = get_es_spec(b_init_web,b_final,node)
        for_rsd.append(ges[1])

    for node in niche.get_basal_ids(webi):
        if node not in es_nodes:
            es_nodes[node] = [get_es_spec(b_init_web,b_final,node)[0]]
        else:
            es_nodes[node].append(get_es_spec(b_init_web,b_final,node)[0])
        if node not in es_nodes_raw:
            es_nodes_raw[node] = [get_es_spec(b_init_web,b_final,node)[1]]
        else:
            es_nodes_raw[node].append(get_es_spec(b_init_web,b_final,node)[1])

    return for_rsd, es_nodes, es_nodes_raw

