"""
class graph(object):
    def __init__(self,Id,mz,rt,mpa):
        self.Id = Id
        self.mz = mz
        self.rt = rt
        self.mpa = mpa
        self.initial_info()
"""
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import correlation_calculation


def construct_LCMS_data_attribute_graph(LCMS_data,peakintensities,RT_distance_cutoff,correlation_cutoff):
    #input -1 this will be fully connected
    #ID mz rt MPA
    graph = nx.Graph()

    MF_ID = []
    for i in range(len(LCMS_data)):
        MF_ID.append(int(LCMS_data[i][0]))
    MF_ID_No_Duplicate = set(MF_ID)
    if len(MF_ID_No_Duplicate) != len(MF_ID):
        raise Exception("duplicate metabolite feature ID detected!")

    RT_list = LCMS_data[:,2]
    for i in range(len(LCMS_data)):
        graph.add_node(int(LCMS_data[i][0]), mz_value = LCMS_data[i][1], RT = LCMS_data[i][2], MPA = LCMS_data[i][3])


    for i in range(len(peakintensities)):
        for j in range(i+1,len(peakintensities)):
            RT1 = LCMS_data[i][2]
            RT2 = LCMS_data[j][2]
            RT_distance = abs(RT1-RT2)
            if RT_distance < RT_distance_cutoff or RT_distance_cutoff == -1:
                Pearson_correlation, p_value = correlation_calculation.pairwise_Pearson_correlation(peakintensities[i], peakintensities[j])
                if Pearson_correlation >= correlation_cutoff:# and p_value < 0.051: #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! I made change here
                    graph.add_edge(int(LCMS_data[i][0]),int(LCMS_data[j][0]),weight = Pearson_correlation)
    print "graph constructed"
    return graph

#def remove_isolated_nodes()

def find_connected_component_set(graph):
    raw_connected_components = list(nx.connected_components(graph))#G[u][v]['weight']
    refined_connected_components = []
    for i in range(len(raw_connected_components)):
        if len(raw_connected_components[i]) > 1:
            refined_connected_components.append(raw_connected_components[i])

    return refined_connected_components


def get_connected_subgraph(graph, refined_connected_component):
    subgraphs = []
    for i in range(len((refined_connected_component))):
        subgraph = graph.subgraph(refined_connected_component[i])
        subgraphs.append(subgraph)

    return subgraphs

def get_fully_connected_graph(LCMS_data,peakintensities,RT_distance_cutoff,correlation_cutoff):
    raw_graph = construct_LCMS_data_attribute_graph(LCMS_data,peakintensities,RT_distance_cutoff,correlation_cutoff)
    connected_sets = find_connected_component_set(raw_graph)
    subgraphs = get_connected_subgraph(raw_graph,connected_sets)
    return subgraphs, raw_graph

def get_raw_graph(LCMS_data,peakintensities,RT_distance_cutoff,correlation_cutoff):
    raw_graph = construct_LCMS_data_attribute_graph(LCMS_data,peakintensities,RT_distance_cutoff,correlation_cutoff)
    return raw_graph
#(3, weight=0.4, UTM=('13S', 382871, 3972649))

if __name__ == '__main__':
    construct_LCMS_data_attribute_graph(12)


