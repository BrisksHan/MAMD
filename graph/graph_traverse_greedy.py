import networkx as nx
import numpy as np
import bisect
from operator import itemgetter





def attribute_graph_traverse(graph, comprehensive, ppm = 3, RT_tol = 5, postive_mode = True):#start point
    node_info = graph.nodes()
    comprehensive_data_list_first_column = map(itemgetter(0), comprehensive)
    matched_index = []
    matched_index_nodes = []
    #print "start calculating start points"
    for i in range(len(node_info)):#store all matched nodes
        #print "peak ID:",i
        mz = graph.node[node_info[i]]['mz_value']
        tol = _calculate_mz_tolerance(mz, ppm = ppm)
        index_small, index_large = find_index_bisect(comprehensive_data_list_first_column,mz,tol)#find matched m/z
        current_match = []
        at_least_one_match = False
        if index_small != -1 and index_large != -1:
            for j in range(index_small,index_large):
                theoretical_mz = comprehensive_data_list_first_column[j]
                tol = _calculate_mz_tolerance(theoretical_mz)
                theoretical_mz_low = theoretical_mz - tol
                theoretical_mz_high = theoretical_mz + tol
                if  theoretical_mz_low <= mz <=  theoretical_mz_high:
                    #print core_data_list[j][3]
                    at_least_one_match = True
                    current_match.append(j)

        if at_least_one_match == True:
            matched_index.append(current_match)
            matched_index_nodes.append(node_info[i])#this store the node name

    start_points = _get_start_node(matched_index, comprehensive)#get the node names
    #print "start points calculation finished"
    traverse_result = []
    explored_combinations = []
    from tqdm import tqdm

    for i in tqdm(range(len(start_points))):
        #traverse_result.append(_BFS(start_points[i],graph,node_info,matched_index, matched_index_nodes,comprehensive,RT_tol))

        start_node_index = start_points[i][0]
        start_node_name = matched_index_nodes[start_node_index]
        start_node_comprehensive_index = start_points[i][5]
        

        #annotated_module, current_explored_combinations = _greedy_search_annotation(start_points[i],graph,node_info,matched_index, matched_index_nodes,comprehensive,RT_tol)
        #traverse_result.append(annotated_module)
        current_start_points_organised = str(start_node_name) + '_' + str(start_node_comprehensive_index)
        if (current_start_points_organised in explored_combinations) == False:
            annotated_module, current_explored_combinations = _greedy_search_annotation(start_points[i],graph,node_info,matched_index, matched_index_nodes,comprehensive,RT_tol)
            traverse_result.append(annotated_module)
            explored_combinations += current_explored_combinations
            explored_combinations = list(set(explored_combinations))

    return traverse_result

def _greedy_search_annotation(start_point, graph, node_info, matched_index, matched_index_nodes, comprehensive, RT_tol):
    """
    start node: 
    """
    primary_ions = ['M+H','M+K','M+Na','M+NH4','M-H','M+Cl','M+Na-2H','M+K-2H']

    selected_MF = start_point[2]
    start_node_index = start_point[0]
    
    start_node_name = matched_index_nodes[start_node_index]

    start_point_RT = graph.node[start_node_name]['RT']
    start_point_MPA = graph.node[start_node_name]['MPA']
    start_point_adduct = start_point[3]
    peak_abundant_indicator = start_point[4]
    start_ponnt_comprehensive_index = start_point[5]

    a_peak_module = peak_module(selected_MF, graph.copy(), RT_tol, comprehensive)
    
    if start_point_adduct in primary_ions:
        start_point_score = 3
    else:
        start_point_score = 1

    a_peak_module.add_peak(start_node_name, start_point_RT, start_point_MPA, peak_abundant_indicator, start_point_adduct, start_point_score, start_ponnt_comprehensive_index)

    module_node_names = [start_node_name]
    with_candidiate_neighbours = True

    while with_candidiate_neighbours == True:
        current_neighbours = []
        for a_node in module_node_names:
            first_order_neighbors = list(graph.neighbors(a_node))
            current_neighbours += first_order_neighbors
        current_neighbours = list(set(current_neighbours)- set(module_node_names))
        candidate_neighbours = []
        candidate_neighbours_scores = []
        if len(current_neighbours) > 0:
            for a_node in current_neighbours:
                #if 
                current_node_name = a_node
                current_node_RT = graph.node[current_node_name]['RT']
                current_node_MPA = graph.node[current_node_name]['MPA']
                try:
                    current_node_match_index = matched_index_nodes.index(current_node_name)
                    current_node_matches = matched_index[current_node_match_index]
                except:
                    continue
                if a_peak_module._check_RT_violation(current_node_RT) != True:#first test
                    continue
                for a_match in current_node_matches:
                    current_adduct = comprehensive[a_match][5]
                    current_MF = comprehensive[a_match][3]
                    current_abundant_indicator = comprehensive[a_match][1]
                    if current_MF != selected_MF:#second test
                        continue
                    if current_abundant_indicator == False:
                        candidate_neighbours.append([current_node_name, current_node_RT, current_node_MPA, current_adduct, current_abundant_indicator, a_match])
                        if current_adduct in primary_ions:
                            candidate_neighbours_scores.append(3)
                        else:
                            candidate_neighbours_scores.append(1)
                    else:
                        istope_checker = a_peak_module._check_isotope(current_node_name, current_node_MPA, current_adduct, current_node_RT)
                        if istope_checker == True:
                            candidate_neighbours.append([current_node_name, current_node_RT, current_node_MPA, current_adduct, current_abundant_indicator, a_match])
                            candidate_neighbours_scores.append(2)
        else:
            with_candidiate_neighbours = False
        if len(candidate_neighbours_scores) > 0:
            max_score = max(candidate_neighbours_scores)
            selected_neighbour_index = candidate_neighbours_scores.index(max_score)
            selected_candidate = candidate_neighbours[selected_neighbour_index]
            selected_node_name = selected_candidate[0]
            selected_node_RT = selected_candidate[1]
            selected_node_MPA = selected_candidate[2]
            selected_node_adduct = selected_candidate[3]
            selected_node_indicator = selected_candidate[4]
            selected_node_comprehensive =selected_candidate[5]
            a_peak_module.add_peak(selected_node_name, selected_node_RT, selected_node_MPA, selected_node_indicator, selected_node_adduct, max_score, selected_node_comprehensive)
            module_node_names.append(selected_node_name)
        else:
            with_candidiate_neighbours = False
    return a_peak_module._organise_output(), a_peak_module._organise_annotations()

    

class peak_module:
    def __init__(self, MF, graph, RT_tol, comprehensive):
        self.RT_dict = {}
        self.peak_names = []
        self.peak_MPA = {}
        self.peak_abundant_indicator = {}
        self.MF = MF
        self.graph = graph
        self.peak_adducts = {}
        self.RT_tol = RT_tol
        self.score = 0
        self.comprehensive_matches = {}
        self.comprehensive = comprehensive
    
    def add_peak(self, peak_name, peak_RT, peak_MPA, peak_abundant_indicator, adduct, score, comprehensive_match):
        self.peak_names.append(peak_name)
        self.RT_dict[peak_name] = peak_RT
        self.peak_MPA[peak_name] = peak_MPA
        self.peak_abundant_indicator = peak_abundant_indicator
        self.peak_adducts[peak_name] = adduct
        self.score += score
        self.comprehensive_matches[peak_name] = comprehensive_match

    def evaludate_peak_impact(self, peak, peak_RT, peak_MPA, peak_abundant_indicator):
        return False

    def _evaluate_node(self, selected_MF, node_matched_index, comprehensive):
        find_match = False
        for i in range(len(node_matched_index)):
            ion_index = node_matched_index[i]
            ion = comprehensive[ion_index][3]
            if ion == selected_MF:
                find_match = True
                break
        return find_match

    def _check_RT_violation(self, new_RT):#check RT
        existing_RTs = list(self.RT_dict.values())
        result = True
        for item in existing_RTs:
            if abs(new_RT-item) > self.RT_tol:
                result = False
        return result

    def _check_isotope(self, new_node_name, new_node_MPA, new_adduct, new_node_RT):
        result = True
        existing_adducts = list(self.peak_adducts.values())
        if (new_adduct in existing_adducts) == False:
            return False
        else:
            existing_nodes = list(self.peak_adducts.keys())
            for node in existing_nodes:
                if self.peak_adducts[node] == new_adduct:
                    existing_MPA = self.peak_MPA[node]
                    if existing_MPA < new_node_MPA:
                        return False
                    else:
                        if self.graph.has_edge(new_node_name, node):
                            distance = abs(self.RT_dict[node] - new_node_RT)
                            if self.graph.get_edge_data(new_node_name, node, default=0)['weight'] > 0.7 and distance < 0.5:
                                return True
                            else:
                                return False
                        else:
                            return False

    def _organise_annotations(self):
        organised_names = []
        for item in self.peak_names:
            comprehensive_match = self.comprehensive_matches[item]
            name = str(item) + '_' + str(comprehensive_match)
            organised_names.append(name)
        return organised_names

    def _organise_output(self):
        #traverse_result = [selected_MF, MF_matched_count, MF_matched_indexes]
        #[max_score, annotations[max_count_index][0][4], annotations[max_count_index]]#First:score Second:MF Third:annotations
        #matched_info = [j] + comprehensive[matched_index[node_matched_index][l]] + [graph.node[node_name]["mz_value"]] + [graph.node[node_name]["RT"]] + [graph.node[node_name]["MPA"]] + [node_name]
        output_score = self.score
        output_MF = self.MF
        all_annotation_results = []
        all_annotation_results.append(output_score)#[output_MF, output_score]
        all_annotation_results.append(output_MF)
        comprehensive_infos = []
        for a_name in self.peak_names:
            mz = self.graph.node[a_name]["mz_value"]
            RT = self.graph.node[a_name]["RT"]
            mpa = self.graph.node[a_name]["MPA"]
            comprehensive_info = self.comprehensive[self.comprehensive_matches[a_name]]
            current_info = ['none'] + comprehensive_info + [mz] + [RT] + [mpa] + [a_name]
            comprehensive_infos.append(current_info)
        all_annotation_results.append(comprehensive_infos)
        return all_annotation_results

def _get_start_node(matched_index, comprehensive):
    start_point = []

    for i in range(len(matched_index)):
        for j in range(len(matched_index[i])):
            current = matched_index[i][j]
            #print(comprehensive[current])
            isotope_flag = comprehensive[current][1]
            #primary_flag = comprehensive[current][2]
            if isotope_flag == False:
                point_info = [i,j,comprehensive[current][3], comprehensive[current][5], comprehensive[current][1], current]#No.1 node ID No.2 chossed node ID, No.3 ion ID 4. adduct type
                start_point.append(point_info)

    return start_point



def _evaluate_node(selected_MF,node_matched_index,comprehensive):
    find_match = False
    for i in range(len(node_matched_index)):
        ion_index = node_matched_index[i]
        ion = comprehensive[ion_index][3]
        if ion == selected_MF:
            find_match = True
            break
    return find_match





def find_index_bisect(sorted_list,target_mz,offset = 0.1):#this need a bit work !!!!!!!!
    #print target_mz-offset
    #print offset
    bisect_index_small = bisect.bisect(sorted_list, (target_mz-offset))
    bisect_index_large = bisect.bisect(sorted_list, (target_mz+offset))
    list_length = len(sorted_list)
    smallest = 0
    largest = 0
    if bisect_index_small == list_length and sorted_list:
        largest_index = list_length - 1
        list_data = sorted_list[largest_index]
        if list_data + offset < target_mz:
            bisect_index_small = -1
            bisect_index_large = -1

    elif bisect_index_large == 0:
        list_data = sorted_list[0]
        if list_data - offset > target_mz:
            bisect_index_small = -1
            bisect_index_large = -1
    return bisect_index_small, bisect_index_large

def _calculate_mz_tolerance(mass, ppm = 5.0):
    tol = mass * 0.000001 * ppm
    return tol

def _evaluate_treavel_result(traverse_result,minimum_match = 2):
    #selected_MF, MF_matched_count, MF_matched_indexes
    #print "traverse_result:",traverse_result[0]
    MFs = [row[0] for row in traverse_result]
    MFs_set = set(MFs)
    MFs_no_duplicate = list(MFs_set)
    result = []
    for i in range(len(MFs_no_duplicate)):
        current_MF = MFs_no_duplicate[i]
        #print "current_MF:",current_MF
        count = []
        index = []
        for j in range(len(traverse_result)):
            if traverse_result[j][0] == current_MF:
                #print "get a match:",traverse_result[j]
                count.append(traverse_result[j][1])
                #replace this part and use score !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                index.append(j)

        max_count = max(count)
        #print "count:",count
        max_count_index = count.index(max_count)
        choosen_index = index[max_count_index]
        #print "max_count_index:",max_count_index
        if max_count >= minimum_match:
            result.append(traverse_result[choosen_index])
    return result

def _show_ion_forms(annotation_result,node_info,matched_index,matched_index_nodes,comprehensive, graph):
    #print(matched_index_nodes
    for i in range(len(annotation_result)):
        #print annotation_result[i]
        MF = annotation_result[i][0]
        nodes = annotation_result[i][2]
        for j in range(len(nodes)):
            #node_index = matched_index_nodes.index(nodes[j])
            node_name = node_info[nodes[j]]
            node_matched_index = matched_index_nodes.index(node_name)
            #find_match = False
            for l in range(len(matched_index[node_matched_index])):
                if comprehensive[matched_index[node_matched_index][l]][3] == MF:
                    #find_match = True
                    print(j," ", comprehensive[matched_index[node_matched_index][l]],"  ",graph.node[node_name]["mz_value"]," ",graph.node[node_name]["RT"]," ",graph.node[node_name]["MPA"])
                    #break


