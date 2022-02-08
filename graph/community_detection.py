from operator import itemgetter
import bisect

def graph_search(graph, comprehensive, ppm = 5, RT_tol = 5, postive_mode = True):
    node_info = graph.nodes()
    comprehensive_data_list_first_column = map(itemgetter(0), comprehensive)
    matched_index = []
    for i in range(len(node_info)):
        mz = graph.node[node_info[i]]['mz_value']
        tol = _calculate_mz_tolerance(mz, ppm = ppm)
        index_small, index_large = find_index_bisect(comprehensive_data_list_first_column,mz,tol)
        current_match = []
        #at_least_one_match = False
        if index_small != -1 and index_large != -1:
            for j in range(index_small,index_large):
                theoretical_mz = comprehensive_data_list_first_column[j]
                tol = _calculate_mz_tolerance(theoretical_mz)
                theoretical_mz_low = theoretical_mz - tol
                theoretical_mz_high = theoretical_mz + tol
                if  theoretical_mz_low <= mz <=  theoretical_mz_high:
                    current_match.append(j)
        matched_index.append(current_match)
    #for i in range(len(matched_index)):
    #    print "matched_index[i]:",matched_index[i]
    #step 1 a candidate set
    #step 2 a selection function to choose the best candidate
    #step 3 A feasibility function, that is used to determine if a candidate can be used to contribute to a solution
    #step 4 An objective function, which assigns a value to a solution, or a partial solution, and
    #step 5 A solution function, which will indicate when we have discovered a complete solution
    raw_result = []
    for i in range(len(matched_index)):
        for j in range(len(matched_index[i])):
            matched_comprehensive_index = matched_index[i][j]
            #print "comprehensive data:",comprehensive[matched_comprehensive_index]
            #[166.08625459669, False, True, 'C9H11NO2', 'C9H11NO2', 'M+H', '1', '1']
            if comprehensive[matched_comprehensive_index][1] == False:
                MF = comprehensive[matched_comprehensive_index][3]
                #community_info [selected_community, selected_community_comprehensive_index, final_score]
                community_info = _search_communities(i, matched_comprehensive_index, MF, graph, node_info, comprehensive,matched_index, RT_tol)#this part is for community detection part approach 1
                raw_result.append(community_info)

    selected_result = _selected_filete(raw_result, node_info, graph, comprehensive)
    return selected_result


def _search_communities(current_index, comprehensive_index, MF, graph, node_info, comprehensive, matched_index, RT_tol):#test function1
    selected_community = []#this will be the index of the metabolite features
    selected_community.append(current_index)
    selected_community_comprehensive_index = []
    selected_community_comprehensive_index.append(comprehensive_index)
    #score = adduct_score + isotope score - mu
    #score = summation(plus adducts + isotopes) - penality

    while True:
        feasible_neighbors_info, feasible_neighbors_score = _feasible_neighbors(MF, selected_community,
                                                                                selected_community_comprehensive_index,
                                                                                graph, node_info, comprehensive,
                                                                                matched_index, RT_tol)
        if len(feasible_neighbors_info) == 0:
            break
        else:
            #print "current neighbor info",feasible_neighbors_info
            #print "current neighbor score",feasible_neighbors_score
            max_score = max(feasible_neighbors_score)
            max_score_index = feasible_neighbors_score.index(max_score)
            selected_community.append(feasible_neighbors_info[max_score_index][0])
            selected_community_comprehensive_index.append(feasible_neighbors_info[max_score_index][1])
    #print "searched info p1:",selected_community
    #print "searched info p2:",selected_community_comprehensive_index
    final_score = _calculate_score(MF, selected_community, selected_community_comprehensive_index, comprehensive)
    return [selected_community, selected_community_comprehensive_index, final_score, MF]

def _calculate_score(MF, selected_community, selected_community_comprehensive_index, comprehensive, mu = 0.1):#test function1
    """
    calculate the score for each selected community
    :param MF:
    :param selected_community:
    :param selected_community_comprehensive_index:
    :param graph:
    :param node_info:
    :param comprehensive:
    :param mu:
    :return:
    """
    #graph.neighbors(start_node_name)
    #adduct _ isotope - mu (this could also be the average weight and some transofrm)
    #primary_ions = ['M+H', 'M+K', 'M+Na', 'M+NH4', 'M-H', 'M+Cl', 'M+Na-2H', 'M+K-2H']
    primary_ions = ['M+H', 'M+Na', 'M-H', 'M+Cl']
    adducts = []
    score = 0
    for i in range(len(selected_community)):
        #[166.08625459669, False, True, 'C9H11NO2', 'C9H11NO2', 'M+H', '1', '1']
        form = comprehensive[selected_community_comprehensive_index[i]][5]
        isotope_flag = comprehensive[selected_community_comprehensive_index[i]][1]
        if (form in primary_ions) == True:
            if isotope_flag == True:
                score +=2#isotope primary
            else:
                score +=3#non isotope primary
        else:
            if isotope_flag == True:
                score += 2
            else:
                score += 1
    #print "score:",score," ",selected_community
    total_count = len(selected_community)
    score -= total_count*mu

    return score


def _feasible_neighbors(MF, selected_community, selected_community_comprehensive_index, graph, node_info, comprehensive, matched_index, RT_tol):#test function1
    """
    select all the neighbors with a positive score
    if there are no positive score return empty list
    this part does take RT into consideration
    the

    :param MF:
    :param selected_community:
    :param selected_community_comprehensive_index:
    :param graph:
    :param node_info:
    :param comprehensive:
    :param matched_index:
    :param RT_tol:
    :return:
    """
    non_isotope_form = []
    RT_min,RT_max = _get_RT_range(selected_community, graph, node_info)
    for i in range(len(selected_community_comprehensive_index)):
        if comprehensive[selected_community_comprehensive_index[i]][1] == False:
            #adduct_form = comprehensive[selected_community_comprehensive_index]
            adduct_form = comprehensive[selected_community_comprehensive_index[i]][5]
            if (adduct_form in non_isotope_form) == False:
                non_isotope_form.append(adduct_form)
    feasible_neighbors_info = []
    feasible_neighbors_score = []
    un_selected_neighbors = []
    #primary_ions = ['M+H', 'M+K', 'M+Na', 'M+NH4', 'M-H', 'M+Cl', 'M+Na-2H', 'M+K-2H']
    primary_ions = ['M+H', 'M+Na', 'M-H', 'M+Cl']
    for i in range(len(selected_community)):
        current_node_index = selected_community[i]
        current_node = node_info[current_node_index]
        neighbors = graph.neighbors(current_node)
        for j in range(len(neighbors)):
            neighbor_index = node_info.index(neighbors[j])
            neighbor_RT = graph.node[neighbors[j]]['RT']
            if ((neighbor_index in selected_community) == False) and (_check_RT_range(RT_min,RT_max,neighbor_RT,RT_tol) == True):
            #this part make sure the when the new one is searched, there will not break the RT tolerence
                if (neighbor_index in un_selected_neighbors) == False:
                    un_selected_neighbors.append(neighbor_index)

    # [166.08625459669, False, True, 'C9H11NO2', 'C9H11NO2', 'M+H', '1', '1']
    # [1110.7824024027198, True, False, 'C58H109NO18', 'C58H109NO18-isotope', '2M+2H', '2', '2']
    for i in range(len(un_selected_neighbors)):
        current_node_index = un_selected_neighbors[i]
        current_node = node_info[current_node_index]
        current_matched_index = matched_index[current_node_index]
        for j in range(len(current_matched_index)):#this part went wrong ? have to figure it out
            comprehensive_index = current_matched_index[j]
            current_MF = comprehensive[comprehensive_index][3]
            adduct_form = comprehensive[comprehensive_index][5]
            if current_MF == MF:
                if comprehensive[comprehensive_index][1] == True :#if this is isotope metabolite feature
                    if (adduct_form in non_isotope_form) == True:
                        matched = [current_node_index, comprehensive_index]
                        feasible_neighbors_info.append(matched)
                        feasible_neighbors_score.append(2)#isotope pattern
                else:
                    #if this is not a isotope metabolite feature
                    ion_form = comprehensive[comprehensive_index][5]
                    if (ion_form in primary_ions) == True:
                        matched = [current_node_index, comprehensive_index]
                        feasible_neighbors_info.append(matched)
                        feasible_neighbors_score.append(3)#non isotope primary
                    else:
                        matched = [current_node_index, comprehensive_index]
                        feasible_neighbors_info.append(matched)
                        feasible_neighbors_score.append(1)#non isotope non primary

    return feasible_neighbors_info, feasible_neighbors_score


def _get_RT_range(selected_community, graph, node_info):
    """
    get the RT range from the selected community
    this is for choosing
    :param selected_community:
    :param graph:
    :param node_info:
    :return:
    """
    #mz = graph.node[node_info[i]]['mz_value']
    first_index = selected_community[0]
    RT = graph.node[node_info[first_index]]['RT']
    RT_min = RT
    RT_max = RT
    for i in range(1,len(selected_community)):
        current_node_index = selected_community[i]
        current_RT = graph.node[node_info[current_node_index]]['RT']
        if current_RT < RT_min:
            RT_min = current_RT
        if current_RT > RT_max:
            RT_max = current_RT
    return RT_min,RT_max

def _check_RT_range(RT_min,RT_max,new_RT,RT_threshold):
    """
    check whether the new one is within the RT tolorence
    :param RT_min:
    :param RT_max:
    :param new_RT:
    :param RT_threshold:
    :return:
    """
    result = False
    if new_RT >= RT_min and new_RT <= RT_max:
        result = True
    elif new_RT > RT_max:
        if (new_RT - RT_min) < RT_threshold:
            result = True
    elif new_RT < RT_min:
        if (RT_max - new_RT) < RT_threshold:
            result = True
    return result

def _calculate_mz_tolerance(mass, ppm = 5.0):
    tol = mass * 0.000001 * ppm
    return tol

def find_index_bisect(sorted_list,target_mz,offset = 0.1):#this need a bit work !!!!!!!!
    """
    set the price without
    :param sorted_list:
    :param target_mz:
    :param offset:
    :return:
    """
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

def _selected_filete(raw_result, node_info, graph, comprehensive):
    #[selected_community, selected_community_comprehensive_index, final_score, MF]
    #[1110.7824024027198, True, False, 'C58H109NO18', 'C58H109NO18-isotope', '2M+2H', '2', '2']

    print "original raw result:",raw_result

    MF_list = []
    MF_column = []

    selected_result = []

    for i in range(len(raw_result)):
        MF_column.append(raw_result[i][3])

    for i in range(len(raw_result)):
        if (raw_result[i][3] in MF_list) == False:
            MF_list.append(raw_result[i][3])

    for list_index in range(len(MF_list)):
        indices = [i for i, x in enumerate(MF_column) if x == MF_list[list_index]]
        max_score = raw_result[indices[0]][2]
        choosen_indices = indices[0]
        for j in range(1,len(indices)):
            new_score = raw_result[indices[j]][2]
            if new_score > max_score:
                choosen_indices = indices[j]
        single_annotation = _get_organized_annotation_result(raw_result[choosen_indices], node_info, graph, comprehensive)
        selected_result.append(single_annotation)

    return selected_result

def _get_organized_annotation_result(raw_result, node_info, graph, comprehensive):
    organized_result = []
    organized_result.append(raw_result[3])
    organized_result.append(raw_result[2])
    communtiy_info = []
    #print "raw_result:",raw_result
    for i in range(len(raw_result[0])):
        #[selected_community, selected_community_comprehensive_index, final_score, MF]
        node_index = raw_result[0][i]
        node_name = node_info[node_index]
        measured_mz = graph.node[node_name]['mz_value']
        measured_RT = graph.node[node_name]['RT']
        comprehensive_index = raw_result[1][i]
        comprehensive_info = comprehensive[comprehensive_index]
        current_info = [node_name, measured_mz, measured_RT, comprehensive_info]
        communtiy_info.append(current_info)
    organized_result.append(communtiy_info)
    return organized_result









