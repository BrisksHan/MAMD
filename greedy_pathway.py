from graph import construct_graph,graph_traverse_greedy
from read_portals import read_LCMS_csv,read_compound_table
import numpy as np
import copy

def run_traverse(MF_list, pos = True, RTcutoff = 4.9999, correlation_cutoff = 0.7, LCMS_file_path = "/home/han/Downloads/data/test1/h5.csv" ,ppm = 3,include_pathway_check = True
                 ,minimum_score = 4):#3-5
    peakdata, peakintensities = read_LCMS_csv.read_csv_file(LCMS_file_path)
    graphs, raw_graph = construct_graph.get_fully_connected_graph(peakdata, peakintensities, RTcutoff, correlation_cutoff)
    print('number of isolated graphs:', len(graphs))
    
    #print(raw_graph.has_edge(1086, 346))
    #print(raw_graph.has_edge(1086, 1))
    #raise Exception('stop')

    primary, comprehensive = read_compound_table.main_user_defined_MF(pos, MF_list)#primary is not used as the information is all included in comprehensive 
    annotation_result = []

    for i in range(len(graphs)):
        #show_graph(graphs[i])
        single_graph_annotation_result = graph_traverse_greedy.attribute_graph_traverse(graphs[i], comprehensive, ppm = ppm)
        #print(single_graph_annotation_result)
        annotation_result += single_graph_annotation_result
        print "finished:",i," of ",len(graphs)
    print('number of community without check:',len(annotation_result))
    #print(annotation_result[0:10])

    print('combine annotation result')
    combined_result = combine_annotation_result(annotation_result)#only interested in peak modules

    #minimum_score_checked_result = get_pathway_annotation(combined_result, include_pathway_check, minimum_score)
    MFs, pathway_score = get_pathway_annotation_boost_score(combined_result, include_pathway_check)
    #print('annotation result:', annotation_result[0])
    #print('combined result:', combined_result[0])
    #print(minimum_score_checked_result[0])
    boosted_result = update_annotation_all_scores(annotation_result, MFs, pathway_score)#update the annotation result
    minimum_score_checked_result = remove_low_confidence_score(boosted_result, minimum_score)
    return minimum_score_checked_result, annotation_result
        

def combine_annotation_result(annotation_result):#get all unique MFs
    MFs = []
    Scores = []
    final_result = []
    for i in range(len(annotation_result)):
        current_MF = annotation_result[i][1]
        current_score = annotation_result[i][0]
        if (current_MF in MFs) == False:
            MFs.append(current_MF)
            Scores.append(current_score)
            final_result.append(annotation_result[i])

    return final_result

def update_annotation_all_scores(annotation_result, MFs, pathway_scores):#the new version
    final_annotation_result = copy.deepcopy(annotation_result)
    for item in final_annotation_result:
        current_MF = item[1]
        index = MFs.index(current_MF)
        item[0] += pathway_scores[index]

    return final_annotation_result

def remove_low_confidence_score(annotations, minimum_score):
    minimum_score_checked_result = []
    for i in range(len(annotations)):
        if annotations[i][0] >= minimum_score:
            minimum_score_checked_result.append(annotations[i])
    return minimum_score_checked_result

def get_kegg_names(annotation_result):
    ct = read_compound_table._read_buildin_compound_table()
    MFs_list = []
    for i in range(len(ct)):
        MFs_list.append(ct[i][0])
    corresponding_KEGG_IDs = []
    for i in range(len(annotation_result)):
        current_KEGG_IDs = []
        if annotation_result[i][1] in MFs_list:
            all_match = [t for t,k in enumerate(MFs_list) if k == annotation_result[i][1]]
        else:
            all_match = []
        #annotation_result[i][1] is the molecular formula
        for j in range(len(all_match)):
            if ct[all_match[j]][4] != "":
                current_KEGG_IDs.append(ct[all_match[j]][4])
        corresponding_KEGG_IDs.append(current_KEGG_IDs)
    return corresponding_KEGG_IDs

def get_pathway_annotation_boost_score(combined_result, include_pathway_check, up_limit = 1.5):
    kegg_IDs = get_kegg_names(combined_result)

    refined_reaction_table = read_compound_table.reaction_table_type_check()
    boost_pair_indexes = []

    MFs = []
    for item in combined_result:
        MFs.append(item[1])
    score = np.zeros(len(combined_result))

    if include_pathway_check == True:
        for i in range(len(kegg_IDs)):
            for j in range(i+1,len(kegg_IDs)):
                find_match = False
                for l in range(len(kegg_IDs[i])):
                    if find_match == True:
                        break
                    for k in range(len(kegg_IDs[j])):
                        if find_match == True:
                            break
                        reaction_combined = [kegg_IDs[i][l],kegg_IDs[j][k]]
                        for m in range(len(refined_reaction_table)):
                            if find_match == True:
                                break
                            intersected = list(set(refined_reaction_table[m]).intersection(reaction_combined))
                            if len(intersected) == 2:
                                current_boost_pair_indexes = [i, j]
                                boost_pair_indexes.append(current_boost_pair_indexes)
                                #print "find match:", current_boost_pair_indexes
                                find_match = True
        score = np.zeros(len(combined_result))
        
        for i in range(len(boost_pair_indexes)):
            indexA = boost_pair_indexes[i][0]
            indexB = boost_pair_indexes[i][1]
            if score[indexA] < up_limit:
                score[indexA] += 1
            if score[indexB] < up_limit:
                score[indexB] += 1

    return MFs, score