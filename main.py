from read_portals import read_LCMS_csv,calculate_theoretical_mz,read_compound_table, read_write
from argparse import ArgumentDefaultsHelpFormatter, ArgumentParser
import greedy_pathway
import csv
import write_result

def parse_args():
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, conflict_handler='resolve')
    parser.add_argument('--LCMS_path', default='demo/std1_pos_0.csv',
                        help='peak feature file path')
    parser.add_argument('--MODE', default='POS', choices = ['POS','NEG'],
                        help='LC-MS mode choose POS or NEG')
    parser.add_argument('--score', default=4,
                        help='the minimum confidence score')
    parser.add_argument('--ppm', default=5,
                        help='minimum tolorance')
    parser.add_argument('--CORR', default=0.7,
                        help='correlation cutoff between peak features')
    parser.add_argument('--RT', default=5,
                        help='retention time cutoff between peak features')
    parser.add_argument('--suspect_library', default='CUSTOM', choices=['KEGG','HMDB','CUSTOM'],
                        help='suspect metabolite library')
    parser.add_argument('--CUSTOM_path', default='demo/custom.csv',
                        help='path of user defined metabolite library')
    parser.add_argument('--delimeter', default='^',
                        help='delimeter for the custom csv file the default is ^')
    parser.add_argument('--output', default='output/annotation_output.csv',
                        help='output path')
    parser.add_argument('--pathway', default='True',
                        help='include kegg pathway analysis note that this will slow down the algorithm significantly')
    args = parser.parse_args()
    return args



def write_final_csv(lines, name):
    with open(name, "w") as f:
        writer = csv.writer(f)
        writer.writerows(lines)
        #return recall, precision, F1, n_recall, n_precision, n_F1
    
def run_annotation(args):
    if args.MODE == 'POS':
        mode = True
    elif args.MODE == 'NEG':
        mode = False
    else:
        raise Exception('LCMS_mode , please input POS or NEG')
    if args.suspect_library == 'KEGG':#not implemented yet
        library_infos = []
        input_MF_list = 'get_kegg'#not implemented yet
    elif args.suspect_library == 'HMDB':
        input_MF_list = 'get_HMDB'#not implemented yet
    elif args.suspect_library == 'CUSTOM':
        library_infos, input_MF_list = read_write.read_suspect_metabolite(args.CUSTOM_path)
    else:
        raise Exception('please input KEGG, HMDB or USER')
    output_path = args.output
    minimum_score = float(args.score)

    annotation_results, raw_results = graph_annotation_results = greedy_pathway.run_traverse(pos=mode, RTcutoff=float(args.RT), correlation_cutoff=float(args.CORR),
                                                LCMS_file_path=args.LCMS_path, MF_list = input_MF_list, ppm = float(args.ppm), include_pathway_check = args.pathway, minimum_score=minimum_score)
    #write_csv(args.output, list)
    print('module number:',len(annotation_results))
    write_result.write_result_to_csv(annotation_results, output_path)




if __name__ == "__main__":
    run_annotation(parse_args())