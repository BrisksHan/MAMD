import csv
def write_result_to_csv(raw_results, output_path):
    print('start to write results')
    final_results = organise_results(raw_results)
    import csv
    with open(output_path, 'wb') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        spamwriter.writerows(final_results)

def organise_results(raw_results):
    final_list = []
    group_id = 1
    first_line = ['peak_id', 'measured_mz', 'measure_rt', 'MF', 'isotope', 'adduct', 'score', 'group_id']
    final_list.append(first_line)
    for item in raw_results:
        score = item[0]
        MF = item[1]
        for an_annotated_peak in item[2]:
            isotope = an_annotated_peak[5]
            adduct = an_annotated_peak[6]
            measured_mz = an_annotated_peak[9]
            measure_rt = an_annotated_peak[10]
            peak_id = an_annotated_peak[12]
            current_list = [peak_id, measured_mz, measure_rt, MF, isotope, adduct, score, group_id]
            final_list.append(current_list)
        group_id += 1
    return final_list
            