import numpy as np
import scipy.stats

def _pairwise_relation(data1, data2, i, j):
    nas = np.logical_or(np.isnan(data1), np.isnan(data2))
    data1normalized = _normalization_sum(data1[~nas])
    data2normalized = _normalization_sum(data2[~nas])
    return _chi_square_distance(data1normalized, data2normalized, i, j)


def _chi_square_distance(row1, row2, a, b):
    nominator = 0
    denominator = 0
    dimension = len(row1)
    for i in range(len(row1)):
        nominator += np.square(row1[i] - row2[i])
        denominator += np.square(row1[i] + row2[i])
    return nominator / denominator, dimension, a, b


def pairwise_Pearson_correlation(data1, data2):#this one will be used
    nas = np.logical_or(np.isnan(data1), np.isnan(data2))
    # print nasq
    tcorr = scipy.stats.pearsonr(data1[~nas], data2[~nas])
    # print "tcorr ",tcorr
    return tcorr


def _normalization_sum(row):
    sum = 0
    result = np.zeros(len(row))
    for i in range(len(row)):
        sum += row[i]
    for i in range(len(row)):
        result[i] = row[i] / sum * 100
    return result