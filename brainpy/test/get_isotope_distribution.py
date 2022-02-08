import unittest

from brainpy import IsotopicDistribution, isotopic_variants, calculate_mass, neutral_mass, _has_c, Peak
from brainpy.composition import _parse_formula as parse_formula

if _has_c:
    from brainpy.brainpy import _IsotopicDistribution


def isotope_distribution(formula):
    #hexnac = {'H': 13, 'C': 8, 'O': 4, 'N': 5}#243.09675,C8H13N5O4
    #243.10077,C13H13N3O2
    #247.04269,C8H10FN3O3S
    composition = parse_formula(formula)
    hexnac = composition.items()
    hexnac = dict((x,y) for x,y in hexnac)
    dist = isotopic_variants(hexnac)
    #print dist[0].mz#intensity charge mz
    #print dist[0].intensity
    return dist


if __name__ == '__main__':
    isotope_distribution("C8H13N5O4")