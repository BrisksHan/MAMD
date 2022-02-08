

from brainpy.composition import _parse_formula as parse_formula, calculate_mass
from brainpy.brainpy import _has_c

"""
if _has_c:
    from brainpy._c.composition import parse_formula as cparse_formula
"""

def test_parse(formula):
    composition = parse_formula(formula)
    print composition


def test_mass(formula):
    composition = parse_formula(formula)
    print composition.viewitems()
    print composition.mass()
    #self.assertAlmostEqual(composition.mass(), calculate_mass({"C": 6, "H": 12, "O": 6}))
#247.04269,C8H10FN3O3S
#247.10875,C10H18ClN3O2
#236.06847,C12H12O5
#360.06273,C13H16N2O8S
if __name__ == '__main__':
   #("C10H18ClN3O2")
   test_mass("C13H16N2O8S")
