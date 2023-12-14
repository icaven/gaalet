"""
Simple example to show the interoperability between the python clifford package and the gaalet library using
Conformal Geometric Algebra
"""

import numpy as  np

# Suppress deprecation warnings from numba, imported indirectly using clifford.
import warnings
import numba.core.errors

# The warnings must be suppressed before importing clifford
warnings.filterwarnings('ignore', category=numba.core.errors.NumbaDeprecationWarning, module='numba')
warnings.filterwarnings('ignore', category=numba.core.errors.NumbaDeprecationWarning, module='clifford.tools.g3c')

# Suppress user warnings from pyganga, sometimes imported indirectly
# warnings.filterwarnings('ignore', category=UserWarning, module='pyganja')

from clifford.tools.g3c import random_euc_mv, normalise_n_minus_1, normalised
from clifford.g3c import layout

# The following module is the example module; you should change this to the name of the module in your application
# that uses the interface to gaalet
import clifford_gaalet_example


def main():
    # Create some random points
    points = []
    for pt_index in range(10):
        points.append(normalise_n_minus_1(layout.up(random_euc_mv())))

    lines_using_clifford = []
    print("Some random conformal points:")
    for p in points:
        print(p)
        # Compute the lines using clifford for later comparison
        lines_using_clifford.append(normalised(layout.eo ^ p ^ layout.einf))

    # Compute lines to the points from the origin
    lines = clifford_gaalet_example.lines_from_origin_to_points(points)
    print(f"Type of lines is {type(lines)}")
    print("Lines from the origin to each point")
    for index, l in enumerate(lines):
        print(l)
        if not np.all(np.isclose(l.value, lines_using_clifford[index].value)):
            print(f"Lines don't match {l} <=> {lines_using_clifford[index]}")

        # The outer product should always be zero because the point is on the line
        assert l ^ points[index] == 0, "Point not on line"


if __name__ == '__main__':
    main()
