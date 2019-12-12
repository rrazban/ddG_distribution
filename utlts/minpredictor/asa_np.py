#!/usr/bin/env python

"""
Routines to calculate the Accessible Surface Area of a set of atoms.
The algorithm is adapted from the Rose lab's chasa.py, which uses
the dot density technique found in:

Shrake, A., and J. A. Rupley. "Environment and Exposure to Solvent
of Protein Atoms. Lysozyme and Insulin." JMB (1973) 79:351-371.
"""

import numpy as np
import math
from scipy.spatial.distance import cdist


def get_radii_from_prody_pdb(pdb):
    '''
    This is a helper function to avoid using the molecule
    class and handle AtomGroup classes from prody after
    a parsePDB loads a PDB from disk.
    '''

    import prody as pr

    atom_radii = { 
    'H': 1.20,
    'N': 1.55,
    'NA': 2.27,
    'CU': 1.40,
    'CL': 1.75,
    'C': 1.70,
    'O': 1.52,
    'I': 1.98,
    'P': 1.80,
    'B': 1.85,
    'BR': 1.85,
    'S': 1.80,
    'SE': 1.90,
    'F': 1.47,
    'FE': 1.80,
    'K':  2.75,
    'MN': 1.73,
    'MG': 1.73,
    'ZN': 1.39,
    'HG': 1.8,
    'XE': 1.8,
    'AU': 1.8,
    'LI': 1.8,
    '.': 1.8
    }

    two_char_elements = [el for el, r in atom_radii.items() if len(el) == 2]

    names = pdb.getNames()

    elements = []

    for name in names:
        element = ''
        for c in name:
            if not c.isdigit() and c != " ":
                element += c
        if element[:2] in two_char_elements:
            element = element[:2]
        else:
            element = element[0]
        elements.append(element)

    radii = []
    for element in elements:
        if element in atom_radii:
            r = atom_radii[element]
        else:
            r = atom_radii['.']
        radii.append(r)

    return(np.array(radii))


def generate_sphere_points(n):
    """
    Returns list of 3d coordinates of points on a sphere using the
    Golden Section Spiral algorithm.
    """
    points = []
    inc = np.pi * (3 - np.sqrt(5))
    offset = 2 / float(n)

    for k in xrange(int(n)):
        y = k * offset - 1 + (offset / 2)
        r = math.sqrt(1 - y*y)
        phi = k * inc
        points.append([np.cos(phi)*r, y, np.sin(phi)*r])

    return np.array(points)


def find_neighbor_indices_np(dist_sq, radii, probe, k):
    """
    Returns list of indices of atoms within probe distance to atom k
    leveraging numpy for accelerated calculations.  This implementation
    is about 10x faster than the original.
    """

    radius = radii[k] + probe + probe
    radii = (radii + radius) **2

    dist_sq = (dist_sq < radii)

    dist_sq[k] = False

    return dist_sq


def calculate_asa_np(pdb, probe, n_sphere_point=960):
    """
    Returns list of accessible surface areas of the atoms, using the probe
    and atom radius to define the surface.
    """
    sphere_points = generate_sphere_points(n_sphere_point)

    points = pdb.getCoords()

    radii = get_radii_from_prody_pdb(pdb)

    radii_plus_probe_sq = (radii + probe)**2

    dist_matrix_sq = cdist(points, points, 'sqeuclidean')

    const = 4.0 * math.pi / n_sphere_point

    areas = np.zeros(len(pdb))

    for i in xrange(0,len(pdb)):

        neighbor_indices = find_neighbor_indices_np(dist_matrix_sq[i,:], radii, probe, i)

        radius = probe + radii[i]

        point_i = points[i, :]

        test_sphere_points = sphere_points*radius + point_i

        neighbor_points = points[neighbor_indices, :]

        neighbor_radii_sq = radii_plus_probe_sq[neighbor_indices]

        diff_sq = cdist(test_sphere_points, neighbor_points, 'sqeuclidean')

        dist_test = (diff_sq < neighbor_radii_sq)

        accessible_points = np.sum(dist_test,1)

        area = (n_sphere_point - np.count_nonzero(accessible_points)) * const * radius * radius
        areas[i] = area

    return areas
