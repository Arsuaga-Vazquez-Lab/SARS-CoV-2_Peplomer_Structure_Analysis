"""
A collection of functions to use in analyzing the SARS-CoV2 spike protein

Throughout this file, groups of points will be stored mostly as matrices with
each row representing a point, and each column representing a dimension (x, y, or z)
"""
import math
import numpy as np
from io import StringIO

def backbone(protein_pdb: StringIO, chain="all chains", CA_only=True):
    # Uses a PDB file to get the coordinates of atoms in the backbone of a protein chain
    pdb_data = protein_pdb.read()
    atoms = []
    backbone_atoms = ['CA'] if CA_only else ['CA', 'C', 'N']
    for line in pdb_data.split('\n'):
        # By default, look at all chains, unless just one is specified
        if line[0:4] == 'ATOM' and (chain is "all chains" or line[21] == chain) and line[13:16].rstrip(' ') in backbone_atoms:
            residue_num = int(line[24:26])
            x = float(line[31:38].strip(' '))
            y = float(line[39:46].strip(' '))
            z = float(line[47:54].strip(' '))
            atoms.append([x, y, z])
    return np.matrix(atoms)

def dump_beads(*components, file_name='beads.txt') -> None:
    # Given a list of atom coordinates, creates a file for knotplot to open
    # Find centroid of all chains, in order to center it
    centroid = np.array([[0., 0., 0.]])
    total_atoms = 0
    for component in components:
        for atom in component:
            centroid += atom[0]
            total_atoms += 1
    centroid /= total_atoms
    for component in components:
        for atom in component:
            atom -= centroid
    knotplot_data = ''
    components = [component.tolist() for component in components]
    for component_num in range(len(components)):
        knotplot_data += f'Component {component_num + 1} of {len(components)}:\n'
        component = components[component_num]
        for atom in component:
            knotplot_data += ' '.join([str(coord) for coord in atom]) + '\n'
    with open(file_name, 'w+') as knotplot_file:
        knotplot_file.write(knotplot_data)

def rmsd(chain1, chain2):
    # TODO: create functions for just getting the centroid or rotation matrix
    n = len(chain1)
    if n != len(chain2): raise Error
    chain1 = np.matrix(chain1)
    chain2 = np.matrix(chain2)
    # Center each group of points by subtracting their respective centroids
    chain1 = chain1 - np.ones((n, n)) * chain1 / n
    chain2 = chain2 - np.ones((n, n)) * chain2 / n
    # Use Kabsch algorithm to line up chain1 and chain2 in order to minimize RMSD
    # The Kabsch algo finds the rotation matrix that best maps chain1 onto chain2
    # https://en.wikipedia.org/wiki/Kabsch_algorithm
    covariance_matrix = chain1.T * chain2
    U, S, VT = np.linalg.svd(covariance_matrix)
    d = np.linalg.det(U * VT)
    rotation_matrix = U * np.diag([1,1,d]) * VT
    chain1 *= rotation_matrix
    total_square_diff = sum([row * row.T for row in chain1 - chain2])[0, 0]
    return math.sqrt(total_square_diff / n)

def unit_vec(vector):
    if np.linalg.norm(vector) == 0:
        raise ZeroDivisionError
    return vector / np.linalg.norm(vector)

def sign(num): return 1 if num > 0 else -1 if num < 0 else 0

def ACN_and_space_writhe(chain) -> (float, float):
    # calculates ACN & space writhe of a piecewise-linear path, which is given as a matrix of points
    # returns ACN and space writhe together, as a tuple of 2 floats
    # This formula, as well as all the variable names, is taken from:
    # https://en.wikipedia.org/wiki/Writhe#Numerically_approximating_the_Gauss_integral_for_writhe_of_a_curve_in_space
    acn = 0
    space_writhe = 0
    # Iterate through every pair of line segments
    for i in range(1, len(chain)):
        for j in range(1, i - 1):
            # p1, p2, p3, & p4 are endpoints of line segments i & j
            p1 = chain[i - 1]
            p2 = chain[i]
            p3 = chain[j - 1]
            p4 = chain[j]
            n1 = unit_vec(np.cross(p3 - p1, p4 - p1))
            n2 = unit_vec(np.cross(p4 - p1, p4 - p2))
            n3 = unit_vec(np.cross(p4 - p2, p3 - p2))
            n4 = unit_vec(np.cross(p3 - p2, p3 - p1))
            omega_star = math.asin(n1.dot(n2.T)[0, 0]) + \
                         math.asin(n2.dot(n3.T)[0, 0]) + \
                         math.asin(n3.dot(n4.T)[0, 0]) + \
                         math.asin(n4.dot(n1.T)[0, 0])
            omega_star /= 2 * math.pi
            omega_signed = omega_star * sign(np.cross(p4 - p3, p2 - p1).dot((p3 - p1).T)[0, 0])
            acn += omega_star
            space_writhe += omega_signed
    return (acn, space_writhe)

def rog(chain) -> float:
    # given a cluster of points (as an n by 3 matrix), calculates the radius of gyration
    # For now, assumes all points have the same weight
    moment = sum(chain)
    centroid = moment / len(chain)
    sum_of_distances_squared = sum([pow(np.linalg.norm(point - centroid), 2) for point in chain])
    return math.sqrt(sum_of_distances_squared / len(chain))
