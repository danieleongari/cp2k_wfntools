from __future__ import print_function

from __future__ import absolute_import
import sys
import os
import scipy.io
import math
import numpy as np
from numpy.linalg import inv  #to obtain the inverse of a cell
import random
import copy  #make duplicate of a class with copy.deepcopy(class)
from atomic_data import atomic_valence_default_dict
from atomic_data import atomic_rad_UFF_dict
from collections import Counter  # to count the atom types
#from pprint import pprint  #for debug: prints all the attributs of an obj: pprint(vars(your_object))
from six.moves import range

# pylint: disable=too-many-arguments,too-many-locals,too-many-branches


class wfn(object):
    """ Wave function object, as stored in a .wfn file """

    def __init__(self,
                 natom=None,
                 nspin=None,
                 nao_tot=None,
                 nset_max=None,
                 nshell_max=None):
        #defined at the initialization
        self.natom = int(natom)  # number of atoms
        self.nspin = int(nspin)  # number of spin (i.e., 1 or 2)
        self.nao_tot = int(
            nao_tot)  # total number of atomic orbitals (=sum(self.nao))
        self.nset_max = int(
            nset_max)  # max number of sets in the basis set (!!!! always 1???)
        self.nshell_max = int(nshell_max)  # mak number of shells in each set
        #spin independent
        self.nset = [0 for _i in range(self.natom)
                     ]  # number of sets for the bs of each atom
        #if self.nset_max==1:
        #    self.nshell = [0 for i in range(self.natom)]                                # number of shells for the bs of each atom
        #    self.nao = [[0 for i in range(self.nshell_max)] for j in range(self.natom)] # number of atomic orbitals for each shell of each atom
        #else: # TODO: to be merged later
        self.nshell = [
            [0 for _i in range(self.nset_max)] for _j in range(self.natom)
        ]  # number of shells for the bs of each atom
        self.nao = [[[0 for _l in range(self.nshell_max)]
                     for _i in range(self.nset_max)]
                    for _j in range(self.natom)
                    ]  # number of atomic orbitals for each shell of each atom
        #spin dependent
        self.nmo = [0, 0]  # number of molecular orbitals
        self.nocc = [0, 0]  # number of occupied molecular orbitals
        self.nvirt = [0, 0]  # number of virtual molecular orbitals
        self.nel = [0, 0]  # number of electrons

    def initialize_lists(self, ispin):
        # Initialize (or reset to zero) alpha or beta molecular orbitals
        if ispin == 0:
            self.eigen = [[0 for _i in range(self.nmo[ispin])]
                          ]  # eigenvalues for each molecular orbital
            self.occup = [[0 for _i in range(self.nmo[ispin])]
                          ]  # occupancies for each molecular orbital
            self.coeff = [[[0 for _j in range(self.nao_tot)]
                           for _i in range(self.nmo[ispin])]
                          ]  #coefficient list for each moleculat orbital
        elif ispin == 1:
            self.eigen.append([0 for _i in range(self.nmo[ispin])])
            self.occup.append([0 for _i in range(self.nmo[ispin])])
            self.coeff.append([[0 for _j in range(self.nao_tot)]
                               for _i in range(self.nmo[ispin])])

    def add_nel(self, nel_geom, charge,
                multiplicity):  #default: charge=0, mult=1
        # Check that nel/mult/charge are consistent
        # If nspin=2 split into alpha and beta
        # If charge!=0, add/remove electrons
        # If multiplicity>1, rearrange alpha and beta electrons
        self.charge = charge
        self.mult = multiplicity
        if self.nspin == 1:
            if self.mult > 1:
                print("WARNING: nspin=1 but multiplicity>1! EXIT")
                sys.exit()
            self.nel[0] = int(nel_geom - self.charge)
            if self.nel[0] % 2 == 1:
                print("WARNING: nspin=1 but odd number of electrons! EXIT")
                sys.exit()
        elif self.nspin == 2:
            if (nel_geom - self.charge % 2 == 1) and (self.mult % 2 == 1):
                print(
                    "WARNING: both number of electrons and multiplicity are odd! EXIT"
                )
            self.nel[0] = int((nel_geom - self.charge) / 2 +
                              (self.mult - 1) / 2)
            self.nel[1] = int((nel_geom - self.charge) / 2 -
                              (self.mult - 1) / 2)


class mol(object):
    """ Molecule object, without the unit cell """

    def __init__(self, natom=None):
        # Initialize the molecule, known the number of atoms
        self.natom = natom  # number of atoms
        self.atom_type = [None for _i in range(self.natom)
                          ]  # atomic type: e.g., Cu1
        self.atom_xyz = [[None for _j in range(3)]
                         for _i in range(self.natom)]  # cartesian coordinates

    def compute_types(self):
        # Count and list the atom types: extract the element name if type = elementXX (with xx being a number)
        self.atom_element = [extract_element(x) for x in self.atom_type
                             ]  # atomic element: e.g., Cu
        self.ntype = len(set(self.atom_type))  # number of different atom types
        self.type_symbol = list(set(self.atom_type))  # list of types
        self.type_element = {}
        for x in self.type_symbol:
            self.type_element.update({
                x: extract_element(x)
            })  # list of elements for each type
        self.type_count = Counter(
            self.atom_type)  # count of atoms for each type
        self.bs = {}  # gaussian basis set
        self.pot = {}  # pseudopotential
        self.q = {}  # valence electrons

    def compute_nel_geom(self):
        # Count the number of electrons, according to the atom types
        self.nel = 0
        for element in self.atom_element:
            self.nel += atomic_valence_default_dict[
                element]  #TODO: user number of valence el

    def compute_fract(self, cell):
        # Given a cell, compute the fractional coordinates of the atoms
        self.atom_fract = np.zeros((self.natom, 3))
        for i in range(self.natom):
            for j in range(3):
                self.atom_fract[i][j] = self.atom_xyz[i][0] * cell.invmatrix[0][j] + \
                                        self.atom_xyz[i][1] * cell.invmatrix[1][j] + \
                                        self.atom_xyz[i][2] * cell.invmatrix[2][j]

    def overwrite_xyz_from_fract(self, cell):
        # Compute and overwrite the cartesian coordinates from the fractional ones
        for i in range(self.natom):
            for j in range(3):
                self.atom_xyz[i][j] = self.atom_fract[i][0] * cell.matrix[0][j] + \
                                      self.atom_fract[i][1] * cell.matrix[1][j] + \
                                      self.atom_fract[i][2] * cell.matrix[2][j]

    def add_bspot(self, kind, bs, pot):
        # If valence electron are specified, remove it and store it separately
        self.bs.update({kind: bs.split("-q")[0]})
        self.pot.update({kind: pot.split("-q")[0]})
        for val in bs, pot:
            if len(val.split("-q")) > 1:
                self.q.update({kind: int(val.split("-q")[1])})


class cell(object):
    """ Unit cell object, you should specify if it is read from lengths+angles or matrix"""

    def __init__(self):
        self.from_lengths_angles = False
        self.from_matrix = False
        self.length = np.zeros(3)
        self.angle_deg = np.zeros(3)
        self.angle_rad = np.zeros(3)
        self.matrix = np.zeros((3, 3))

    def compute_all_info(self):
        # Compute several thing after the unit cell is specified
        if self.from_lengths_angles:
            if hasattr(self, 'deg'):
                self.angle_rad = [math.radians(i) for i in self.angle_deg]
            else:
                self.angle_deg = [math.degrees(i) for i in self.angle_rad]
            self.matrix = np.zeros((3, 3))
            self.matrix[0][0] = self.length[0]
            self.matrix[0][1] = 0.0
            self.matrix[0][2] = 0.0
            self.matrix[1][0] = self.length[1] * math.cos(self.angle_rad[2])
            self.matrix[1][1] = self.length[1] * math.sin(self.angle_rad[2])
            self.matrix[1][2] = 0.0
            self.matrix[2][0] = self.length[2] * math.cos(self.angle_rad[1])
            self.matrix[2][1] = self.length[2] * (
                math.cos(self.angle_rad[0]) - math.cos(self.angle_rad[2]) *
                math.cos(self.angle_rad[1])) / math.sin(self.angle_rad[2])
            self.matrix[2][2] = self.length[2] * math.sqrt(
                1 - (math.cos(self.angle_rad[1]))**2 -
                ((math.cos(self.angle_rad[0]) -
                  math.cos(self.angle_rad[2]) * math.cos(self.angle_rad[1])) /
                 math.sin(self.angle_rad[2]))**2)
        if self.from_matrix:
            self.angle_rad = np.zeros(3)
            self.length[0] = math.sqrt(self.matrix[0][0] * self.matrix[0][0] +
                                       self.matrix[0][1] * self.matrix[0][1] +
                                       self.matrix[0][2] * self.matrix[0][2])
            self.length[1] = math.sqrt(self.matrix[1][0] * self.matrix[1][0] +
                                       self.matrix[1][1] * self.matrix[1][1] +
                                       self.matrix[1][2] * self.matrix[1][2])
            self.length[2] = math.sqrt(self.matrix[2][0] * self.matrix[2][0] +
                                       self.matrix[2][1] * self.matrix[2][1] +
                                       self.matrix[2][2] * self.matrix[2][2])
            self.angle_rad[0] = math.acos(
                (self.matrix[1][0] * self.matrix[2][0] + self.matrix[1][1] *
                 self.matrix[2][1] + self.matrix[1][2] * self.matrix[2][2]) /
                self.length[1] / self.length[2])  #alpha=B^C
            self.angle_rad[1] = math.acos(
                (self.matrix[0][0] * self.matrix[2][0] + self.matrix[0][1] *
                 self.matrix[2][1] + self.matrix[0][2] * self.matrix[2][2]) /
                self.length[0] / self.length[2])  #beta=A^C
            self.angle_rad[2] = math.acos(
                (self.matrix[0][0] * self.matrix[1][0] + self.matrix[0][1] *
                 self.matrix[1][1] + self.matrix[0][2] * self.matrix[1][2]) /
                self.length[0] / self.length[1])  #gamma=A^B
            self.angle_deg = [math.degrees(i) for i in self.angle_rad]
        self.invmatrix = inv(self.matrix)


def extract_element(element_type):
    """ Clear numberical values at the end of a string, for each string in a list """
    for x in element_type[::-1]:  # Reverse string
        if x.isnumeric():
            element_type = element_type[:-1]
    return element_type


def advreadline(inpfile):
    """ Avanced readline: return the splitted line or the string 'EoF' """
    line = inpfile.readline()
    if line == "":
        data = "EoF"
    else:
        data = line.split()
    return data


def read_wfn_file(wfn_file):
    """ Parser for the .wfn file, that returns the wfn object """
    inpfile = scipy.io.FortranFile(wfn_file, "r")
    natom, nspin, nao_tot, nset_max, nshell_max = inpfile.read_ints()
    w = wfn(natom, nspin, nao_tot, nset_max, nshell_max)
    nset = inpfile.read_ints()
    for i, val in enumerate(nset):
        w.nset[i] = val
    nshell = inpfile.read_ints()
    for i, val in enumerate(nshell):
        iatom, iset = divmod(i, w.nset_max)
        w.nshell[iatom][iset] = val
    nao = inpfile.read_ints()
    for i, val in enumerate(nao):
        iatom, isetshell = divmod(
            i, w.nset_max * w.nshell_max
        )  # this splits the N of atomic orbitals in one list for each atom
        iset, ishell = divmod(
            isetshell, w.nshell_max
        )  # this splits the N of atomic orbitals in one list for eash set
        w.nao[iatom][iset][ishell] = val
    for ispin in range(w.nspin):
        nmo, nocc, nvirt, nel = inpfile.read_ints()
        w.nmo[ispin] = nmo
        w.nocc[ispin] = nocc
        w.nvirt[ispin] = nvirt
        w.nel[ispin] = nel
        w.initialize_lists(ispin)
        eigen_and_occup = inpfile.read_reals()  #this line contains both
        for i in range(w.nmo[ispin]):
            w.eigen[ispin][i] = eigen_and_occup[i]
        for i in range(w.nmo[ispin]):
            w.occup[ispin][i] = eigen_and_occup[i + w.nmo[ispin]]
        for imo in range(w.nmo[ispin]):
            coeff = inpfile.read_reals()
            for i, val in enumerate(coeff):
                w.coeff[ispin][imo][i] = val
    inpfile.close()
    return w


def make_clean_wfn(w):
    """ Given a wfn, duplicate it and set to zero the coeffs """
    x = copy.deepcopy(w)
    for ispin in range(x.nspin):
        for i in range(x.nmo[ispin]):
            for j in range(x.nao_tot):
                x.coeff[ispin][i][j] = 0.0
    return x


def combine_wfn(a, b):  # noqa: MC0001
    """ Combine wfns a and b into a single ab wfn object """
    natom = a.natom + b.natom
    nspin = max(a.nspin, b.nspin)
    nao_tot = a.nao_tot + b.nao_tot
    nset_max = max(a.nset_max, b.nset_max)
    nshell_max = max(a.nshell_max, b.nshell_max)
    w = wfn(natom, nspin, nao_tot, nset_max, nshell_max)
    # combine and assign nset, for each atom
    nset = a.nset + b.nset
    for i, val in enumerate(nset):
        w.nset[i] = val
    nshell = a.nshell + b.nshell
    # combine and assign nshell, for each set, for each atom
    for iatom in range(w.natom):
        if iatom < a.natom:
            for iset in range(a.nset_max):
                w.nshell[iatom][iset] = a.nshell[iatom][iset]
        else:
            for iset in range(b.nset_max):
                w.nshell[iatom][iset] = b.nshell[iatom - a.natom][iset]
    # combine and assign nao, for each shell, for each set, for each atom
    for iatom in range(w.natom):
        if iatom < a.natom:
            for iset in range(a.nset_max):
                for ishell in range(a.nshell_max):
                    w.nao[iatom][iset][ishell] = a.nao[iatom][iset][ishell]
        else:
            for iset in range(b.nset_max):
                for ishell in range(b.nshell_max):
                    w.nao[iatom][iset][ishell] = b.nao[iatom -
                                                       a.natom][iset][ishell]
    # work with spin dependent properties
    if a.nspin == 1 and w.nspin == 2: a = makeopenshell(a)
    if b.nspin == 1 and w.nspin == 2: b = makeopenshell(b)
    for ispin in range(w.nspin):
        w.nmo[ispin] = a.nmo[ispin] + b.nmo[ispin]
        w.nocc[ispin] = a.nocc[ispin] + b.nocc[ispin]
        w.nvirt[ispin] = a.nvirt[ispin] + b.nvirt[
            ispin]  #not sure this is right!
        w.nel[ispin] = a.nel[ispin] + b.nel[ispin]
        w.initialize_lists(ispin)
        w.eigen[ispin] = a.eigen[ispin] + b.eigen[ispin]
        w.occup[ispin] = a.occup[ispin] + b.occup[ispin]
        #Note: you need to leave zeros for mixed AO/MO coeffs
        for imo in range(w.nmo[ispin]):
            if imo < a.nmo[ispin]:
                for iao in range(a.nao_tot):
                    w.coeff[ispin][imo][iao] = a.coeff[ispin][imo][iao]
            else:
                for iao in range(b.nao_tot):
                    w.coeff[ispin][imo][a.nao_tot + iao] = b.coeff[ispin][
                        imo - a.nmo[ispin]][iao]
    return w


def split_wfn(ab, A, B, a_charge, b_charge, a_mult, b_mult):  # noqa: MC0001
    """ Split the ab wfn into a and b wfns """
    nab = ab.natom
    na = A.natom
    nb = B.natom
    if nab != (na + nb):
        print("WARNING: in split_wfn the atoms in ab (from the wfn file)" +
              " are not the sum of A and B (from the geometry file)! EXIT")
        sys.exit()
    # Count the number of atomic orbitals for A and B and initialize the a and b wfns
    a_nao_tot = 0
    b_nao_tot = 0
    for i in range(nab):
        if i < na:
            for j in range(ab.nset_max):
                a_nao_tot += sum(ab.nao[i][j])
        else:
            for j in range(ab.nset_max):
                b_nao_tot += sum(ab.nao[i][j])
    #Note1: for semplicity I keep the same n_shell max even if the single wfn has less
    #Note2: if nspin=2 you want to keep uks for both in the CP calculation
    a = wfn(na, ab.nspin, a_nao_tot, ab.nset_max, ab.nshell_max)
    b = wfn(nb, ab.nspin, b_nao_tot, ab.nset_max, ab.nshell_max)
    for i in range(nab):
        if i < na:
            a.nset[i] = ab.nset[i]
            a.nshell[i] = ab.nshell[i]
            a.nao[i] = ab.nao[i]
        else:
            b.nset[i - na] = ab.nset[i]
            b.nshell[i - na] = ab.nshell[i]
            b.nao[i - na] = ab.nao[i]
    # Count the electrons and the occupied orbitals, and initialize the lists
    A.compute_nel_geom()
    B.compute_nel_geom()
    a.add_nel(A.nel, a_charge, a_mult)
    b.add_nel(B.nel, b_charge, b_mult)
    if ab.nspin == 1: el_per_orb = 2
    elif ab.nspin == 2: el_per_orb = 1
    for w in [a, b]:
        for ispin in range(ab.nspin):
            w.nmo[ispin] = int(w.nel[ispin] / el_per_orb)
            w.nocc[ispin] = w.nmo[ispin]  #all the orbitals are occupied
            w.nvirt[ispin] = 0  #no virtual orbitals
            w.nel[ispin] = a.nel[ispin]
            w.initialize_lists(ispin)
            for i in range(w.nmo[ispin]):
                w.eigen[ispin][i] = 0.0  #left to zero
                w.occup[ispin][i] = float(el_per_orb)
    # Split the coefficients
    for ispin in range(ab.nspin):
        for i in range(a.nmo[ispin]):
            for j in range(a.nao_tot):
                a.coeff[ispin][i][j] = ab.coeff[ispin][i][j]
        for i in range(b.nmo[ispin]):
            for j in range(b.nao_tot):
                b.coeff[ispin][i][j] = ab.coeff[ispin][i + a.nmo[ispin]][
                    j + a.nao_tot]
    return a, b


def write_wfn_file(wfn_object, wfn_file):
    """ Print a wfn object as a binary fortran file """
    outfile = scipy.io.FortranFile(wfn_file, "w")
    w = wfn_object
    outfile.write_record(
        np.array([w.natom, w.nspin, w.nao_tot, w.nset_max, w.nshell_max],
                 dtype=np.int32))
    outfile.write_record(np.array(w.nset, dtype=np.int32))
    outfile.write_record(np.array(w.nshell, dtype=np.int32))
    outfile.write_record(np.array(w.nao, dtype=np.int32))
    for ispin in range(w.nspin):
        outfile.write_record(
            np.array(
                [w.nmo[ispin], w.nocc[ispin], w.nvirt[ispin], w.nel[ispin]],
                dtype=np.int32))
        outfile.write_record(
            np.array(
                np.concatenate([w.eigen[ispin], w.occup[ispin]]),
                dtype=np.dtype("f8")))
        for imo in range(w.nmo[ispin]):
            outfile.write_record(
                np.array(w.coeff[ispin][imo], dtype=np.dtype("f8")))
    outfile.close()
    return


def write_fwfn_file(wfn_object, fwfn_file):
    """ Print a wfn object as a formatted text file """
    outfile = open(fwfn_file, "w")
    w = wfn_object
    print("Number of atoms: %d" % w.natom, file=outfile)
    print("Number of spin: %d" % w.nspin, file=outfile)
    print(
        "Total number of atomic orbitals (AOs): %d" % w.nao_tot, file=outfile)
    print("Maximum number of sets: %d" % w.nset_max, file=outfile)
    print("Maximum number of shells: %d" % w.nshell_max, file=outfile)
    print("Number of sets for each atom:", end="", file=outfile)
    for val in w.nset:
        print(" %d" % val, end="", file=outfile)
    print("", file=outfile)
    print("Number of shells for each [atom] and set: ", end="", file=outfile)
    for iset in range(len(w.nshell)):
        print("[ ", end="", file=outfile)
        for val in w.nshell[iset]:
            print("%d " % val, end="", file=outfile)
        print("] ", end="", file=outfile)
    print("", file=outfile)
    print(
        "Number of atomic orbitals for each atom, [set] and shell:",
        file=outfile)
    for i, thelist in enumerate(w.nao):
        print("atom %d:\t" % i, end="", file=outfile)
        for iset in range(len(thelist)):
            print("[ ", end="", file=outfile)
            for val in thelist[iset]:
                print("%d " % val, end="", file=outfile)
            print("] ", end="", file=outfile)
        print("", file=outfile)
    for ispin in range(w.nspin):
        print("", file=outfile)
        print("Information for spin %d:" % ispin, file=outfile)
        print(
            "Number of molecular orbitals (MOs): %d" % w.nmo[ispin],
            file=outfile)
        print("Number of occupied MOs: %d" % w.nocc[ispin], file=outfile)
        print("Number of virtual MOs: %d" % w.nvirt[ispin], file=outfile)
        print("Number of electrons: %d" % w.nel[ispin], file=outfile)
        print("Eigenvalue for each MO: [", end="", file=outfile)
        for val in w.eigen[ispin]:
            print(" %.3f" % val, end="", file=outfile)
        print(" ]", file=outfile)
        print("Electron occupancy for each MO: [", end="", file=outfile)
        for val in w.occup[ispin]:
            print(" %.3f" % val, end="", file=outfile)
        print(" ]", file=outfile)
        print("Coefficients for the MO/AO matrix:", file=outfile)
        for i, thelist in enumerate(w.coeff[ispin]):
            print("MO %d:\t[" % i, end="", file=outfile)
            for j in thelist:
                print(" %+.3f" % j, end="", file=outfile)
            print(" ]", file=outfile)
    outfile.close()
    return


def makeopenshell(v):
    """ Makes a restricted (e.g., nspin=1) wfn open shell (e.g., nspin=2) """
    if v.nspin == 2:
        return v
    elif v.nspin == 1:
        w = wfn(v.natom, 2, v.nao_tot, v.nset_max, v.nshell_max)
        w.nset = v.nset
        w.nshell = v.nshell
        w.nao = v.nao
        for ispin in range(2):
            w.nmo[ispin] = v.nmo[0]
            w.nocc[ispin] = v.nocc[0]
            w.nvirt[ispin] = v.nvirt[0]
            w.nel[ispin] = v.nel[0] / 2
            w.initialize_lists(ispin)
            w.eigen[ispin] = [x / 2. for x in v.eigen[0]]
            w.occup[ispin] = [x / 2. for x in v.occup[0]]
            w.coeff[ispin] = v.coeff[0]
        return w


def read_xyz_file(xyz_file):
    """ Read an .xyz file, containing only one fragment. Neglect _A and _B labels if present """
    inpfile = open(xyz_file, "r")
    data = advreadline(inpfile)
    natom = int(data[0])
    M = mol(natom)
    advreadline(inpfile)
    for i in range(natom):
        data = advreadline(inpfile)
        if data[0][-2:] == "_A" or data[0][-2:] == "_B":
            M.atom_type[i] = data[0][:-2]
        else:
            M.atom_type[i] = data[0]
        M.atom_xyz[i][0] = float(data[1])
        M.atom_xyz[i][1] = float(data[2])
        M.atom_xyz[i][2] = float(data[3])
    inpfile.close()
    M.compute_types()
    return M


def read_xyzlabel_file(xyz_file):
    """ Read an .xyz labelled file containing both A and B fragments.
    NB: it does not matter the order of A and B in the file """
    inpfile = open(xyz_file, "r")
    nab = int(advreadline(inpfile)[0])
    advreadline(inpfile)
    na = 0
    nb = 0
    for _f in range(nab):
        atom_label = advreadline(inpfile)[0][-2:]
        if atom_label == "_A":
            na += 1
        elif atom_label == "_B":
            nb += 1
        else:
            print(
                "WARNING: AB.xyz does not contain labels. Use labelAB first! EXIT"
            )
            sys.exit()
    if not nab == (na + nb):
        print("WARNING: the count of A and B atoms in AB.xyz is weird! EXIT")
        sys.exit()
    inpfile.close()  #TODO: use rewind
    inpfile = open(xyz_file, "r")
    advreadline(inpfile)
    advreadline(inpfile)
    outfileA = open("tmp_A.xyz", "w")
    outfileB = open("tmp_B.xyz", "w")
    print(na, file=outfileA)
    print(nb, file=outfileB)
    print("Fragment A", file=outfileA)
    print("Fragment B", file=outfileB)
    for _i in range(na):
        line = inpfile.readline()
        print(line, end="", file=outfileA)
    for _i in range(nb):
        line = inpfile.readline()
        print(line, end="", file=outfileB)
    outfileA.close()
    outfileB.close()
    A = read_xyz_file("tmp_A.xyz")
    B = read_xyz_file("tmp_B.xyz")
    os.remove("tmp_A.xyz")
    os.remove("tmp_B.xyz")
    return A, B, na, nb


def write_xyz_file(xyz_file, A, B, label):
    """ Write .xyz file, given two fragments. With label True prints 'Cu1_A', with False 'Cu' """
    outfile = open(xyz_file, "w")
    print("%d" % (A.natom + B.natom), file=outfile)
    print("Printed using cp2k_wfntool", file=outfile)
    for i in range(A.natom):
        if label: atom_label = A.atom_type[i] + "_A"
        else: atom_label = A.atom_element[i]
        print(
            "%s %10.5f %10.5f %10.5f" % (atom_label, A.atom_xyz[i][0],
                                         A.atom_xyz[i][1], A.atom_xyz[i][2]),
            file=outfile)
    for i in range(B.natom):
        if label: atom_label = B.atom_type[i] + "_B"
        else: atom_label = B.atom_element[i]
        print(
            "%s %10.5f %10.5f %10.5f" % (atom_label, B.atom_xyz[i][0],
                                         B.atom_xyz[i][1], B.atom_xyz[i][2]),
            file=outfile)
    outfile.close()
    return


def read_cell_file(cell_file):
    """ Read the &CELL section in CP2K format.

    TODO: this function can be more and more flexible in the parsing 
    """
    c = cell()
    nlines = sum(1 for line in open(cell_file))
    inpfile = open(cell_file, "r")
    for _l in range(nlines):
        data = advreadline(inpfile)
        if len(data) >= 4 and data[0] == "A":
            c.from_matrix = True
            if data[1][0] == "[":
                c.matrix[0] = data[2:5]
            else:
                c.matrix[0] = data[1:4]
        if len(data) >= 4 and data[0] == "B":
            if data[1][0] == "[":
                c.matrix[1] = data[2:5]
            else:
                c.matrix[1] = data[1:4]
        if len(data) >= 4 and data[0] == "C":
            if data[1][0] == "[":
                c.matrix[2] = data[2:5]
            else:
                c.matrix[2] = data[1:4]
        if len(data) >= 4 and data[0] == "ABC":
            c.from_lengths_angles = True
            if data[1][0] == "[":
                c.length = data[2:5]
            else:
                c.length = data[1:4]
        if len(data) >= 4 and data[0] == "ALPHA_BETA_GAMMA":
            if data[1][0] == "[":
                c.angle = data[2:5]
            else:
                c.angle = data[1:4]
    c.compute_all_info()
    return c


def read_kind_file(kind_file, A, B):  # noqa: MC0001
    """ Read CP2K's &KIND section.

    Use A or B = None to use only one fragment, so that labels are not needed """
    inpfile = open(kind_file, "r")
    data = None
    while data != "EoF":
        data = advreadline(inpfile)
        if len(data) >= 2 and data[0] == "&KIND":
            kind = data[1]
            # Decide if the kind belongs to A or B
            if kind[-2:] == "_A" or B is None:
                M = A
            if kind[-2:] == "_B" or A is None:
                M = B
            if kind[-2:] == "_A" or kind[-2:] == "_B":
                kind = kind[:-2]
            # Check if the kind is in the molecule, otherwise skip it
            kind_in_mol = False
            if kind in M.type_symbol:
                kind_in_mol = True
            # Look for BASIS_SET and PSEUDO
            bs = None
            pot = None
            while data != "EoF" and kind_in_mol:
                data = advreadline(inpfile)
                if len(data) >= 1 and data[0] == "&END":
                    break
                if len(data) >= 2 and data[0] == "BASIS_SET":
                    bs = data[1]
                if len(data) >= 2 and data[0] == "POTENTIAL":
                    pot = data[1]
                if len(data) >= 1 and data[0] == "GHOST":
                    print(
                        "WARNING: reading .kind and not expecting a GHOST here! EXIT"
                    )
                    sys.exit()
            if kind_in_mol:
                M.add_bspot(kind, bs, pot)
    # Final check
    for M in A, B:
        if M is not None:
            if len(M.bs) != M.ntype:
                print(
                    "WARNING: the basis set is read from a kind file, but has a different number of types!"
                )
    return


def write_kind_file(kind_file, A, B, A_is_ghost, B_is_ghost, bs_choice,
                    pot_choice):
    """ Write CP2K's &KIND section """
    outfile = open(kind_file, "w")
    for M in A, B:
        for x in M.type_symbol:
            if M == A:
                label = "_A"
                is_ghost = A_is_ghost
            if M == B:
                label = "_B"
                is_ghost = B_is_ghost
            e = M.type_element[x]
            q = atomic_valence_default_dict[e]
            if bs_choice == "read" and (x in M.q):
                q = M.q[x]
            print("    &KIND %s" % (x + label), file=outfile)
            print("      ELEMENT %s" % e, file=outfile)
            if bs_choice == "read":
                bs = M.bs[x] + "-q" + str(q)
            else:
                bs = bs_choice + "-q" + str(q)
            print("      BASIS_SET %s" % bs, file=outfile)
            if is_ghost:
                print("      GHOST", file=outfile)
            else:
                if pot_choice == "read":
                    pot = M.pot[x] + "-q" + str(q)
                else:
                    pot = pot_choice + "-q" + str(q)
                print("      POTENTIAL %s" % pot, file=outfile)
            print("    &END", file=outfile)
    outfile.close()
    return


def rotate_rand(B):
    """ Make a rotated B molecule, randomly and fixed the centroid

    Note: if B is very long, and the UC not cubic, strange things may happen! 
    """
    M = copy.deepcopy(B)
    # Translate the centroid to the origin
    centroid_fract = np.zeros(3)
    for k in range(3):
        centroid_fract[k] = np.mean(
            [M.atom_fract[i][k] for i in range(M.natom)])
    for i in range(M.natom):
        for k in range(3):
            M.atom_fract[i][j] = M.atom_fract[i][j] - centroid_fract[k]
    # ********* TODO : rotation ***********************************
    M.overwrite_xyz_from_fract(cell)
    return M


def translate_rand(B, cell):
    """ Make a translated B molecule """
    M = copy.deepcopy(B)
    rand_transl_fract = [random.random(), random.random(), random.random()]
    for i in range(M.natom):
        for k in range(3):
            M.atom_fract[i][k] = M.atom_fract[i][k] + rand_transl_fract[k]
            # get in the cell if the coordinate is negative or more than unit.
            # Note: this algorithm works well for all the cases, 0<a<1, a>1 and a<0
            M.atom_fract[i][j] -= math.floor(M.atom_fract[i][j])
    M.overwrite_xyz_from_fract(cell)
    return M


def has_overlap(A, B, cell, scaledrad, debug):
    """ Checks if A and B fragments are overlapping """
    DEBUG = debug
    ovlp_found = False
    dist_fract = np.zeros(3)
    dist_xyz = np.zeros(3)
    if DEBUG: mindist_store = 100
    for i in range(A.natom):
        for j in range(B.natom):
            for k in range(3):
                dist_fract[k] = A.atom_fract[i][k] - B.atom_fract[j][k]
                dist_fract[k] = dist_fract[k] - int(round(dist_fract[k]))  #pbc
            for k in range(3):
                dist_xyz[k] = dist_fract[0] * cell.matrix[0][k] + \
                              dist_fract[1] * cell.matrix[1][k] + \
                              dist_fract[2] * cell.matrix[2][k]
            mindist = math.sqrt(sum([x**2 for x in dist_xyz]))
            vdwdist = scaledrad * (atomic_rad_UFF_dict[A.atom_type[i]] +
                                   atomic_rad_UFF_dict[B.atom_type[j]])
            if DEBUG and mindist < mindist_store:
                i_store = i
                j_store = j
                dist_fract_store = dist_fract
                dist_xyz_store = dist_xyz
                mindist_store = mindist
                vdwdist_store = vdwdist
            if mindist < vdwdist:
                ovlp_found = True
                break
        if ovlp_found:
            break
    if DEBUG and not ovlp_found:  #TODO: improve the printing
        print()
        print(i_store, j_store)
        print(dist_fract_store)
        print(dist_xyz_store)
        print(mindist_store)
        print(vdwdist_store)
    return ovlp_found
