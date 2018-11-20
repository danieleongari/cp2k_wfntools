#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import

import sys
import os
import argparse
from argparse import RawTextHelpFormatter  #needed to go next line in the help text
from six.moves import range
#from pprint import pprint  #for debug: prints all the attributs of an obj: pprint(vars(your_object))

from atomic_data import atomic_symbol
from utilities import mol
from utilities import read_wfn_file
from utilities import combine_wfn, split_wfn, make_clean_wfn
from utilities import write_fwfn_file, write_wfn_file
from utilities import read_xyz_file, read_xyzlabel_file
from utilities import read_cell_file, read_kind_file
from utilities import write_xyz_file, write_kind_file
from utilities import rotate_rand, translate_rand, has_overlap

parser = argparse.ArgumentParser(
    description="Program to parse and combine CP2K's .wfn wavefunction files",
    formatter_class=RawTextHelpFormatter)

parser.add_argument(
    "function",
    choices=[
        "parse", "combine", "cp", "makeAB", "labelAB", "swapAB", "check",
        "debug_kind", "debug_split"
    ],
    help="Choose the function to execute.\n" + "parse: print formatted wfn\n" +
    "combine: combine A and B geometries and wfn\n" +
    "cp: print files for a CounterPoise correction calculation\n" +
    "makeAB: insert B in A, in a random position\n" +
    "labelAB: given unlabelled AB.xyz, make _A and _B labels\n" +
    "swapAB: given labelled AB.xyz, swap _A with _B\n" +
    "check: check that a given geometry and wfn have the same atoms and types")

parser.add_argument(
    "-a",
    "--wfna",
    type=str,
    dest="wfna",
    help="Binary .wfn function file for fragment A")

parser.add_argument(
    "-b",
    "--wfnb",
    type=str,
    dest="wfnb",
    default=None,
    help="Binary .wfn function file for fragment A")

parser.add_argument(
    "-ab",
    "--wfnab",
    type=str,
    dest="wfnab",
    default=None,
    help="Binary .wfn function file for both fragments together")

parser.add_argument(
    "-A", "--geoa", type=str, dest="geoa", help="Geometry of fragment A")

parser.add_argument(
    "-B",
    "--geob",
    type=str,
    dest="geob",
    default=None,
    help="Geometry of fragment B")

parser.add_argument(
    "-AB",
    "--geoab",
    type=str,
    dest="geoab",
    default=None,
    help="Geometry of both fragments together (NB: should contain labels!)")

parser.add_argument(
    "-kindA",
    "--kindfileA",
    type=str,
    dest="kindfilea",
    default=None,
    help=
    "File containing the BS and PSEUDO in CP2K &KIND format, for fragment A")

parser.add_argument(
    "-kindB",
    "--kindfileB",
    type=str,
    dest="kindfileb",
    default=None,
    help=
    "File containing the BS and PSEUDO in CP2K &KIND format, for fragment B")

parser.add_argument(
    "-kindAB",
    "--kindfileAB",
    type=str,
    dest="kindfileab",
    default=None,
    help=
    "File containing the BS and PSEUDO in CP2K &KIND format, for both fragments"
)

parser.add_argument(
    "-cell",
    "--unitcelldimensions",
    type=str,
    dest="cell",
    default=None,
    help="File containing the unit cell in CP2K format, &CELL")

parser.add_argument(
    "-o",
    "--outfilename",
    type=str,
    dest="outfilename",
    default="wfntools_output",
    help="Name of the output file(s)")

parser.add_argument(
    "-bs",
    "--basisiset",
    type=str,
    dest="bs",
    default="DZVP-MOLOPT-SR-GTH",
    help="CP2K basis set choice for the .kind file.\n" +
    "Use 'read' if these are read with the -kind option.")

parser.add_argument(
    "-pot",
    "--potential",
    type=str,
    dest="pot",
    default="GTH-PBE",
    help="CP2K pseudopotential choice for the .kind file\n" +
    "Use 'read' if these are read with the -kind option.")

parser.add_argument(
    "-nA",
    "--numberofatomsinA",
    type=int,
    dest="na",
    default=None,
    help="Number of atoms of fragment A")

parser.add_argument(
    "-nB",
    "--numberofatomsinB",
    type=int,
    dest="nb",
    default=None,
    help="Number of atoms of fragment B")

parser.add_argument(
    "-qA",
    "--chargeofA",
    type=int,
    dest="qa",
    default=0,
    help="Charge of fragment A")

parser.add_argument(
    "-qB",
    "--chargeofB",
    type=int,
    dest="qb",
    default=0,
    help="Charge of fragment B")

parser.add_argument(
    "-multA",
    "--multiplicityofA",
    type=int,
    dest="multa",
    default=1,
    help="Multiplicity of fragment A")

parser.add_argument(
    "-multB",
    "--multiplicityofB",
    type=int,
    dest="multb",
    default=1,
    help="Multiplicity of fragment B")

parser.add_argument(
    "-swapAB",
    "--swapAwithB",
    action="store_true",
    dest="swapab",
    default=False,
    help="In the .xyz file read B first and A later")

parser.add_argument(
    "-sbs",
    "--singlebasisset",
    action="store_true",
    dest="sbs",
    default=False,
    help="In a CounterPoise calculation, print also Aa and Bb")

parser.add_argument(
    "-maxfail",
    "--maxfailures",
    type=int,
    dest="maxfail",
    default=100,
    help=
    "When using makeAB, the program is stopped after maxfail \
    overlapping configurations"
)

parser.add_argument(
    "-srad",
    "--scaledradius",
    type=float,
    dest="srad",
    default=1.0,
    help=
    "When using makeAB, scaling factor for the atomic radii, \
    to evaluate the overlap"
)

parser.add_argument(
    "-nout",
    "--numberofoutputAB",
    type=int,
    dest="nout",
    default=1,
    help="When using makeAB, number of printed output combinations of A and B")

parser.add_argument(
    "-pf",
    "--printformatted",
    action="store_true",
    dest="printformatted",
    default=False,
    help="When a .wfn file is printed, it prints also a formatted .fwfn file")

parser.add_argument(
    "-pmd",
    "--printmindist",
    action="store_true",
    dest="printmindist",
    default=False,
    help="When using makeAB, prints the minimum A-B distance")

args = parser.parse_args()

if args.function == "parse":
    # Note: any -a -b -ab tag can be used to specify the wfn, but not more than one
    if args.wfna is None and args.wfnb is None and args.wfnab is None:
        print("WARNING: wfn file not specified! EXIT")
        sys.exit()
    elif args.wfna is not None and args.wfnb is None and args.wfnab is None:
        wfnfile = args.wfna
    elif args.wfna is None and args.wfnb is not None and args.wfnab is None:
        wfnfile = args.wfnb
    elif args.wfna is None and args.wfnb is None and args.wfnab is not None:
        wfnfile = args.wfnab
    else:
        print("WARNING: more than one wfn file specified! EXIT")
        sys.exit()
    if not os.path.isfile(wfnfile):
        print("WARNING: wfn file not existing! EXIT")
        sys.exit()
    print("Converting %s binary file in %s.fwfn formatted file." %
          (wfnfile, args.outfilename))
    a = read_wfn_file(wfnfile)
    write_fwfn_file(a, args.outfilename + ".fwfn")

if args.function == "combine":
    #if wfargs.na and wfnb are provided, combine the two wfn and print fwfn if asked
    if (args.wfna is None or args.wfnb is None):
        print("WARNING: wavefunctions for A or B were not provided!")
    if (args.wfna is not None and args.wfnb is not None):
        print("Converting wfns %s and %s into %s.wfn" % (args.wfna, args.wfnb,
                                                         args.outfilename))
        a = read_wfn_file(args.wfna)
        b = read_wfn_file(args.wfnb)
        ab = combine_wfn(a, b)
        write_wfn_file(ab, args.outfilename + ".wfn")
        if args.printformatted:
            write_fwfn_file(ab, args.outfilename + ".fwfn")
    #if geoa and geob are provided, combine the two geometries and make the kind file
    if (args.geoa is None or args.geob is None):
        print("WARNING: geometries for A or B were not provided!")
    if (args.geoa is not None and args.geob is not None):
        print("Converting geometries %s and %s into %s_label.xyz" %
              (args.geoa, args.geob, args.outfilename))
        A = read_xyz_file(args.geoa)
        B = read_xyz_file(args.geob)
        if args.kindfilea is not None:
            read_kind_file(args.kindfilea, A, None)
        if args.kindfileb is not None:
            read_kind_file(args.kindfileb, None, B)
        write_xyz_file(args.outfilename + "_nolabel.xyz", A, B, label=False)
        write_xyz_file(args.outfilename + "_label.xyz", A, B, label=True)
        write_kind_file(args.outfilename + ".kind", A, B, False, False,
                        args.bs, args.pot)

if args.function == "cp":
    if args.geoab is None or not os.path.isfile(args.geoab):
        print("WARNING: -AB geometry not provided or not existing! EXIT")
        sys.exit()
    A, B, na, nb = read_xyzlabel_file(args.geoab)
    if args.kindfileab is not None:
        read_kind_file(args.kindfileab, A, B)
    if args.wfnab is None or not os.path.isfile(args.wfnab):
        print("WARNING: -ab wave function not provided or not existing! EXIT")
        sys.exit()
    ab = read_wfn_file(args.wfnab)
    a, b = split_wfn(ab, A, B, args.qa, args.qb, args.multa, args.multb)
    a_clean = make_clean_wfn(a)
    b_clean = make_clean_wfn(b)

    # Print Aab
    a_ab = combine_wfn(a, b_clean)
    write_xyz_file(args.outfilename + "_Aab.xyz", A, B, label=True)
    write_kind_file(args.outfilename + "_Aab.kind", A, B, False, True, args.bs,
                    args.pot)
    write_wfn_file(a_ab, args.outfilename + "_Aab.wfn")
    if args.printformatted:
        write_fwfn_file(a_ab, args.outfilename + "_Aab.fwfn")
    # Print Bab
    b_ab = combine_wfn(a_clean, b)
    write_xyz_file(args.outfilename + "_Bab.xyz", A, B, label=True)
    write_kind_file(args.outfilename + "_Bab.kind", A, B, True, False, args.bs,
                    args.pot)
    write_wfn_file(b_ab, args.outfilename + "_Bab.wfn")
    if args.printformatted:
        write_fwfn_file(b_ab, args.outfilename + "_Bab.fwfn")
    if args.sbs:
        # Initialize an empty fragment
        X = mol(0)
        X.compute_types()
        # Print Aa
        write_xyz_file(args.outfilename + "_Aa.xyz", A, X, label=True)
        write_kind_file(args.outfilename + "_Aa.kind", A, X, False, False,
                        args.bs, args.pot)
        write_wfn_file(a, args.outfilename + "_Aa.wfn")
        if args.printformatted:
            write_fwfn_file(a, args.outfilename + "_Aa.fwfn")
        # Print Bb
        write_xyz_file(args.outfilename + "_Bb.xyz", A, X, label=True)
        write_kind_file(args.outfilename + "_Bb.kind", A, X, False, False,
                        args.bs, args.pot)
        write_wfn_file(b, args.outfilename + "_Bb.wfn")
        if args.printformatted:
            write_fwfn_file(a, args.outfilename + "_Bb.fwfn")

if args.function == "labelAB":
    if args.geoab is None or not os.path.isfile(args.geoab):
        print("WARNING: -AB geometry not provided or not existing! EXIT")
        sys.exit()
    # Read the number of atoms an check if everything makes sense: na+nb=nab
    inpfile = open(args.geoab, "r")
    nab = int(inpfile.readline().split()[0])
    if args.na is None and args.nb is None:
        print("WARNING: number of atoms for A (or B) not provided! EXIT")
        sys.exit()
    if args.na is not None and args.nb is not None:
        if nab != (args.na + args.nb):
            print(
                "WARNING: number of atoms for A and B both specified, \
                but they don't match with nAB! EXIT"
            )
            sys.exit()
    if args.na is None:
        args.na = nab - args.nb
    else:
        args.nb = nab - args.na
    if (args.na <= 0) or (args.nb <= 0):
        print("WARNING: something is wrong with the number of atoms! EXIT")
        sys.exit()

    # Print two temporary files for A and B geometries
    outfileA = open("tmp_A.xyz", "w")
    outfileB = open("tmp_B.xyz", "w")
    print(args.na, file=outfileA)
    print(args.nb, file=outfileB)
    print("Fragment A", file=outfileA)
    print("Fragment B", file=outfileB)
    junk = inpfile.readline()
    if not args.swapab:
        for i in range(args.na):
            line = inpfile.readline()
            print(line, end="", file=outfileA)
        for i in range(args.nb):
            line = inpfile.readline()
            print(line, end="", file=outfileB)
    else:
        for i in range(args.nb):
            line = inpfile.readline()
            print(line, end="", file=outfileB)
        for i in range(args.na):
            line = inpfile.readline()
            print(line, end="", file=outfileA)
    outfileA.close()
    outfileB.close()
    A = read_xyz_file("tmp_A.xyz")
    B = read_xyz_file("tmp_B.xyz")
    os.remove("tmp_A.xyz")
    os.remove("tmp_B.xyz")
    write_xyz_file(args.outfilename + "_label.xyz", A, B, label=True)

if args.function == "swapAB":
    if args.geoab is None or not os.path.isfile(args.geoab):
        print("WARNING: -AB geometry file not provided or not existing! EXIT")
        sys.exit()
    A, B, na, nb = read_xyzlabel_file(args.geoab)
    write_xyz_file(args.outfilename + "_label.xyz", B, A, label=True)

if args.function == "makeAB":
    if args.geoa is None or not os.path.isfile(args.geoa):
        print("WARNING: -A geometry file not provided or not existing! EXIT")
        sys.exit()
    if args.geob is None or not os.path.isfile(args.geob):
        print("WARNING: -B geometry file not provided or not existing! EXIT")
        sys.exit()
    if args.cell is None or not os.path.isfile(args.cell):
        print("WARNING: -cell file not provided or not existing! EXIT")
        sys.exit()
    A = read_xyz_file(args.geoa)
    B = read_xyz_file(args.geob)
    cell = read_cell_file(args.cell)
    cell.compute_all_info()
    A.compute_fract(cell)
    B.compute_fract(cell)
    for i in range(args.nout):
        print("makeAB: combination #%d " % i, end="")
        for j in range(args.maxfail):
            print(".", end="")
            B_rand = translate_rand(rotate_rand(B), cell)
            if not has_overlap(A, B_rand, cell, args.srad, args.printmindist):
                print(" not-overlapping configuration found!")
                write_xyz_file(
                    args.outfilename + "_" + str(i) + "_nolabel.xyz",
                    A,
                    B_rand,
                    label=False)
                write_xyz_file(
                    args.outfilename + "_" + str(i) + "_label.xyz",
                    A,
                    B_rand,
                    label=True)
                break
            if j == args.maxfail - 1:
                print("STOP")
                print("WARNING: -maxfail threshold reached. EXIT")
                print("TIP 1: increase -maxfail")
                print("TIP 2: lower -srad")
                print("TIP 3: check if your B molecule is too big \
                      for your A framework")
                sys.exit()

if args.function == "debug_kind":
    """ DEBUG utility to parse a .kind file with a geometry
    and print a new .kind file"""
    print("DEBUG: Reading %s geometry as A" % args.geoa)
    A = read_xyz_file(args.geoa)
    print("DEBUG: Reading %s kind file for A" % args.kindfile)
    read_kind_file(args.kindfilea, A, None)
    X = mol(0)
    X.compute_types()
    print("DEBUG: Printing %s kind file for A, after parsing" %
          args.outfilename + ".kind")
    write_kind_file(args.outfilename + ".kind", A, X, False, False, args.bs,
                    args.pot)

if args.function == "debug_split":
    """ DEBUG utility to read ab.wfn and AB_label.xyz and print a.wfn + b.wfn"""
    print("DEBUG: Reading %s geometry as AB" % args.geoab)
    A, B, na, nb = read_xyzlabel_file(args.geoab)
    print("DEBUG: Reading %s wavefunction as ab" % args.wfnab)
    ab = read_wfn_file(args.wfnab)
    print("DEBUG: Splitting the wfn")
    a, b = split_wfn(ab, A, B, args.qa, args.qb, args.multa, args.multb)
    print("DEBUG: Printing %s and %s" % (args.outfilename + "_Aa.fwfn",
                                         args.outfilename + "_Bb.fwfn"))
    write_fwfn_file(a, args.outfilename + "_Aa.fwfn")
    write_fwfn_file(b, args.outfilename + "_Bb.fwfn")
