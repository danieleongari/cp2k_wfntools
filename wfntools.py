#!/usr/bin/env python3

import sys
import os
from atomic_data import atomic_symbols, atomic_names
from utilities import wfn, read_wfn_file, write_wfn_file, write_fwfn_file, combine_wfn
from utilities import mol, read_xyz_file, write_xyz_file, write_kind_file
import argparse
from argparse import RawTextHelpFormatter #needed to go next line in the help text

parser = argparse.ArgumentParser(description="Program to parse and combine CP2K's .wfn wavefunction files", formatter_class=RawTextHelpFormatter)

parser.add_argument("function",
                      choices=['parse','combine','cp','makeAB','labelAB'],
                      help='Function to execute')

parser.add_argument("-a","--wfna",
                      type=str,
                      dest="wfna",
                      help="Binary .wfn function file for fragment A")

parser.add_argument("-b","--wfnb",
                      type=str,
                      dest="wfnb",
                      default=None,
                      help="Binary .wfn function file for fragment A")

parser.add_argument("-ab","--wfnab",
                      type=str,
                      dest="wfnab",
                      default=None,
                      help="Binary .wfn function file for both fragments together")

parser.add_argument("-A","--geoa",
                      type=str,
                      dest="geoa",
                      help="Geometry of fragment A")

parser.add_argument("-B","--geob",
                      type=str,
                      dest="geob",
                      default=None,
                      help="Geometry of fragment B")

parser.add_argument("-AB","--geoab",
                      type=str,
                      dest="geoab",
                      default=None,
                      help="Geometry of both fragments together (NB: should contain labels!)")

parser.add_argument("-o","--outfilename",
                      type=str,
                      dest="outfilename",
                      default="wfntools_output",
                      help="Name of the output file(s)")

parser.add_argument("-bs","--basisiset",
                      type=str,
                      dest="bs",
                      default="DZVP-MOLOPT-SR-GTH",
                      help="CP2K basis set choice for the .kind file")

parser.add_argument("-pot","--potential",
                      type=str,
                      dest="pot",
                      default="GTH-PBE",
                      help="CP2K pseudopotential choice for the .kind file")

parser.add_argument("-sbs","--singlebasisset",
                      action="store_true",
                      dest="sbs",
                      default=False,
                      help="In the 'combine' calculation, print also single basis set files for A and B")

parser.add_argument("-na","--numberatomsa",
                      type=int,
                      dest="na",
                      default=None,
                      help="Number of atoms of fragment A")

parser.add_argument("-nb","--numberatomsb",
                      type=int,
                      dest="nb",
                      default=None,
                      help="Number of atoms of fragment B")

parser.add_argument("-swapab","--swapAwithB",
                      action="store_true",
                      dest="swapab",
                      default=False,
                      help="In the .xyz file read B first and A later")

parser.add_argument("-f","--printformatted",
                      action="store_true",
                      dest="printformatted",
                      default=False,
                      help="When a .wfn file is printed, it prints also a formatted .fwfn file")

args = parser.parse_args()

if args.function=='parse':
    print('Converting %s binary file in %s.fwfn formatted file.' %(args.wfna,args.outfilename))
    a=read_wfn_file(args.wfna)
    write_fwfn_file(a,args.outfilename+'.fwfn')

if args.function=='combine':
    #if wfargs.na and wfnb are provided, combine the two wfn and print fwfn if asked
    if (args.wfargs.na == None or args.wfnb == None):
        print('WARNING: wavefunctions for A or B were not provided!')
    if (args.wfargs.na != None and args.wfnb != None):
        print('Converting wfns %s and %s into %s.wfn' %(args.wfna,args.wfnb,args.outfilename))
        a=read_wfn_file(args.wfna)
        b=read_wfn_file(args.wfnb)
        ab=combine_wfn(a,b)
        write_wfn_file(ab,args.outfilename+'.wfn')
        if args.printformatted:
            write_fwfn_file(ab,args.outfilename+'.fwfn')
    #if geoa and geob are provided, combine the two geometries and make the kind file
    if (args.geoa == None or args.geob == None):
        print('WARNING: geometries for A or B were not provided!')
    if (args.geoa != None and args.geob != None):
        print('Converting geometries %s and %s into %s_label.xyz' %(args.geoa,args.geob,args.outfilename))
        A=read_xyz_file(args.geoa)
        B=read_xyz_file(args.geob)
        write_xyz_file(args.outfilename+'_nolabel.xyz',A,B,label=False)
        write_xyz_file(args.outfilename+'_label.xyz',A,B,label=True)
        write_kind_file(args.outfilename+'.kind',A,B,False,False,args.bs,args.pot)

if args.function=='labelAB':
    # Read the number of atoms an check if everything makes sense: na+nb=nab
    inpfile = open(args.geoab,'r')
    data = inpfile.readline()
    nab = int(data[0])
    if args.na== None and args.nb == None:
        print("WARNING: number of atoms for A (or B) not provided! EXIT")
        sys.exit()
    if args.na!= None and args.nb != None:
        if nab != (args.na + args.nb):
             print("WARNING: number of atoms for A and B both specified, but they don't match with nAB! EXIT")
             sys.exit()
    if args.na== None:
        args.na= nab - args.nb
    else:
        args.nb = nab - args.na
    if (args.na <= 0) or (args.nb <= 0):
        print("WARNING: something is wrong with the number of atoms! EXIT")
        sys.exit()

    # Print two temporary files for A and B geometries
    outfileA = open('tmp_A.xyz','w')
    outfileB = open('tmp_B.xyz','w')
    print(args.na, file=outfileA)
    print(args.nb, file=outfileB)
    print('Fragment A', file=outfileA)
    print('Fragment B', file=outfileB)
    junk = inpfile.readline()
    if not args.swapab:
        for i in range(args.na):
            line = inpfile.readline()
            print(line, end='', file=outfileA)
        for i in range(args.nb):
            line = inpfile.readline()
            print(line, end='', file=outfileB)
    else:
        for i in range(args.nb):
            line = inpfile.readline()
            print(line, end='', file=outfileB)
        for i in range(args.na):
            line = inpfile.readline()
            print(line, end='', file=outfileA)
    outfileA.close()
    outfileB.close()
    A=read_xyz_file('tmp_A.xyz')
    B=read_xyz_file('tmp_B.xyz')
    os.remove('tmp_A.xyz')
    os.remove('tmp_B.xyz')
    write_xyz_file(args.outfilename+'_label.xyz',A,B,label=True)
