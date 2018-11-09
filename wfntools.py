#!/usr/bin/env python3

from atomic_data import atomic_symbols, atomic_names
from utilities import wfn, read_wfn_file, write_wfn_file, write_fwfn_file, combine_wfn
from utilities import mol, read_xyz_file, write_xyz_file, write_kind_file
import argparse
from argparse import RawTextHelpFormatter #needed to go next line in the help text

parser = argparse.ArgumentParser(description="Program to parse and combine CP2K's .wfn wavefunction files", formatter_class=RawTextHelpFormatter)

parser.add_argument("function",
                      choices=['parse','combine','cp'],
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

parser.add_argument("-A","--geoa",
                      type=str,
                      dest="geoa",
                      help="Geometry of fragment A")

parser.add_argument("-B","--geob",
                      type=str,
                      dest="geob",
                      default=None,
                      help="Geometry of fragment B")

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

parser.add_argument("-f","--printformatted",
                      action="store_true",
                      dest="printformatted",
                      default="False",
                      help="When a .wfn file is printed, it prints also a formatted .fwfn file")

parser.add_argument("-makeAB",
                      action="store_true",
                      dest="silent",
                      default=False,
                      help="Combine the geometries, inserting B in A in a random but non-overlapping position")

args = parser.parse_args()

if args.function=='parse':
    print('Converting %s binary file in %s.fwfn formatted file.' %(args.wfna,args.outfilename))
    a=read_wfn_file(args.wfna)
    write_fwfn_file(a,args.outfilename+'.fwfn')

if args.function=='combine':
    #if wfna and wfnb are provided, combine the two wfn and print fwfn if asked
    if (args.wfna == None or args.wfnb == None):
        print('WARNING: wavefunctions for A or B were not provided!')
    if (args.wfna != None and args.wfnb != None):
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
