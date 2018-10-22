#!/usr/bin/env python3

from atomic_data import atomic_symbols, atomic_names
from utilities import wfn, read_wfn_file, write_wfn_file, write_fwfn_file, combine_wfn
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
                      help="Binary .wfn function file for fragment A")

parser.add_argument("-o","--outfilename",
                      type=str,
                      dest="outfilename",
                      default="wfntools_output",
                      help="Name of the output file(s)")

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
    print('Converting wfns %s and %s into %s.wfn' %(args.wfna,args.wfna,args.outfilename))
    a=read_wfn_file(args.wfna)
    b=read_wfn_file(args.wfnb)
    ab=combine_wfn(a,b)
    write_wfn_file(ab,args.outfilename+'.wfn')
    if args.printformatted:
        write_fwfn_file(ab,args.outfilename+'.fwfn')
