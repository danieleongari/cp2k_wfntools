#!/usr/bin/env python

from __future__ import print_function  #python3 like print()
from __future__ import absolute_import
import os

import sys
import subprocess
import shutil #to removed not-empty directories
import filecmp #to compare two files

shutil.rmtree("./output")
os.makedirs("./output")

######################################################################## Test 01
print("Test 01: parse wfn")
cmd = ["../wfntools.py",
       "parse",
       "-a",
       "./input/water2_ANO+MOLOPT.wfn",
       "-o",
       "./output/test01",
       ]
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
output, error = process.communicate()

expected_files=["test01.fwfn"
                ]
for f in expected_files:
    if os.path.exists("./output/"+f) \
       and filecmp.cmp("./output/"+f, "./reference/"+f):
        print("\tOutput file %s: OK" %f )
    else:
        print("\tOutput file %s: MISSING or DIFFERENT" %f )

######################################################################## Test 02
print("Test 02: combine wfn")
cmd = ["../wfntools.py",
       "combine",
       "-a",
       "./input/water_DZ-ANO-ALL.wfn",
       "-b",
       "./input/water_DZVP-MOLOPT-SR.wfn",
       "-A",
       "./input/water_A.xyz",
       "-B",
       "./input/water_B.xyz",
       "-kindA",
       "./input/water_A.kind",
       "-kindB",
       "./input/water_B.kind",
       "-bs",
       "read",
       "-pot",
       "read",
       "-o",
       "./output/test02",
       "-pf",
       ]
process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
output, error = process.communicate()

expected_files=["test02.fwfn",
                "test02.kind", #can have a different order
                "test02_label.xyz",
                "test02_nolabel.xyz",
                "test02.wfn",
                ]
for f in expected_files:
    if os.path.exists("./output/"+f) \
       and filecmp.cmp("./output/"+f, "./reference/"+f):
        print("\tOutput file %s: OK" %f )
    else:
        print("\tOutput file %s: MISSING or DIFFERENT" %f )
