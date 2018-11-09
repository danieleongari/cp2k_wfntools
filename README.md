# wfnTools
### Tools to parse and combine CP2K's .wfn wavefunction files

Consider:
* the wavefunction of the fragment 1, a
* the wavefunction of the fragment 2, b 
* the geometry of the fragment 1, A
* the geometry of the fragment 2, B
 
### Examples:
#### 1) Print the formatted wavefunction:
```
wfntools parse -a fragment.wfn -o outputname
```
Prints the file `outputname.fwfn`

#### 2) Combine two wavefunctions:
```
wfntools combine -a fragment1.wfn -b fragment2.wfn -o outputname
```
Prints the combined `outputname.wfn` of the two fragments

#### 3) Combine two wavefunctions and geometries:
```
wfntools combine -a fragment1.wfn -b fragment2.wfn -A fragment1.xyz -B fragment2.xyz -o outputname
```
* Prints the combined `outputname.wfn` of the two fragments.
* Prints the combined `outputname_nolabel.xyz` of the two fragments, for visualization.
* Prints the combined `outputname_label.xyz` of the two fragments, to be used for CP2K with "_A" and "_B" labels for the atoms.
* Prints the combined `outputname.kind` h the CP2K's &KIND section. Use `-bs` and `-pot` options to specify the basis set and the pseudopotential choice (TODO: read it from a CP2K file


#### 4) TODO: Insert the fragment B in a pore of fragment A
```
wfntools makeAB -A fragment1.xyz -B fragment2.xyz -cell unitcell1.cell -o outputname -nout 10
```
This is an usefull utility to generate a number of `-nout` starting configurations of combined geometries, possibly to later adjust in Avogadro. Fragment B is inserted in a random not-overlapping position, conisdering the atoms as rigid spheres with a defined radius (from UFF) possibly scaled by a factor `-scalerad`.
Prints `outputname_AB_1.xyz` `outputname_B_1.xyz` `outputname_AB_2.xyz` `outputname_B_2.xyz` ... `outputname_AB_nout.xyz` `outputname_B_nout.xyz` ... 

#### 5) TODO: Prepare the inputs for a Counterpoise Corrected calculation (recycling the wavefunctions)
```
wfntools cp -AB bothfragments_label.xyz -o outputname
```

