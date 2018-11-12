# wfnTools
### Tools to parse and combine CP2K's .wfn wavefunction files

Consider:
* a sistem with two fragments, we want to compute the interaction
* the wave function / basis set of the fragmentA, a
* the wave function / basis set of the fragmentB, b
* the geometry of the fragmentA, A
* the geometry of the fragmentB, B
* by convention, A is the framework (or surface, bigger molecule) and B is the adsorbate
* Aa and Aab, means the wfn of fragmentA with the basis set of A (a), or both A and B (ab)
* the non-counterpoise-corrected interaction energy is: E=E(ABab)-E(Aa)-E(Bb)
* the counterpoise-corrected interaction energy is: E=E(ABab)-E(Aab)-E(Bab)
* E=E(Aab)-E(Aa) is the difference in the energy due to a better minimization of the wave function of the fragment A when the basis set of B is added. This value should be negative.
* Different atom types can be specified for the same element: use the convention "type = elementX" with X being a scalar. For example, in an antiferromagnetic calculation one can use Cu1 and Cu2 to specify atoms with opposite spin. It is possible to use any number (e.g., Cu10, Cu3333), but NOT letters or symbols (e.g., Cu_1, CuUP, CuSpin1, ..., are not read correctly)

### Examples:
#### 1) Print the formatted wavefunction:
```
wfntools parse -a fragment.wfn -o outputname
```
Prints the file `outputname.fwfn`

#### 2) Combine two wavefunctions:
```
wfntools combine -a fragmentA.wfn -b fragmentB.wfn -o outputname
```
Prints the combined `outputname.wfn` of the two fragments

#### 3) Combine two wavefunctions and geometries:
```
wfntools combine -a fragmentA.wfn -b fragmentB.wfn -A fragmentA.xyz -B fragmentB.xyz -o outputname
```
* Prints the combined `outputname.wfn` of the two fragments.
* Prints the combined `outputname_nolabel.xyz` of the two fragments, for visualization.
* Prints the combined `outputname_label.xyz` of the two fragments, to be used for CP2K with "A" and "B" labels for the atoms.
* Prints the combined `outputname.kind` h the CP2K's &KIND section. Use `-bs` and `-pot` options to specify the basis set and the pseudopotential choice (TODO: read it from a CP2K input or output file)
TODO: use the `-check` keyword to check if the number of atoms are consistent in the geometry and wfn files.

#### 4) TODO: Prepare the inputs for a Counterpoise Corrected calculation (recycling the wavefunctions)
```
wfntools cp -ab bothfragments.wfn -AB bothfragments_label.xyz -o outputname
```
* Prints outputname_A.xyz (no label), outputname_Aab.wfn and outputname_Aab.kind
* Prints outputname_B.xyz (no label), outputname_Bab.wfn and outputname_Bab.kind
and, if `-sbs` option is used:
* Prints outputname_Aa.wfn and outputname_Aa.kind
* Prints outputname_Bb.wfn and outputname_Bb.kind
**NB: if the ab.wfn was not computed with labelled geometry, the method does not work, because the atomic &KIND are mixed!**
Advanced options:
* `-qA` and `-qB` to specify the charge of the fragments (default: 0)
* `multA` and `multB` to specify the multiplicity of the fragments (default: 1)

### Utilities:

#### 5) TODO: Insert the fragment B in a pore of fragment A:
```
wfntools makeAB -A fragment1.xyz -B fragment2.xyz -cell unitcell1.cell -o outputname -nout 10
```
This is an usefull utility to generate a number of `-nout` starting configurations of combined geometries, possibly to later adjust in Avogadro. Fragment B is inserted in a random not-overlapping position, conisdering the atoms as rigid spheres with a defined radius (from UFF) possibly scaled by a factor `-scalerad`.
Prints `outputname_AB_1.xyz` `outputname_B_1.xyz` `outputname_AB_2.xyz` `outputname_B_2.xyz` ... `outputname_AB_nout.xyz` `outputname_B_nout.xyz` ...

#### 6) Label A and B in the AB geometry file:
```
wfntools labelAB -AB both_fragments.xyz -nA 156 -o outputname
```
or
```
wfntools labelAB -AB both_fragments.xyz -nB 3 -o outputname
```
 Print a new `outputname_label.xyz` file where A and B are labelled according to the number of `-nA` atoms in A (or specifying `-nB`). Use the `-swapAB` option to read first B and then A in the .xyz file.

#### 7) Swap A and B labels in the AB geometry file:
```
wfntools swaplAB -AB both_fragments.xyz -o outputname
```
Print a new `outputname_label.xyz` file where A and B are labels are swapped.
