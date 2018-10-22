import scipy.io
import numpy as np

class wfn:
    """ Wave function object, as stored in a .wfn file """
    def __init__(self,
                 natom = None,
                 nspin = None,
                 nao_tot = None,
                 nset_max = None,
                 nshell_max = None):
        #defined at the initialization
        self.natom = natom              # number of atoms
        self.nspin = nspin              # number of spins (i.e., 1 or 2)
        self.nao_tot = nao_tot          # total number of atomic orbitals (=sum(self.nao))
        self.nset_max = nset_max        # max number of sets in the basis set (!!!! always 1???)
        self.nshell_max = nshell_max    # mak number of shells in each set
        #spin independent
        self.nset = [0 for i in range(self.natom)]                                  # number of sets for the bs of each atom
        self.nshell = [0 for i in range(self.natom)]                                # number of shells for the bs of each atom
        self.nao = [[0 for i in range(self.nshell_max)] for j in range(self.natom)] # number of atomic orbitals for each shell of each atom
        #spin dependent
        self.nmo = [0,0]    # number of molecular orbitals
        self.nocc = [0,0]   # number of occupied molecular orbitals
        self.nvirt = [0,0]  # number of virtual molecular orbitals
        self.nel = [0,0]    # number of electrons
    def initialize_lists(self,ispin):
        if ispin==0:
            self.eigen = [[0 for i in range(self.nmo[ispin])]]                                  # eigenvalues for each molecular orbital
            self.occup = [[0 for i in range(self.nmo[ispin])]]                                  # occupancies for each molecular orbital
            self.coeff = [[[0 for j in range(self.nao_tot)] for i in range(self.nmo[ispin])]]   #coefficient list for each moleculat orbital
        elif ispin==1:
            self.eigen.append([0 for i in range(self.nmo[ispin])])
            self.occup.append([0 for i in range(self.nmo[ispin])])
            self.coeff.append([[0 for j in range(self.nao_tot)] for i in range(self.nmo[ispin])])


def read_wfn_file(wfn_file):
    """ Parser for the .wfn file, that returns the wfn object """

    inpfile = scipy.io.FortranFile(wfn_file, 'r')

    natom, nspin, nao_tot, nset_max, nshell_max = inpfile.read_ints()
    w=wfn(natom,nspin,nao_tot,nset_max,nshell_max)

    nset = inpfile.read_ints()
    for i, val in enumerate(nset):
        w.nset[i]=val

    nshell = inpfile.read_ints()
    for i, val in enumerate(nshell):
        w.nshell[i]=val

    nao = inpfile.read_ints()
    for i, val in enumerate(nao):
        iatom, ishell = divmod(i,w.nshell_max) # this splits the atomic orbitals in one list for each atom
        w.nao[iatom][ishell]=val

    for ispin in range(w.nspin):
        nmo, nocc, nvirt, nel = inpfile.read_ints()
        w.nmo[ispin]=nmo
        w.nocc[ispin]=nocc
        w.nvirt[ispin]=nvirt
        w.nel[ispin]=nel
        w.initialize_lists(ispin)

        eigen_and_occup = inpfile.read_reals() #this line contains both
        for i in range(w.nmo[ispin]):
            w.eigen[ispin][i]=eigen_and_occup[i]
        for i in range(w.nmo[ispin]):
            w.occup[ispin][i]=eigen_and_occup[i+w.nmo[ispin]]

        for imo in range(w.nmo[ispin]):
            coeff = inpfile.read_reals()
            for i, val in enumerate(coeff):
                w.coeff[ispin][imo][i]=val
    inpfile.close()
    return w

def write_wfn_file(wfn_object,wfn_file):
    """ Print a wfn object as a binary fortran file """
    outfile = scipy.io.FortranFile(wfn_file, 'w')
    w=wfn_object
    outfile.write_record(np.array([w.natom, w.nspin, w.nao_tot, w.nset_max, w.nshell_max], dtype=np.int32))
    outfile.write_record(np.array(w.nset, dtype=np.int32))
    outfile.write_record(np.array(w.nshell, dtype=np.int32))
    outfile.write_record(np.array(np.concatenate(w.nao), dtype=np.int32))
    for ispin in range(w.nspin):
        outfile.write_record(np.array([w.nmo[ispin], w.nocc[ispin], w.nvirt[ispin], w.nel[ispin]], dtype=np.int32))
        outfile.write_record(np.array(np.concatenate([w.eigen[ispin],w.occup[ispin]]), dtype=np.dtype('f8')))
        for imo in range(w.nmo[ispin]):
            outfile.write_record(np.array(w.coeff[ispin][imo], dtype=np.dtype('f8')))
    outfile.close()
    return

def write_fwfn_file(wfn_object,fwfn_file):
    """ Print a wfn object as a formatted text file """
    outfile = open(fwfn_file, 'w')
    w=wfn_object
    print("Number of atoms: %d" %w.natom, file=outfile)
    print("Number of spins: %d" %w.nspin, file=outfile)
    print("Total number of atomic orbitals (AOs): %d" %w.nao_tot, file=outfile)
    print("Maximum number of sets: %d" %w.nset_max, file=outfile)
    print("Maximum number of shells: %d" %w.nshell_max, file=outfile)
    print("Number of sets for each atom: [", end='', file=outfile)
    for val in w.nset:
        print(" %d" %val, end='',file=outfile)
    print(" ]",file=outfile)
    print("Number of shells for each atom: [", end='', file=outfile)
    for val in w.nshell:
        print(" %d" %val, end='',file=outfile)
    print(" ]",file=outfile)
    print("Number of atomic orbitals for each shell for each atom:", file=outfile)
    for i, list in enumerate(w.nao):
        print("atom %d:\t[" %i, end='', file=outfile)
        for j in list:
            print(" %d" %j, end='',file=outfile)
        print(" ]",file=outfile)

    for ispin in range(w.nspin):
        print("", file=outfile)
        print("Information for spin %d:" %ispin, file=outfile)
        print("Number of molecular orbitals (MOs): %d" %w.nmo[ispin], file=outfile)
        print("Number of occupied MOs: %d" %w.nocc[ispin], file=outfile)
        print("Number of virtual MOs: %d" %w.nvirt[ispin], file=outfile)
        print("Number of electrons: %d" %w.nel[ispin], file=outfile)
        print("Eigenvalue for each MO: [", end='', file=outfile)
        for val in w.eigen[ispin]:
            print(" %.3f" %val, end='', file=outfile)
        print(" ]",file=outfile)
        print("Electron occupancy for each MO: [", end='', file=outfile)
        for val in w.occup[ispin]:
            print(" %.3f" %val, end='', file=outfile)
        print(" ]",file=outfile)
        print("Coefficients for the MO/AO matrix:", file=outfile)
        for i, list in enumerate(w.coeff[ispin]):
            print("MO %d:\t[" %i, end='', file=outfile)
            for j in list:
                print(" %+.3f" %j, end='', file=outfile)
            print(" ]",file=outfile)
    outfile.close()
    return


def combine_wfn(a,b):
    """ Combine wfns a and b into a single ab wfn object """

    natom=a.natom+b.natom
    nspin=max(a.nspin,b.nspin)
    nao_tot=a.nao_tot+b.nao_tot
    nset_max=max(a.nset_max,b.nset_max)
    nshell_max=max(a.nshell_max,b.nshell_max)
    w=wfn(natom,nspin,nao_tot,nset_max,nshell_max)

    nset = a.nset+b.nset
    for i, val in enumerate(nset):
        w.nset[i]=val

    nshell = a.nshell+b.nshell
    for i, val in enumerate(nshell):
        w.nshell[i]=val

    #BC: now it can be that w.nshell_max>a.shell_max and I need to leave some zeros
    for iatom in range(w.natom):
        if iatom<a.natom:
            for ishell in range(a.nshell_max):
                w.nao[iatom][ishell]=a.nao[iatom][ishell]
        else:
            for ishell in range(b.nshell_max):
                w.nao[iatom][ishell]=b.nao[iatom-a.natom][ishell]

    if a.nspin==1 and w.nspin==2: a=makeopenshell(a)
    if b.nspin==1 and w.nspin==2: b=makeopenshell(b)

    for ispin in range(w.nspin):
        w.nmo[ispin]=a.nmo[ispin]+b.nmo[ispin]
        w.nocc[ispin]=a.nocc[ispin]+b.nocc[ispin]
        w.nvirt[ispin]=a.nvirt[ispin]+b.nvirt[ispin] #not sure this is right!
        w.nel[ispin]=a.nel[ispin]+b.nel[ispin]
        w.initialize_lists(ispin)

        w.eigen[ispin]=a.eigen[ispin]+b.eigen[ispin]
        w.occup[ispin]=a.occup[ispin]+b.occup[ispin]

        #BC: you need to leave zeros for mixed AO/MO coeffs
        for imo in range(w.nmo[ispin]):
            if imo<a.nmo[ispin]:
                for iao in range(a.nao_tot):
                    w.coeff[ispin][imo][iao]=a.coeff[ispin][imo][iao]
            else:
                for iao in range(b.nao_tot):
                    w.coeff[ispin][imo][a.nao_tot+iao]=b.coeff[ispin][imo-a.nmo[ispin]][iao]
    return w

def makeopenshell(v):
    """ Makes a restricted (e.g., nspin=1) wfn open shell (e.g., nspin=2) """
    if v.nspin==2:
        return v
    elif v.nspin==1:
        w=wfn(v.natom,2,v.nao_tot,v.nset_max,v.nshell_max)
        w.nset=v.nset
        w.nshell=v.nshell
        w.nao=v.nao
        for ispin in range(2):
            w.nmo[ispin]=v.nmo[0]
            w.nocc[ispin]=v.nocc[0]
            w.nvirt[ispin]=v.nvirt[0]
            w.nel[ispin]=v.nel[0]/2
            w.initialize_lists(ispin)
            w.eigen[ispin]=[x/2. for x in v.eigen[0]]
            w.occup[ispin]=[x/2. for x in v.occup[0]]
            w.coeff[ispin]=v.coeff[0]
        return w
