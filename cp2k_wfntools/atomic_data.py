atomic_symbols = ("X", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
                  "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
                  "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
                  "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
                  "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
                  "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
                  "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
                  "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
                  "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
                  "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
                  "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
                  "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og")

atomic_names = ("DUMMY", "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron",
                "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon", "Sodium",
                "Magnesium", "Aluminium", "Silicon", "Phosphorus", "Sulfur",
                "Chlorine", "Argon", "Potassium", "Calcium", "Scandium",
                "Titanium", "Vanadium", "Chromium", "Manganese", "Iron",
                "Cobalt", "Nickel", "Copper", "Zinc", "Gallium", "Germanium",
                "Arsenic", "Selenium", "Bromine", "Krypton", "Rubidium",
                "Strontium", "Yttrium", "Zirconium", "Niobium", "Molybdenum",
                "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver",
                "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine",
                "Xenon", "Caesium", "Barium", "Lanthanum", "Cerium",
                "Praseodymium", "Neodymium", "Promethium", "Samarium",
                "Europium", "Gadolinium", "Terbium", "Dysprosium", "Holmium",
                "Erbium", "Thulium", "Ytterbium", "Lutetium", "Hafnium",
                "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium",
                "Platinum", "Gold", "Mercury", "Thallium", "Lead", "Bismuth",
                "Polonium", "Astatine", "Radon", "Francium", "Radium",
                "Actinium", "Thorium", "Protactinium", "Uranium", "Neptunium",
                "Plutonium", "Americium", "Curium", "Berkelium", "Californium",
                "Einsteinium", "Fermium", "Mendelevium", "Nobelium",
                "Lawrencium", "Rutherfordium", "Dubnium", "Seaborgium",
                "Bohrium", "Hassium", "Meitnerium", "Darmstadtium",
                "Roentgenium", "Copernicium", "Nihonium", "Flerovium",
                "Moscovium", "Livermorium", "Tennessine", "Oganesson")

# According to the default in CP2K MOLOPT basis set
atomic_valence_default_dict = {
    "H": 1,
    "He": 2,
    "Li": 3,
    "Be": 4,
    "B": 3,
    "C": 4,
    "N": 5,
    "O": 6,
    "F": 7,
    "Ne": 8,
    "Na": 9,
    "Mg": 2,  #default: 10, but problematic
    "Al": 3,
    "Si": 4,
    "P": 5,
    "S": 6,
    "Cl": 7,
    "Ar": 8,
    "K": 9,
    "Ca": 10,
    "Sc": 11,
    "Ti": 12,
    "V": 13,
    "Cr": 14,
    "Mn": 15,
    "Fe": 16,
    "Co": 17,
    "Ni": 18,
    "Cu": 11,
    "Zn": 12,
    "Ga": 13,
    "Ge": 4,
    "As": 5,
    "Se": 6,
    "Br": 7,
    "Kr": 8,
    "Rb": 9,
    "Sr": 10,
    "Y": 11,
    "Zr": 12,
    "Nb": 13,
    "Mo": 14,
    "Tc": 15,
    "Ru": 16,
    "Rh": 9,
    "Pd": 18,
    "Ag": 11,
    "Cd": 12,
    "In": 13,
    "Sn": 4,
    "Sb": 5,
    "Te": 6,
    "I": 7,
    "Xe": 8,
    "Cs": 9,
    "Ba": 10,
    "La": 11,
    "Ce": 12,
    "Pr": 13,
    "Nd": 14,
    "Pm": 15,
    "Sm": 16,
    "Eu": 17,
    "Gd": 18,
    "Tb": 19,
    "Dy": 20,
    "Ho": 21,
    "Er": 22,
    "Tm": 23,
    "Yb": 24,
    "Lu": 25,
    "Hf": 12,
    "Ta": 13,
    "W": 14,
    "Re": 15,
    "Os": 16,
    "Ir": 17,
    "Pt": 18,
    "Au": 11,
    "Hg": 12,
    "Tl": 13,
    "Pb": 4,
    "Bi": 5,
    "Po": 6,
    "At": 7,
    "Rn": 8,
}

# From UFF: L-J's sigma * 0.5 (doi: 10.1021/ja00051a040)
atomic_rad_UFF_dict = {
    "H": 1.286,
    "He": 1.052,
    "Li": 1.092,
    "Be": 1.223,
    "B": 1.819,
    "C": 1.715,
    "N": 1.630,
    "O": 1.559,
    "F": 1.498,
    "Ne": 1.445,
    "Na": 1.329,
    "Mg": 1.346,
    "Al": 2.004,
    "Si": 1.913,
    "P": 1.847,
    "S": 1.797,
    "Cl": 1.758,
    "Ar": 1.723,
    "K": 1.698,
    "Ca": 1.514,
    "Sc": 1.468,
    "Ti": 1.414,
    "V": 1.400,
    "Cr": 1.347,
    "Mn": 1.319,
    "Fe": 1.297,
    "Co": 1.279,
    "Ni": 1.262,
    "Cu": 1.557,
    "Zn": 1.231,
    "Ga": 1.952,
    "Ge": 1.907,
    "As": 1.884,
    "Se": 1.873,
    "Br": 1.866,
    "Kr": 1.845,
    "Rb": 1.833,
    "Sr": 1.622,
    "Y": 1.490,
    "Zr": 1.392,
    "Nb": 1.410,
    "Mo": 1.360,
    "Tc": 1.335,
    "Ru": 1.320,
    "Rh": 1.305,
    "Pd": 1.291,
    "Ag": 1.402,
    "Cd": 1.269,
    "In": 1.988,
    "Sn": 1.956,
    "Sb": 1.969,
    "Te": 1.991,
    "I": 2.005,
    "Xe": 1.962,
    "Cs": 2.012,
    "Ba": 1.649,
    "La": 1.569,
    "Ce": 1.584,
    "Pr": 1.606,
    "Nd": 1.592,
    "Pm": 1.580,
    "Sm": 1.568,
    "Eu": 1.556,
    "Gd": 1.500,
    "Tb": 1.537,
    "Dy": 1.527,
    "Ho": 1.519,
    "Er": 1.511,
    "Tm": 1.503,
    "Yb": 1.494,
    "Lu": 1.621,
    "Hf": 1.399,
    "Ta": 1.412,
    "W": 1.367,
    "Re": 1.316,
    "Os": 1.390,
    "Ir": 1.265,
    "Pt": 1.227,
    "Au": 1.467,
    "Hg": 1.205,
    "Tl": 1.936,
    "Pb": 1.914,
    "Bi": 1.947,
    "Po": 2.098,
    "At": 2.116,
    "Rn": 2.123,
    "Fr": 2.183,
    "Ra": 1.638,
    "Ac": 1.549,
    "Th": 1.513,
    "Pa": 1.525,
    "U": 1.512,
    "Np": 1.525,
    "Pu": 1.525,
    "Am": 1.506,
    "Cm": 1.482,
    "Bk": 1.487,
    "Cf": 1.476,
    "Es": 1.470,
    "Fm": 1.464,
    "Md": 1.458,
    "No": 1.447,
    "Lr": 1.441,
}