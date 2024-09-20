from structure_clustering import Atom, Structure


def element_to_atomic_number(symbol):
    # fmt:off
    elements_to_atomic_numbers = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20,
    "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30,
    "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36, "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40,
    "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60,
    "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78, "Au": 79, "Hg": 80,
    "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90,
    "Pa": 91, "U": 92, "Np": 93, "Pu": 94, "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100,
    "Md": 101, "No": 102, "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109,
    "Ds": 110, "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118
    }
    # fmt:on
    return elements_to_atomic_numbers.get(symbol)


def import_multi_xyz(file_path):
    """
    Reads a multi-XYZ file and returns the data as a list of dictionaries.

    Each dictionary represents a frame and contains:
    - 'atom_count': the number of atoms in the frame
    - 'atoms': a list of tuples, each containing (element, x, y, z)

    Args:
    - file_path (str): Path to the multi-XYZ file

    Returns:
    - List of dictionaries, where each dictionary represents a frame
    """
    structures = []
    current_structure = None

    with open(file_path, "r") as file:
        lines = file.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()

        # check if line is a number indicating the number of atoms
        if line.isdigit():
            if current_structure:
                structures.append(current_structure)
            current_structure = Structure(i)

        # otherwise, process atom lines
        else:
            if current_structure:
                parts = line.split()
                if len(parts) == 4:
                    x, y, z = map(float, parts[1:])
                    current_structure.addAtom(
                        Atom(element_to_atomic_number(parts[0]), x, y, z)
                    )

        i += 1

    # append the last structure
    if current_structure:
        structures.append(current_structure)

    return structures
