import hartree_fock as hf

orbitals = [(1, 0, 0), (2, 0, 0), (2, 1, 0), (3, 0, 0),
            (3, 1, 0), (3, 2, 0), (4, 0, 0), (4, 1, 0),]

hf.create_atom(36, 0)
hf.add_orbitals(orbitals)
hf.set_grid()
hf.hartree_fock()
