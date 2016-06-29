import atomman as am

try:
    import ase
    has_ase = True
except:
    has_ase = False
    
def from_ase_Atoms(ase_atoms):
    """Convert an ase.Atoms into an atomman.System and list of elements."""
    assert has_ase, 'ase not imported'
    box = am.Box(vects = ase_atoms.get_cell())
    atoms = am.Atoms(natoms=len(ase_atoms))
    atoms.prop(key = 'pos', value=ase_atoms.get_positions())

    all_elements = np.array(ase_atoms.get_chemical_symbols())
    elements, atype = np.unique(all_elements, return_inverse)
    atype += 1
    
    atoms.prop(key='atype', value=atype)
    
    return am.System(atoms=atoms, box=box), elements

def to_ase_Atoms(system, elements):
    """Convert an atomman.System and list of elements into an ase.Atoms."""
    assert has_ase, 'ase not imported'
    
    positions = system.atoms_prop(key='pos')
    cell = system.box.vects
    atypes = system.atoms_prop(key='atype')
    symbols = elements[atype-1]
    pbc = system.pbc
    
    return ase.Atoms(symbols=symbols, positions=positions, pbc=pbc, cell=cell)
    
