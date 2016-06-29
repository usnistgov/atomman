import atomman as am
import numpy as np
import atomman.unitconvert as uc

try:
    import pymatgen as pmg
    has_pmg = True
except:
    has_pmg = False
    
def from_pymatgen_Struct(pmg_struct, scale=True):
    """Convert a pymatgen.Structure into an atomman.System and list of elements."""
    assert has_pmg, 'pymatgen not imported'

    box = am.Box(vects=pmg_struct.lattice.matrix)

    prop = pmg_struct.site_properties
    prop['pos'] = pmg_struct.cart_coords

    all_elements =  np.array([str(symbol) for symbol in pmg_struct.species])
    elements, atype = np.unique(all_elements, return_inverse=True)
    prop['atype'] = atype + 1

    atoms = am.Atoms(natoms=pmg_struct.num_sites, prop=prop)
    
    return am.System(box=box, atoms=atoms, scale=scale), elements

    
    
def to_pymatgen_Struct(system, elements):
    """Convert an atomman.System and list of elements into a pymatgen.Structure."""
    assert has_pmg, 'ase not imported'

    elements = np.asarray(elements)
    atype = system.atoms_prop(key='atype')
    
    latt = pmg.Lattice(system.box.vects)
    species = elements[atype-1]
    sites = system.atoms_prop(key='pos', scale=True)
    
    properties = {}
    for prop in system.atoms_prop():
        if prop != 'atype' and prop != 'pos':
            properties[prop] = system.atoms_prop(key=prop)

    return pmg.Structure(latt, species, sites, site_properties=properties)
