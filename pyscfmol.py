from pyscf import gto
import numpy as np
from ase.build import molecule


def ase_atoms_to_pyscf(ase_atoms):
    '''
    Convert ASE atoms to PySCF atom.
    Note: ASE atoms always use A.
    '''
    return [[atom.symbol, atom.position] for atom in ase_atoms]


if __name__ == "__main__":
    atoms = molecule("H2O")
    calc_molcell = gto.M()
    calc_molcell.atom = ase_atoms_to_pyscf(atoms)
    calc_molcell.a = np.asarray(atoms.cell)
    if hasattr(atoms, "info"):
        if isinstance(atoms.info, dict):
            charge = atoms.info.get("charge", None)
            uhf = atoms.info.get("uhf", None)
    if charge is None:  # Read from initial charges
        charge = np.sum(atoms.get_initial_charges())
    if uhf is None:  # Read from initial magnetic moments
        uhf = np.sum(atoms.get_initial_magnetic_moments())
    calc_molcell.charge = charge
    calc_molcell.spin = uhf
    print(calc_molcell)
