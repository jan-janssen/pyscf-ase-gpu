#!/usr/bin/env python3
# 
# Copied from https://github.com/heini-phys-chem/ASE_calculators/blob/master/calculators/pyscf_simple.py

import numpy as np
from ase.build import molecule
from ase.calculators.calculator import Calculator
from ase.atoms import Atoms

from pyscf import gto, scf, grad, mp, dft


convert_energy    =  27.2114  # hartree to eV
convert_forces    =  -27.2114 / 0.529177  # Bohr to Angstrom
convert_positions =  0.529177 # Bohr to Angstrom 


def ase_atoms_to_pyscf(ase_atoms):
  return [ [ase_atoms.get_chemical_symbols()[i], ase_atoms.get_positions()[i]] for i in range(len(ase_atoms.get_positions()))]


class PySCF(Calculator):
	name = 'PySCF'
	implemented_properties = ['energies', 'forces']

	def __init__(self, atoms, basis='def2-tzvpp', xc='LDA'):
		self.basis = basis
		self.xc = xc
		self.results = {}

	def get_potential_energy(self, atoms=None, force_consistent=False):
		mf = dft.RKS(gto.M(atom=ase_atoms_to_pyscf(atoms), basis=self.basis, charge=0, verbose=0))
		mf.xc = self.xc
		energy = mf.kernel()
		energy *= convert_energy

		self.results['energy'] = energy
		return energy

	def get_forces(self, atoms=None):
		mf = dft.RKS(gto.M(atom=ase_atoms_to_pyscf(atoms), basis=self.basis, charge=0, verbose=0)).run()
		mf.xc = self.xc
		gradient = mf.nuc_grad_method().kernel()
		forces = gradient * convert_forces

		self.results['forces'] = forces
		return forces


if __name__ == "__main__":
    h2o = molecule("H2O")
    h2o.set_calculator(PySCF(atoms=h2o, basis='def2-tzvpp', xc='LDA'))
    print(h2o.get_potential_energy())
    print(h2o.get_forces())
