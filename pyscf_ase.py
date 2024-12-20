#!/usr/bin/env python

# Copyright 2014-2020 The PySCF Developers. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Garnet Chan <gkc1000@gmail.com>
#         Xing Zhang <zhangxing.nju@gmail.com>
#

'''
ASE package interface
'''

import scipy
import numpy as np
from ase.calculators.calculator import Calculator
import ase.dft.kpoints

BOHR_TO_ANGSTROM = (
    scipy.constants.physical_constants["Bohr radius"][0] / scipy.constants.angstrom
)
HARTREE_TO_EV = scipy.constants.physical_constants["Hartree energy in eV"][0]

convert_energy = HARTREE_TO_EV
convert_forces = -HARTREE_TO_EV / BOHR_TO_ANGSTROM


def pyscf_to_ase_atoms(cell):
    '''
    Convert PySCF Cell/Mole object to ASE Atoms object
    '''
    from ase import Atoms
    from pyscf.lib import param
    from pyscf.pbc import gto

    symbols = cell.elements
    positions = cell.atom_coords() * param.BOHR
    if isinstance(cell, gto.Cell):
        a = cell.lattice_vectors() * param.BOHR
        return Atoms(symbols, positions, cell=a, pbc=True)
    else:
        return Atoms(symbols, positions, pbc=False)

def ase_atoms_to_pyscf(ase_atoms):
    '''Convert ASE atoms to PySCF atom.

    Note: ASE atoms always use A.
    '''
    return [[atom.symbol, atom.position] for atom in ase_atoms]
atoms_from_ase = ase_atoms_to_pyscf

class PySCF(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, restart=None, label='PySCF', atoms=None, scratch=None, **kwargs):
        """Construct PySCF-calculator object.

        Parameters
        ==========
        label: str
            Prefix to use for filenames (label.in, label.txt, ...).
            Default is 'PySCF'.

        mfclass: PySCF mean-field class
        molcell: PySCF :Mole: or :Cell:
        """
        Calculator.__init__(self, restart=restart, label=label, atoms=atoms, scratch=scratch, **kwargs)

        # TODO
        # This explicitly refers to "cell". How to refer
        # to both cell and mol together?

        self.mf=None
        self.initialize(**kwargs)

    def initialize(self, molcell, mf_class, mf_dict):
        if not molcell.unit.startswith(('A','a')):
            raise RuntimeError("PySCF unit must be A to work with ASE")

        self.molcell=molcell
        self.mf_class=mf_class
        self.mf_dict=mf_dict

    def set(self, **kwargs):
        changed_parameters = Calculator.set(self, **kwargs)
        if changed_parameters:
            self.reset()

    def calculate(self, atoms=None, properties=['energy', 'forces'],
                  system_changes=['positions', 'numbers', 'cell',
                                  'pbc', 'charges','magmoms']):

        Calculator.calculate(self, atoms)

        calc_molcell = self.molcell.copy()
        calc_molcell.atom = ase_atoms_to_pyscf(atoms)
        calc_molcell.a = np.asarray(atoms.cell)
        # Begin - Set charges and spins
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
        calc_molcell.verbose = 0
        # End - Set charges and spins
        calc_molcell.build(None,None)
        self.mf = self.mf_class(calc_molcell)
        for key in self.mf_dict:
            self.mf.__dict__[key] = self.mf_dict[key]
        # Begin - Calculate Forces and Energy
        self.results['energy'] = self.mf.scf() * convert_energy
        self.results['forces'] = self.mf.nuc_grad_method().kernel() * convert_forces
        # End - Calculate Forces and Energy
        # self.results['energy']=self.mf.scf()
        self.results['mf']=self.mf

def make_kpts(cell, nks):
    raise DeprecationWarning('Use cell.make_kpts(nks) instead.')
