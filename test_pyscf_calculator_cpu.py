from pyscf_ase import PySCF
from ase.build import molecule
from pyscf import gto, scf, grad, mp, dft


if __name__ == "__main__":
    h2o = molecule("H2O")
    h2o.set_calculator(PySCF(atoms=h2o, molcell=gto.M(basis='def2-tzvpp'), mf_class=dft.RKS, mf_dict={"xc": "LDA"}))
    print(h2o.get_potential_energy())
    print(h2o.get_forces())
