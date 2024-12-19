# Benchmark

## CPU 
```
python  pyscfasecpu.py
```

Output:
```
-2047.482028073281
[[-2.16207850e-14 -4.47763055e-13  1.66039055e+00]
 [-7.08441892e-15  1.20453596e+00 -8.30458046e-01]
 [ 6.20610842e-14 -1.20453596e+00 -8.30458046e-01]]
```

## GPU
```
python pyscfasegpu.py
```

Output: 
```
-2047.4820280735028
[[ 2.07487556e-13 -3.37737131e-12  1.66038073e+00]
 [-1.09305783e-13  1.20452916e+00 -8.30453138e-01]
 [-9.50707678e-14 -1.20452916e+00 -8.30453138e-01]]
```

## References:
* [pyscf/gpu4pyscf](https://github.com/pyscf/gpu4pyscf)
* [pyscf/pyscf_ase](https://github.com/pyscf/pyscf/blob/master/pyscf/pbc/tools/pyscf_ase.py)
* [heini-phys-chem/ASE_calculators](https://github.com/heini-phys-chem/ASE_calculators)
* [chelleorc/pyscf_ase_github](https://github.com/chelleorc/pyscf_ase_github)
