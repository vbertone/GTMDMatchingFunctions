# GTMD matching functions at one loop

This repository accompanies the paper in which the one-loop GTMD
matching functions on GPDs are presented. The sample codes in the
folder `code` can be used to reproduce the plots in the paper.

In order to run the  codes, a few dependencies need to be
installed prior to compilation:

- [PARTONSv3](https://partons.cea.fr/partons/doc/html/index.html),
- [APFEL++](https://github.com/vbertone/apfelxx),
- [NangaParbat](https://github.com/MapCollaboration/NangaParbat).

## Compilation

In order to compile the code, you need to edit the `Makefile` in the
`code` folder to adjust the compiler (`CXX`) and set the correct path
to the PARTONS libraries (`PARTONSLIBS`) and includes
(`PARTONSINCS`). Once this is done, typing `make` will compile the
codes producing a number of executables. A more detailed
description of the output of each executable can be found in the
corresponding source code.

## References

If you use this code, please make sure to cite the paper of the
calculation:

- V. Bertone, *Matching generalised transverse-momentum-dependent distributions onto generalised parton distributions at one loop*, [arXiv:2207.09526](https://arxiv.org/pdf/2207.09526.pdf),

as well the publications of the dependency softwares:

- B. Berthou et al., *PARTONS: PARtonic Tomography Of Nucleon
  Software : A computing framework for the phenomenology of
  Generalized Parton Distributions*
  [arXiv:1512.06174](https://arxiv.org/pdf/1512.06174.pdf) for
  `PARTONS`,
- V. Bertone, S. Carrazza, and J. Rojo, *APFEL: A PDF Evolution
  Library with QED corrections*,
  [arXiv1310.1394](https://arxiv.org/pdf/1310.1394.pdf) and
  V. Bertone, *APFEL++: A new PDF evolution library in C++*
  [arXiv:1708.00911](https://arxiv.org/pdf/1708.00911.pdf) for
  `APFEL++`
- A. Bacchetta et al., *Transverse-momentum-dependent parton
  distributions up to N3LL from Drell-Yan data*,
  [arXiv:1912.07550](https://arxiv.org/pdf/1912.07550.pdf) for
  `NangaParbat`.

## Contact

Questions and/or suggestions can be addressed to:

- Valerio Bertone: valerio.bertone@cern.ch
