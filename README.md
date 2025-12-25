# snake-dmrg
Snake is a highly optimized DMRG program. It can be used to calculate the ground state, real-time dynamics and finite temperature properties of quantum 1D or impurity systems. Abelian symmetry can be used to speed up the calculation.

A Matlab script is used to define the problem and generate operators and Hamiltonian. The c++ DMRG code can be used as a black box for simple problems.

The program has been tested on Linux, Unix and Mac OS X, but should also be possible to migrate to other platforms. 

This program has been developed under the supervision of Prof. Tao Xiang and Prof. Jan von Delft. Thank Shijie Hu, Honggang Luo, Shaojing Qin, Jize Zhao, Hantao Lu, Zhihui Wang and Shuming Li for either direct contribution or discussions during the development of the program.

If you use this code or build on these ideas, please cite:

```
@article{PhysRevB.79.115137,
  title = {Density matrix renormalization group study of a quantum impurity model with Landau-Zener time-dependent Hamiltonian},
  author = {Guo, Cheng and Weichselbaum, Andreas and Kehrein, Stefan and Xiang, Tao and von Delft, Jan},
  journal = {Phys. Rev. B},
  volume = {79},
  issue = {11},
  pages = {115137},
  numpages = {6},
  year = {2009},
  month = {Mar},
  publisher = {American Physical Society},
  doi = {10.1103/PhysRevB.79.115137},
  url = {https://link.aps.org/doi/10.1103/PhysRevB.79.115137}
}
```
