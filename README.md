# Snake DMRG (Snake)

![](1.png)

Snake is a research DMRG code (C++) plus MATLAB helpers for generating model files.
It supports infinite DMRG (iDMRG), finite DMRG (fDMRG), and adaptive real-time evolution (tDMRG) for 1D / impurity-style problems.

The intended workflow is:
1) use MATLAB to generate a `model/` directory (binary operators + Hamiltonian terms),
2) run the `Snake` executable in the same folder,
3) read results from `results/`.

## Repository layout

- `cppsrc/`: C++ implementation and build scripts
- `*.m`: MATLAB scripts for generating inputs and plotting outputs
- `README.md`, `1.png`: this documentation and a screenshot/figure

## Build the `Snake` executable

### Dependencies

The C++ code is written against:
- BLAS + LAPACK
- ARPACK
- LAPACK++ (`lapackpp`) headers and library
- A C++ compiler and a Fortran compiler (CMake enables Fortran for LAPACK/ARPACK linkage)

### CMake build (recommended)

From the repo root:

```bash
mkdir -p build
cd build
cmake ..
cmake --build . -j
```

This should produce `build/cppsrc/Snake`.

#### LAPACK++ include/library paths

`cppsrc/CMakeLists.txt` currently assumes LAPACK++ is installed under:
- headers: `$HOME/usr/include/lapackpp`
- libs: `$HOME/usr/lib`

If your LAPACK++ install lives elsewhere, either:
- edit `cppsrc/CMakeLists.txt`, or
- set up symlinks so the above paths exist.

## Quick start (run the included MATLAB example)

The scripts in this repo generate a particular impurity + spinless-fermion chain example (see `gen_model.m`).

1) Build `Snake` (see above).
2) In MATLAB (from the repo root), generate a run folder:

```matlab
% Example parameters: gen_model(Omega, omegaratio, Lambda, BathL, z)
gen_model(0.05, 1.0, 2.0, 50, 1);
```

This creates a directory named like `ratio...Lambda...BathL...z.../` containing:
- `Snake` (copied from `./build/cppsrc/Snake`)
- `model/` (binary input files for the C++ code)

3) Run Snake inside the generated folder:

```bash
cd ratio1Lambda2BathL50z1   # name depends on your parameters
./Snake
```

4) Inspect outputs in `results/` (see below).

Tip: `Snake` reads `./model/...` relative to your current working directory, so you can also run a binary from elsewhere (e.g. `../../build/cppsrc/Snake`) as long as you `cd` into the run folder first.

## Input contract: the `model/` directory

The C++ executable uses fixed relative paths and expects a `./model/` directory next to the executable.
The files are written by `savepara.m` and read by the C++ code:

- `model/problemparmeters.dat` (note the spelling)
  - `int32`: chain length `L`
  - `int32`: target good quantum number `TGQN` (used as a total particle-number-like constraint)
- `model/site_base.dat`: local basis and good-quantum-number information per site (see `sitebase2file.m`)
- `model/site_operators.dat`: on-site operators per site (see `savepara.m`)
- `model/Hfac.dat`: real-valued model parameters (see `savepara.m`)
  - hopping terms, on-site energies, nearest-neighbor interaction terms
- `model/HC.dat`: two-site Hamiltonians for building the starting state (`H0` terms in `gen_model.m`)
- `model/rt_T0.dat`: time-evolution operators for the “bulk” two-site gates
- `model/rt_H1_T0.dat`: time-dependent evolution operators for the impurity + first site

Binary format notes (from `mat2file.m` and the C++ readers in `cppsrc/public.h`):
- integers are written as `int32`
- real scalars are written as `real*8` (IEEE-754 double)
- matrices are written with dimensions followed by raw column-major data (Fortran/MATLAB order)
- complex matrices are written as `(real, imag)` double pairs per element

## Outputs

`Snake` creates (and overwrites) directories on each run:
- `results/`: kept after the program exits
- `data/`: intermediate blocks; deleted on exit (and also cleared at startup)

Common output files produced by tDMRG (see `cppsrc/SuperChain_tDMRG.cpp`):
- `results/sigmaz_t.dat`: time trace of a local observable (named `sigma_z(t)` in the code)
- `results/vonneumannentropy_t.dat`: von Neumann entropies during evolution
- `results/steperror_t.dat`: placeholder for step-by-step diagnostics (currently created but not written to)
- `results/rdm.dat`: reduced density matrix diagnostics (used by `plotspectral.m`)

Some provided plotting scripts (run from the parent directory that contains your `ratio...` folders):
- `plotsigmaz.m`
- `plotspectral.m`
- `plotvisibility.m`


## Acknowledgements and citation

This program has been developed under the supervision of Prof. Tao Xiang and Prof. Jan von Delft.
Thanks to Shijie Hu, Honggang Luo, Shaojing Qin, Jize Zhao, Hantao Lu, Zhihui Wang and Shuming Li for either direct contribution or discussions during the development of the program.

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
