# Possible Future Improvements #

Here are the possible improvements I wish to see. However, as I no longer work in this field, I probably won't be able to implement them by myself alone. Please contact me if you are interested to implement them or have other ideas.

  * Redesign the class structure and better encapsulation.
  * Integrate Lapack++ into the code for easier installation process, or use Boost C++ Libraries instead.
  * Use the same input and output formats as those used in [ALPS](http://alps.comp-phys.org/mediawiki/index.php/Main_Page) (XML etc.).
  * Write test and examples.
  * Documentation.

The following is my sketch of the new class structure:

  * GoodQuantumNumber
  * Base
  * BlockMatrix

  * Operator:BlockMatrix
  * ChainOperator:Operator
  * SiteOperator:ChainOperator
  * Hamiltonian:Operator
  * ChainHamiltonian:Hamiltonian, ChainOperator
  * SiteHamiltonian:ChainHamiltonian
  * Chain
  * Site:Chain
  * WaveFunction:BlockMatrix
  * DensityMatrix:BlockMatrix
  * TruncationMatrix:BlockMatrix (or TransferMatrix)
  * MatrixProductState

  * MRG
  * DMRG:DMRG
  * DMRG:DMRG
  * DMRG:DMRG
  * TDMRG:DMRG