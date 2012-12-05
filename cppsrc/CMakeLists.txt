FIND_PACKAGE(BLAS)
FIND_PACKAGE(LAPACK)

ENABLE_LANGUAGE(Fortran)
INCLUDE(FortranCInterface)

include_directories($ENV{HOME}/usr/include/lapackpp)
link_directories($ENV{HOME}/usr/lib)

add_executable(Snake 
snake.cpp
InfiniteDMRG.cpp
FiniteDMRG.cpp
AdaptiveTimeDependentDMRG.cpp
FiniteTemperatureDMRG.cpp
ChainHamiltonian.cpp        
public.cpp          
dmrg.cpp        
gqnmat.cpp              
Chain.cpp                   
dtmat.cpp           
site.cpp    
SuperChain.cpp
SuperChain_iDMRG.cpp
SuperChain_tDMRG.cpp
SuperChain_FTDMRG.cpp
gqnbase.cpp
gqn.cpp
setting.h
)

target_link_libraries(Snake
lapackpp
arpack
${BLAS_LIBRARIES}
${LAPACK_LIBRARIES}
)



