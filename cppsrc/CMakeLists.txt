FIND_PACKAGE(BLAS)
FIND_PACKAGE(LAPACK)

ENABLE_LANGUAGE(Fortran)
INCLUDE(FortranCInterface)

include_directories($ENV{HOME}/usr/include/lapackpp)
link_directories($ENV{HOME}/usr/lib)

add_executable(Snake 
snake.cpp
blocham.cpp        
public.cpp          
dmrg.cpp        
gqnmat.cpp              
block.cpp                   
dtmat.cpp           
site.cpp    
supblock.cpp  
supblock_iDMRG.cpp
supblock_tDMRG.cpp
supblock_FTDMRG.cpp
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



