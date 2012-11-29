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
gqnbase.cpp  
gqnmat.cpp              
block.cpp                   
dtmat.cpp       
gqn.cpp      
site.cpp    
supblock.cpp  
supblock_iDMRG.cpp
supblock_tDMRG.cpp
supblock_FTDMRG.cpp
)

target_link_libraries(Snake
lapackpp
arpack
${BLAS_LIBRARIES}
${LAPACK_LIBRARIES}
)



