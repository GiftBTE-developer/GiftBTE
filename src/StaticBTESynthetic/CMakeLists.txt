

add_library(StaticBTESynthetic StaticBTESynthetic.cpp StaticBTESyntheticSolver.cpp)


target_link_libraries(StaticBTESynthetic PUBLIC BTEBand BTEBoundaryCondition BTEMesh BTEAngle BTEutility
        PRIVATE MPI::MPI_CXX "${PETSC_LIBRARIES}" "${HYPRE_LIBRARIES}")
   

