add_library(TransientFourier   TransientFourier.cpp )
target_include_directories(TransientFourier PUBLIC "/home/yuehu/Downloads" "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESolver"
        PRIVATE "${PETSC_INCLUDES}" "${HYPRE_INC_PATH}")

target_link_libraries(TransientFourier PUBLIC BTEBand BTEBoundaryCondition BTEMesh BTEAngle BTEutility
        PRIVATE "${PETSC_LIBRARIES}" MPI::MPI_CXX "${HYPRE_LIB_PATH}")