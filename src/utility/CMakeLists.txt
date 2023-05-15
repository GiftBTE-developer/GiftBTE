#add_library(BTEutility utility.cpp "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESynthetic/utility.h")
add_library(BTEutility  utility.cpp )
# target_include_directories(BTEutility PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes"   "${PETSC_INCLUDES}" "${HYPRE_INC_PATH}")
target_link_libraries(BTEutility PUBLIC  MPI::MPI_CXX "${HYPRE_LIBRARIES}")
