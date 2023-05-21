find_package(MPI REQUIRED)
find_package(PETSc REQUIRED)

if(USE_GPU)
    find_package(CUDAToolkit REQUIRED)
    enable_language(CUDA)
    set_source_files_properties(BTEMain.cpp PROPERTIES LANGUAGE CUDA)   
endif()

include_directories("${StaticBTESolver_SOURCE_DIR}/includes" "${PETSC_INCLUDES}")

add_subdirectory(BTEAngle)
add_subdirectory(BTEBand)
add_subdirectory(BTEBoundaryCondition)
add_subdirectory(BTEMesh)

add_subdirectory(StaticBTESynthetic)
add_subdirectory(StaticFourier)
add_subdirectory(TransientFourier)
add_subdirectory(TransientBTE)
add_subdirectory(Solution)
add_subdirectory(utility)

if(USE_GPU)
    add_executable(BTE_GPU BTEMain.cpp)
else()
    add_executable(BTE_CPU BTEMain.cpp)
endif()

if(USE_GPU)
    target_link_libraries(BTE_GPU PUBLIC StaticBTESynthetic MPI::MPI_CXX  "${PETSC_LIBRARIES}")

        if(USE_TIME)
            set_target_properties(BTE_GPU PROPERTIES COMPILE_DEFINITIONS "USE_GPU;USE_TIME")
        else()
            set_target_properties(BTE_GPU PROPERTIES COMPILE_DEFINITIONS "USE_GPU")
        endif()

else()
    target_link_libraries(BTE_CPU PUBLIC StaticBTESynthetic StaticFourier TransientFourier TransientBTE Solution MPI::MPI_CXX "${PETSC_LIBRARIES}")

    if(USE_TIME)
        if(USE_MRM)
            set_target_properties(BTE_CPU PROPERTIES COMPILE_DEFINITIONS "USE_TIME;USE_MRM")
        else()
            set_target_properties(BTE_CPU PROPERTIES COMPILE_DEFINITIONS "USE_TIME")
        endif()
    else()
        if(USE_MRM)
            set_target_properties(BTE_CPU PROPERTIES COMPILE_DEFINITIONS "USE_MRM")     
        endif()
    endif()

endif()



