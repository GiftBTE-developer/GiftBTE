#add_library(BTEBand BTEBand.cpp "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESynthetic/BTEBand.h")
add_library(BTEBand  BTEBand.cpp )
# target_include_directories(BTEBand PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes"  )
target_link_libraries(BTEBand PUBLIC BTEutility)
# target_compile_features(BTEBand PUBLIC cxx_std_14)


