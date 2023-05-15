#add_library(BTEAngle BTEAngle.cpp "${StaticBTESolver_SOURCE_DIR}/includes/StaticBTESynthetic/BTEAngle.h")
add_library(BTEAngle  BTEAngle.cpp)
# target_include_directories(BTEAngle PUBLIC "${StaticBTESolver_SOURCE_DIR}/includes")
target_link_libraries(BTEAngle PUBLIC BTEutility)
# target_compile_features(BTEAngle PUBLIC cxx_std_14)



