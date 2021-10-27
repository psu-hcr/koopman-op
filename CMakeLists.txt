cmake_minimum_required(VERSION 2.8.3)
project(iiwa_koopman)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Armadillo REQUIRED)
find_package(OpenCV REQUIRED)
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
set(EXECUTABLE_OUTPUT_PATH build/)
include_directories(${ARMADILLO_INCLUDE_DIRS} ${OpenCV_INCLUDE_DIRS} ${OpenMP_INCLUDE_DIRS} "/home/kzf5356/iiwa_ros/koopman-op")

set(EXECUTABLE_OUTPUT_PATH build/)

add_executable(errsac_cp errsac_cartpend.cpp)
target_link_libraries(errsac_cp ${ARMADILLO_LIBRARIES})

add_executable(errsac_koopcp errsac_koopman_cp.cpp)
target_link_libraries(errsac_koopcp ${ARMADILLO_LIBRARIES})

add_executable(sac_koopquad errsac_koopquad.cpp)
target_link_libraries(sac_koopquad ${ARMADILLO_LIBRARIES})

#add_executable(ergsac_cartpend ergsac_cartpend.cpp)
#target_link_libraries(ergsac_cartpend ${ARMADILLO_LIBRARIES})

#add_executable(sac_doubleint sac_doubleint.cpp)
#target_link_libraries(sac_doubleint ${ARMADILLO_LIBRARIES})

#add_executable(ergsac_doubleint ergsac_doubleint.cpp)
#target_link_libraries(ergsac_doubleint ${ARMADILLO_LIBRARIES} ${OpenCV_LIBS} ${OpenMP_LIBS})

#add_executable(ergsac_dilincoln ergsac_dilincoln.cpp)
#target_link_libraries(ergsac_dilincoln ${ARMADILLO_LIBRARIES} ${OpenCV_LIBS} ${OpenMP_LIBS})

#add_executable(ergsac_mda ergsac_mda.cpp)
#target_link_libraries(ergsac_mda ${ARMADILLO_LIBRARIES} ${OpenCV_LIBS} ${OpenMP_LIBS})

#add_executable(dklsac_doubleint dklsac_doubleint.cpp)
#target_link_libraries(dklsac_doubleint ${ARMADILLO_LIBRARIES} ${OpenMP_LIBS})

#add_executable(dklsac_di_image dklsac_di_image.cpp)
#target_link_libraries(dklsac_di_image ${ARMADILLO_LIBRARIES} ${OpenMP_LIBS} ${OpenCV_LIBS})

#add_executable(dklmda_di_image dklmda_di_image.cpp)
#target_link_libraries(dklmda_di_image ${ARMADILLO_LIBRARIES} ${OpenMP_LIBS} ${OpenCV_LIBS})

#add_executable(dklmda_si_image dklmda_si_image.cpp)
#target_link_libraries(dklmda_si_image ${ARMADILLO_LIBRARIES} ${OpenMP_LIBS} ${OpenCV_LIBS})

#add_executable(edgedetect edgedetect.cpp)
#target_link_libraries(edgedetect ${ARMADILLO_LIBRARIES} ${OpenCV_LIBS} ${OpenMP_LIBS})

#add_executable(cleanfiles filesplit.cpp)
#target_link_libraries(cleanfiles ${ARMADILLO_LIBRARIES} ${OpenCV_LIBS} ${OpenMP_LIBS})
