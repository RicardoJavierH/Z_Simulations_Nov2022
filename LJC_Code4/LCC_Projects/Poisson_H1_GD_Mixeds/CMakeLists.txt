file(GLOB header *.h)
file(GLOB sources *.cpp)

set(main_SRCS ${sources} ${headers} ${LCCCommon_SRC} ${LCCMaterial_SRC} ${LCCAdaptive_SRC} ${LCCEstimator_SRC})

add_executable(Simulator_H1_GD_Mixeds ${main_SRCS})

target_link_libraries(Simulator_H1_GD_Mixeds ${PZ_LIB})

if(WIN32)
	target_link_libraries(Simulator_H1_GD_Mixeds ${PTHREAD_LIB})
endif()