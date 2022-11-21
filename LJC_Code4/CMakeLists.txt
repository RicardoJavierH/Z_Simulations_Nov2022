# @file SBFem/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required (VERSION 3.0.0)

project (LJC_Codes)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#add_definitions(-D_AUTODIFF)

IF(APPLE)
	# Mac OS X specific code
	add_definitions(-DMACOSX)
	set(MACOSX MACOSX)
	FIND_LIBRARY(ACCELERATE_LIB Accelerate)
	link_libraries(${ACCELERATE_LIB})
	#   SET(OperatingSystem "Mac OS X")
ENDIF(APPLE)

IF (WIN32)
	# Prevents timespec redefinition problem with Visual Studio 2015
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_STRUCT_TIMESPEC")
    if (NOT CMAKE_GENERATOR MATCHES "Ninja")
    	set(CMAKE_INSTALL_PREFIX "C:/Program Files/PZ" CACHE PATH "Installation path, prepended onto install directories." FORCE) 
	endif()
ELSE (WIN32)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")#this flag is for preventing errors when casting from double to float
	if ( NOT APPLE AND NOT CMAKE_BUILD_TYPE )
		set (CMAKE_BUILD_TYPE Release CACHE STRING "Options: Release Debug None MinSizeRel RelWithDebInfo" FORCE) 
	endif()
ENDIF()

option(USING_MKL "Whether the MKL library will be linked in" OFF)
if(USING_MKL)
	if (WIN32)
		set(MKL_ROOT "C:/Program Files (x86)/IntelSWTools/compilers_and_libraries/windows")
		find_path(MKL_INCLUDE NAMES mkl.h PATHS ${MKL_ROOT}/mkl/include)
		find_library(MKL_LIB_INTEL NAMES mkl_rt.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
		find_library(MKL_LIB_CORE NAMES mkl_core.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
		find_library(MKL_LIB_THREAD NAMES mkl_intel_thread.lib PATHS ${MKL_ROOT}/mkl/lib/intel64_win)
                find_library(COMPOSER_OMP NAMES libiomp5md.lib PATHS ${MKL_ROOT}/compiler/lib/intel64_win)
	else()
		find_path(MKL_INCLUDE NAMES mkl.h PATHS ${SEARCH_DIRS} /opt/intel/mkl/include /softwares/intel/mkl/include)
		find_library(MKL_LIB_INTEL NAMES libmkl_intel_lp64.so libmkl_intel_lp64.dylib PATHS
				${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64)
		find_library(MKL_LIB_CORE NAMES libmkl_core.so libmkl_core.dylib PATHS ${SEARCH_DIRS}
				/opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
		find_library(MKL_LIB_THREAD NAMES libmkl_intel_thread.so libmkl_intel_thread.dylib PATHS
				${SEARCH_DIRS} /opt/intel/mkl/lib /opt/intel/mkl/lib/intel64/ /softwares/intel/mkl/lib /softwares/intel/mkl/lib/intel64/)
		find_library(COMPOSER_OMP NAMES libiomp5.so libiomp5.dylib PATHS ${SEARCH_DIRS}
				/opt/intel/composer_xe/compiler/lib
				/opt/intel/composerxe/lib/intel64
				/opt/intel/lib /softwares/intel/lib/intel64
                                /opt/intel/compilers_and_libraries/linux/lib/intel64_lin )
	endif()

	if(MKL_INCLUDE-NOTFOUND)
                set (MKL_INCLUDE "" CACHE PATH "Directory where mkl.h can be found")
        else()
                #include_directories(${MKL_INCLUDE})
        endif()
        if(MKL_LIB-NOTFOUND)
                set (MKL_LIB_INTEL "" CACHE PATH "Directory where the mkl library can be found")
        else()
                link_libraries(${MKL_LIB_INTEL})
                link_libraries(${MKL_LIB_CORE})
                link_libraries(${MKL_LIB_THREAD})
                link_libraries(${COMPOSER_OMP})
        endif()

	include_directories(${MKL_INCLUDE})

endif(USING_MKL)


include_directories(${PZ_INCLUDE_DIRS})

# Location to search files or libraries
set(SEARCH_DIRS ENV PATH /usr/local/opt/ /usr/local/opt/qt5/lib/cmake/ /usr/include /usr/local/include /usr/local/include/boost /usr/local/include/log4cxx /usr/local/include/tbb /usr/lib /usr/lib64 /usr/local/lib ./../Libraries ./../Libraries/pz/include ./../Libraries/pz/lib )


#-----
#enabling to use LOG4CXX library
option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
if(USING_LOG4CXX)
	#Adding a variable to hold this definition
	set(LOG4CXX LOG4CXX)
	find_path(LOG4CXX_INCLUDE log4cxx PATHS ${SEARCH_DIRS} ${PROJECT_SOURCE_DIR}/../log4cxx ${PROJECT_SOURCE_DIR}/../log4cxx_src/include ${PROJECT_SOURCE_DIR}/../Libraries/log4cxx_src ${PROJECT_SOURCE_DIR}/../Libraries/log4cxx_src/include)
	find_library(LOG4CXX_LIB_RELEASE
			NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a
			PATHS ${SEARCH_DIRS} ./../log4cxx_src ./../log4cxx_src/lib ./../Libraries/log4cxx_src ./../Libraries/log4cxx_src/lib ./../Libraries/lib ./../Libraries/lib/Release
			)
	find_library(LOG4CXX_LIB_DEBUG
			NAMES liblog4cxx.dylib log4cxx.lib liblog4cxx.so liblog4cxx.a       # Change name here if debug lib has different name.
			PATHS ${SEARCH_DIRS} ./../log4cxx_src ./../log4cxx_src/lib ./../Libraries/log4cxx_src ./../Libraries/log4cxx_src/lib ./../Libraries/lib ./../Libraries/lib/Debug
			)
	set(LOG4CXX_LIB
			debug ${LOG4CXX_LIB_DEBUG}
			optimized ${LOG4CXX_LIB_RELEASE}
			)
	if(LOG4CXX_INCLUDE-NOTFOUND)
		set (LOG4CXX_INCLUDE "" CACHE PATH "Directory where log4cxx header files can be found")
	else()
		include_directories(${LOG4CXX_INCLUDE})
	endif()

	if(LOG4CXX_LIB-NOTFOUND)
		set (LOG4CXX_LIB "" CACHE PATH "Directory where the log4cxx library can be found")
	endif()
endif(USING_LOG4CXX)



#enabling pz library
find_path(PZ_INCLUDEINI NAMES pz_config.h PATHS ${SEARCH_DIRS} C:/Z_Simulations/externallibs/pzlib/include D:/Z_Simulations/externallibs/pzlib/include)
set(PZ_INCLUDE "${PZ_INCLUDEINI}" "${PZ_INCLUDEINI}/Analysis" "${PZ_INCLUDEINI}/Common" "${PZ_INCLUDEINI}/Frontal" "${PZ_INCLUDEINI}/External" "${PZ_INCLUDEINI}/External/FAD" "${PZ_INCLUDEINI}/External/FAD/fad" "${PZ_INCLUDEINI}/External/FAD/TinyFadET" "${PZ_INCLUDEINI}/External/FAD/TinyFad" "${PZ_INCLUDEINI}/Geom" "${PZ_INCLUDEINI}/Integral" "${PZ_INCLUDEINI}/LinearSolvers" "${PZ_INCLUDEINI}/Material" "${PZ_INCLUDEINI}/Material/REAL" "${PZ_INCLUDEINI}/Material/Complex" "${PZ_INCLUDEINI}/Matrix" "${PZ_INCLUDEINI}/Mesh" "${PZ_INCLUDEINI}/Multigrid" "${PZ_INCLUDEINI}/Optimization" "${PZ_INCLUDEINI}/Post" "${PZ_INCLUDEINI}/Pre" "${PZ_INCLUDEINI}/Refine" "${PZ_INCLUDEINI}/Save" "${PZ_INCLUDEINI}/Shape" "${PZ_INCLUDEINI}/StrMatrix" "${PZ_INCLUDEINI}/SubStruct" "${PZ_INCLUDEINI}/Topology" "${PZ_INCLUDEINI}/Util")
find_library(PZ_LIB NAMES libpz.lib pz.lib libpz.a libpz-debug.dylib PATHS ${SEARCH_DIRS}  ${PROJECT_SOURCE_DIR}/../Libraries/pz/lib C:/Z_Simulations/externallibs/pzlib/lib D:/Z_Simulations/externallibs/pzlib/lib)
if(PZ_INCLUDE-NOTFOUND)
	set(PZ_INCLUDE "" CACHE PATH "Directory where pz.h can be found")
endif()
if(PZ_LIB-NOTFOUND)
	set(PZ_LIBexit "" CACHE PATH "Directory where the pz library can be found")
endif()
include_directories(${PZ_INCLUDE})

#enabling pthread library
if(WIN32)
find_path(PTHREAD_INCLUDE NAMES PThread.h pthread.h PATHS ${SEARCH_DIRS} C:/Z_Simulations/externallibs/pthread/include D:/Z_Simulations/externallibs/pthread/include)
find_library(PTHREAD_LIB NAMES libpthread.dylib pthread.lib pthreadVC2.lib libPThread.so libPThread.a libpthread.so libpthread.a PATHS ${SEARCH_DIRS} C:/Z_Simulations/externallibs/pthread/lib D:/Z_Simulations/externallibs/pthread/lib)
if(PTHREAD_INCLUDE-NOTFOUND)
	set(PTHREAD_INCLUDE "" CACHE PATH "Directory where pthread.h can be found")
endif()
if(PTHREAD_LIB-NOTFOUND)
	set(PTHREAD_LIBexit "" CACHE PATH "Directory where the pthread library can be found")
endif()
include_directories(${PTHREAD_INCLUDE})
endif()

if(WIN32)
	#add a define saying it is a VC compiler
	set(VC "#Will add a define VC on pz_config.h")
	add_definitions(-DVC)
	#define _USE_MATH_DEFINES for example: M_PI constant definitions.
	add_definitions(-D_USE_MATH_DEFINES)
	#disabling VC warnings
	add_definitions(-D_SCL_SECURE_NO_WARNINGS)
	add_definitions(-D_CRT_SECURE_NO_WARNINGS)
	#disabling All warnings
	#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /W0")
	#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /W0")
	#disabling STL warning for std::_Vector_alloc when LOG4CXX is enabled
	if (USING_LOG4CXX)
		set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /wd4251")
		set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /wd4251")
	endif()
	#define use of pthread static lib.
	#enabling /bigobj
	add_definitions("/bigobj")

	FOREACH(FLAG_TYPE EXE MODULE SHARED)
		STRING (REPLACE "INCREMENTAL:YES" "INCREMENTAL" FLAG_TMP "${CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG}")
		STRING (REPLACE "INCREMENTAL:NO" "INCREMENTAL" FLAG_TMP ${FLAG_TMP})
		STRING (REPLACE "INCREMENTAL" "INCREMENTAL:NO" FLAG_TMP ${FLAG_TMP})
		STRING (REPLACE "/EDITANDCONTINUE" "" FLAG_TMP ${FLAG_TMP})
		SET(CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG "${FLAG_TMP}" CACHE STRING "Overriding default debug ${FLAG_TYPE} linker flags." FORCE)
		MARK_AS_ADVANCED (CMAKE_${FLAG_TYPE}_LINKER_FLAGS_DEBUG)
	ENDFOREACH ()

endif()


add_subdirectory(LCCCommon)
add_subdirectory(LCCMaterials)
#add_subdirectory(LCCAnalysis)
add_subdirectory(LCCAdaptive)
add_subdirectory(LCCEstimator)

include_directories(LCCCommon)
include_directories(LCCMaterials)
include_directories(LCCAnalysis)
include_directories(LCCAdaptive)
include_directories(LCCEstimator)

add_subdirectory(LCC_Projects)
