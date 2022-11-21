# Error Estimator source files used by targets of project ErrorEstimation
#

file(GLOB all_sources *.cpp *.h)

set(LCCEstimator_SRC ${all_sources} PARENT_SCOPE)
   
include_directories(.)