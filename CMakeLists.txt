cmake_minimum_required( VERSION 2.8 )

project( SPRUCE )

set( merge_src
  src/merge.cpp
  src/utils.cpp
  src/matrix.cpp
  src/realmatrix.cpp
  src/tensor.cpp
  src/realtensor.cpp
  src/perfectphylomatrix.cpp
  src/perfectphylotree.cpp
  src/perfectphylograph.cpp
  src/rootedcladisticancestrygraph.cpp
  src/statetree.cpp
  src/solution.cpp
  src/solutionset.cpp
)

set( merge_hdr
  src/utils.h
  src/matrix.h
  src/realmatrix.h
  src/tensor.h
  src/realtensor.h
  src/perfectphylomatrix.h
  src/perfectphylotree.h
  src/perfectphylograph.h
  src/rootedcladisticancestrygraph.h
  src/statetree.h
  src/solution.h
  src/solutionset.h
)

set( ccf_src
  src/ccf.cpp
  src/utils.cpp
  src/matrix.cpp
  src/realmatrix.cpp
  src/tensor.cpp
  src/realtensor.cpp
  src/perfectphylomatrix.cpp
  src/perfectphylotree.cpp
  src/perfectphylograph.cpp
  src/rootedcladisticancestrygraph.cpp
  src/statetree.cpp
  src/solution.cpp
  src/solutionset.cpp
)

set( ccf_hdr
  src/utils.h
  src/matrix.h
  src/realmatrix.h
  src/tensor.h
  src/realtensor.h
  src/perfectphylomatrix.h
  src/perfectphylotree.h
  src/perfectphylograph.h
  src/rootedcladisticancestrygraph.h
  src/statetree.h
  src/solution.h
  src/solutionset.h
)

set( visualize_src
  src/visualize.cpp
  src/utils.cpp
  src/matrix.cpp
  src/realmatrix.cpp
  src/tensor.cpp
  src/realtensor.cpp
  src/perfectphylomatrix.cpp
  src/perfectphylotree.cpp
  src/perfectphylograph.cpp
  src/rootedcladisticancestrygraph.cpp
  src/statetree.cpp
  src/solution.cpp
  src/solutiongraph.cpp
  src/solutionset.cpp
)

set( visualize_hdr
  src/utils.h
  src/matrix.h
  src/realmatrix.h
  src/tensor.h
  src/realtensor.h
  src/perfectphylomatrix.h
  src/perfectphylotree.h
  src/perfectphylograph.h
  src/rootedcladisticancestrygraph.h
  src/statetree.h
  src/solution.h
  src/solutiongraph.h
  src/solutionset.h
)

set( concordance_src
  src/concordance.cpp
  src/utils.cpp
  src/matrix.cpp
  src/realmatrix.cpp
  src/tensor.cpp
  src/realtensor.cpp
  src/perfectphylomatrix.cpp
  src/perfectphylotree.cpp
  src/perfectphylograph.cpp
  src/rootedcladisticancestrygraph.cpp
  src/statetree.cpp
  src/solution.cpp
  src/solutionset.cpp
)

set( concordance_hdr
  src/utils.h
  src/matrix.h
  src/realmatrix.h
  src/tensor.h
  src/realtensor.h
  src/perfectphylomatrix.h
  src/perfectphylotree.h
  src/perfectphylograph.h
  src/rootedcladisticancestrygraph.h
  src/statetree.h
  src/solution.h
  src/solutionset.h
)

set( enumerate_src
  src/enum.cpp
  src/bronkerbosch.cpp
  src/compatibilitygraph.cpp
  src/stategraph.cpp
  src/character.cpp
  src/charactermatrix.cpp
  src/utils.cpp
  src/matrix.cpp
  src/realmatrix.cpp
  src/tensor.cpp
  src/realtensor.cpp
  src/rootedcladisticancestrygraph.cpp
  src/rootedcladisticnoisyancestrygraph.cpp
  src/perfectphylomatrix.cpp
  src/perfectphylotree.cpp
  src/perfectphylograph.cpp
  src/statetree.cpp
  src/rootedcladisticenumeration.cpp
  src/rootedcladisticnoisyenumeration.cpp
  src/solution.cpp
  src/solutionset.cpp
  src/cnaenumerate.cpp
  src/noisycnaenumerate.cpp
)

set( enumerate_hdr
  src/bronkerbosch.h
  src/compatibilitygraph.h
  src/config.h
  src/stategraph.h
  src/character.h
  src/charactermatrix.h
  src/utils.h
  src/matrix.h
  src/realmatrix.h
  src/tensor.h
  src/realtensor.h
  src/rootedcladisticancestrygraph.h
  src/rootedcladisticnoisyancestrygraph.h
  src/perfectphylomatrix.h
  src/perfectphylotree.h
  src/perfectphylograph.h
  src/statetree.h
  src/rootedcladisticenumeration.h
  src/rootedcladisticnoisyenumeration.h
  src/solution.h
  src/solutionset.h
  src/cnaenumerate.h
  src/noisycnaenumerate.h
)

set( compatiblestatetrees_src
  src/compatiblestatetrees.cpp
  src/stategraph.cpp
  src/character.cpp
  src/charactermatrix.cpp
  src/utils.cpp
  src/matrix.cpp
  src/realmatrix.cpp
  src/tensor.cpp
  src/realtensor.cpp
  src/rootedcladisticancestrygraph.cpp
  src/rootedcladisticnoisyancestrygraph.cpp
  src/perfectphylomatrix.cpp
  src/perfectphylotree.cpp
  src/perfectphylograph.cpp
  src/statetree.cpp
  src/rootedcladisticenumeration.cpp
  src/rootedcladisticnoisyenumeration.cpp
  src/solution.cpp
  src/solutionset.cpp
)

set( compatiblestatetrees_hdr
  src/config.h
  src/stategraph.h
  src/character.h
  src/charactermatrix.h
  src/utils.h
  src/matrix.h
  src/realmatrix.h
  src/tensor.h
  src/realtensor.h
  src/rootedcladisticancestrygraph.h
  src/rootedcladisticnoisyancestrygraph.h
  src/perfectphylomatrix.h
  src/perfectphylotree.h
  src/perfectphylograph.h
  src/statetree.h
  src/rootedcladisticenumeration.h
  src/rootedcladisticnoisyenumeration.h
  src/solution.h
  src/solutionset.h
)

set( cliques_src
  src/cliques.cpp
  src/compatibilitygraph.cpp
  src/bronkerbosch.cpp
  src/stategraph.cpp
  src/character.cpp
  src/charactermatrix.cpp
  src/utils.cpp
  src/matrix.cpp
  src/realmatrix.cpp
  src/tensor.cpp
  src/realtensor.cpp
  src/rootedcladisticancestrygraph.cpp
  src/rootedcladisticnoisyancestrygraph.cpp
  src/perfectphylomatrix.cpp
  src/perfectphylotree.cpp
  src/perfectphylograph.cpp
  src/statetree.cpp
  src/rootedcladisticenumeration.cpp
  src/rootedcladisticnoisyenumeration.cpp
  src/solution.cpp
  src/solutionset.cpp
)

set( cliques_hdr
  src/config.h
  src/compatibilitygraph.h
  src/bronkerbosch.h
  src/stategraph.h
  src/character.h
  src/charactermatrix.h
  src/utils.h
  src/matrix.h
  src/realmatrix.h
  src/tensor.h
  src/realtensor.h
  src/rootedcladisticancestrygraph.h
  src/rootedcladisticnoisyancestrygraph.h
  src/perfectphylomatrix.h
  src/perfectphylotree.h
  src/perfectphylograph.h
  src/statetree.h
  src/rootedcladisticenumeration.h
  src/rootedcladisticnoisyenumeration.h
  src/solution.h
  src/solutionset.h
)

set( randomtree_src
  src/randomtree.cpp
  src/utils.cpp
  src/matrix.cpp
  src/realmatrix.cpp
  src/tensor.cpp
  src/realtensor.cpp
  src/rootedcladisticfullancestrygraph.cpp
  src/rootedcladisticsampler.cpp
  src/perfectphylomatrix.cpp
  src/perfectphylotree.cpp
  src/perfectphylograph.cpp
  src/statetree.cpp
  src/stategraph.cpp
  src/statetreessampler.cpp
)

set( randomtree_hdr
  src/utils.h
  src/matrix.h
  src/realmatrix.h
  src/tensor.h
  src/realtensor.h
  src/rootedcladisticfullancestrygraph.h
  src/rootedcladisticsampler.h
  src/perfectphylomatrix.h
  src/perfectphylotree.h
  src/perfectphylograph.h
  src/statetree.h
  src/stategraph.h
  src/statetreessampler.h
)

set( rank_hdr
  src/utils.h
  src/matrix.h
  src/realmatrix.h
  src/tensor.h
  src/realtensor.h
  src/perfectphylomatrix.h
  src/perfectphylotree.h
  src/perfectphylograph.h
  src/rootedcladisticancestrygraph.h
  src/statetree.h
  src/solution.h
  src/solutionset.h
)

set( rank_src
  src/rank.cpp
  src/utils.cpp
  src/matrix.cpp
  src/realmatrix.cpp
  src/tensor.cpp
  src/realtensor.cpp
  src/perfectphylomatrix.cpp
  src/perfectphylotree.cpp
  src/perfectphylograph.cpp
  src/rootedcladisticancestrygraph.cpp
  src/statetree.cpp
  src/solution.cpp
  src/solutionset.cpp
)

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

set( GitCommand1 "symbolic-ref HEAD 2> /dev/null | cut -b 12-" )
execute_process( COMMAND git symbolic-ref HEAD  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
COMMAND tail -c -7 COMMAND tr -d "\\n" OUTPUT_VARIABLE GitBranch )
set( GitCommand2 "log --pretty=format:\\\"%h\\\" -1" )
execute_process( COMMAND git log --pretty=format:%h WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
COMMAND head -n 1  COMMAND tr -d "\\n" OUTPUT_VARIABLE GitHashtag )
set( SPRUCE_VERSION \"${GitBranch}-${GitHashtag}\" )

configure_file( src/config.h.in ${PROJECT_SOURCE_DIR}/src/config.h )

add_definitions( -DIL_STD )

set( Boost_USE_STATIC_LIBS ON )
find_package( Boost REQUIRED COMPONENTS thread system filesystem )

include_directories( "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} )
link_directories( "${LIBLEMON_ROOT}/lib" )

if(NOT CMAKE_BUILD_TYPE)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set( CommonLibs
  m
  emon
  pthread
  ${Boost_LIBRARIES}
)

set( CplexLibs
  ilocplex
  cplex
  concert
)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

add_executable( enumerate ${enumerate_src} ${enumerate_hdr} )
target_link_libraries( enumerate ${CommonLibs} )

add_executable( compatiblestatetrees EXCLUDE_FROM_ALL ${compatiblestatetrees_src} ${compatiblestatetrees_hdr} )
target_link_libraries( compatiblestatetrees ${CommonLibs} )

add_executable( cliques ${cliques_src} ${cliques_hdr} )
target_link_libraries( cliques ${CommonLibs} )

add_executable( merge ${merge_src} ${merge_hdr} )
target_link_libraries( merge ${CommonLibs} )

add_executable( ccf ${ccf_src} ${ccf_hdr} )
target_link_libraries( ccf ${CommonLibs} )

add_executable( visualize ${visualize_src} ${visualize_hdr} )
target_link_libraries( visualize ${CommonLibs} )

add_executable( recall EXCLUDE_FROM_ALL ${concordance_src} ${concordance_hdr} )
target_link_libraries( recall ${CommonLibs} )

add_executable( rank ${rank_src} ${rank_hdr} )
target_link_libraries( rank ${CommonLibs} )

add_executable( randomtree EXCLUDE_FROM_ALL ${randomtree_src} ${randomtree_hdr} )
target_link_libraries( randomtree ${CommonLibs})
