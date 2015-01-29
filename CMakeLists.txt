cmake_minimum_required( VERSION 2.8 )

project( Heinz )

set( BK_MaxFlow_Hdr
  src/externlib/maxflow-v3.01/graph.h
  src/externlib/maxflow-v3.01/block.h
)

set( Heinz_Hdr
  src/parser/parser.h
  src/parser/mwcsparser.h
  src/parser/stpparser.h
  src/parser/stppcstparser.h
  src/parser/identityparser.h
  src/parser/dimacsparser.h
  src/preprocessing/rule.h
  src/preprocessing/negdeg01.h
  src/preprocessing/negedge.h
  src/preprocessing/posedge.h
  src/preprocessing/negcircuit.h
  src/preprocessing/negdiamond.h
  src/preprocessing/negmirroredhubs.h
  src/preprocessing/negdominatedhubs.h
  src/preprocessing/posdeg01.h
  src/preprocessing/posdiamond.h
  src/preprocessing/shortestpath.h
  src/solver/spqrtree.h
  src/solver/blockcuttree.h
  src/mwcs.h
  src/analysis.h
  src/solver/enumsolverunrooted.h
  src/solver/impl/solverimpl.h
  src/solver/impl/solverunrootedimpl.h
  src/solver/impl/solverrootedimpl.h
  src/solver/solver.h
  src/solver/solverunrooted.h
  src/solver/solverrooted.h
  src/solver/impl/treesolverimpl.h
  src/solver/impl/treesolverunrootedimpl.h
  src/solver/impl/treesolverrootedimpl.h
  src/solver/impl/cplex_incumbent/incumbent.h
  src/solver/impl/cplex_incumbent/pcstincumbent.h
  src/solver/impl/cplex_cut/backoff.h
  src/solver/impl/cplex_cut/nodecut.h
  src/solver/impl/cplex_cut/nodecutuser.h
  src/solver/impl/cplex_cut/nodecutlazy.h
  src/solver/impl/cplex_cut/nodecutrooted.h
  src/solver/impl/cplex_cut/nodecutunrooted.h
  src/solver/impl/cplex_cut/bk_alg.h
  src/solver/impl/cplexsolverimpl.h
  src/solver/impl/cutsolverrootedimpl.h
  src/solver/impl/cutsolverunrootedimpl.h
  src/solver/impl/cplex_heuristic/heuristicrooted.h
  src/solver/impl/cplex_heuristic/heuristicunrooted.h
  src/solver/impl/cplex_branch/branch.h
  src/mwcsgraph.h
  src/mwcspreprocessedgraph.h
  src/utils.h
  src/mwcsgraphparser.h
)

set( Heinz_Src
  src/utils.cpp
)

set( Heinz_Monte_Carlo_Hdr
  src/solver/impl/treeheuristicsolverimpl.h
  src/solver/impl/treeheuristicsolverrootedimpl.h
  src/solver/impl/treeheuristicsolverunrootedimpl.h
)

set( Heinz_Monte_Carlo_Src 
  src/utils.cpp
)

set( Heinz_Print_Graph_Src
  src/printgraph.cpp
  src/utils.cpp
)

set( GitCommand1 "symbolic-ref HEAD 2> /dev/null | cut -b 12-" )
execute_process( COMMAND git symbolic-ref HEAD  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                 COMMAND tail -c -7 COMMAND tr -d "\\n" OUTPUT_VARIABLE GitBranch )
set( GitCommand2 "log --pretty=format:\\\"%h\\\" -1" )
execute_process( COMMAND git log --pretty=format:%h WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                 COMMAND head -n 1  COMMAND tr -d "\\n" OUTPUT_VARIABLE GitHashtag )
set( HEINZ_VERSION \"${GitBranch}-${GitHashtag}\" )

set( LIBOGDF_ROOT "$ENV{HOME}/ogdf" CACHE PATH "Additional search directory for OGDF library" )
set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set( MEMORY_LIMIT 10240 CACHE STRING "Memory limit in MB" )

if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
	set( CPLEX_AUTODETECTION_PREFIXES "$ENV{HOME}/ILOG;/opt/ibm/ILOG;$ENV{HOME}/Applications/IBM/ILOG" CACHE STRING "CPLEX automatic lookup directories" )
	foreach( CplexInstallPrefix ${CPLEX_AUTODETECTION_PREFIXES} )
		file( GLOB_RECURSE CplexLibraryPaths ${CplexInstallPrefix}/*/cplex/lib/*/libcplex.a )
		if( CplexLibraryPaths )
			list( LENGTH CplexLibraryPaths NumPath )
			if( NumPath GREATER 1 )
				message( WARNING "Found multiple CPLEX installations in ${CplexInstallPrefix}, automatically selecting the first one, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually if this is not the right selection" )
			else()
				message( "Found CPLEX's installation in ${CPlexInstallPrefix}" )
			endif()
			list( GET CplexLibraryPaths 0 CplexLibraryPath )

			set( MyRegex "${CplexInstallPrefix}/(.*)/cplex/?/lib/(.*)/libcplex.a" )
			string( REGEX MATCH "${MyRegex}" ValidatedCplexLibPath ${CplexLibraryPath} )
			if( NOT ValidatedCplexLibPath )
				message( WARNING "Unable to recognize ${CplexLibraryPath} as a valid CPLEX path" )
			else()
				string( REGEX REPLACE "${MyRegex}" "\\1" CplexHome ${ValidatedCplexLibPath} )
				string( REGEX REPLACE "${MyRegex}" "\\2" CplexArch ${ValidatedCplexLibPath} )

				set( CPLEX_INC_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/include"            CACHE PATH "CPLEX's include dir" )
				set( CPLEX_LIB_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/lib/${CplexArch}"   CACHE PATH "CPLEX's library dir" )
				set( CONCERT_INC_DIR "${CplexInstallPrefix}/${CplexHome}/concert/include"          CACHE PATH "Concert's include dir" )
        set( CONCERT_LIB_DIR "${CplexInstallPrefix}/${CplexHome}/concert/lib/${CplexArch}" CACHE PATH "Concert's library dir" )

				break()
			endif()
		else()
			message( "Unable to find CPLEX and Concert installation in ${CplexInstallPrefix}" )
		endif()
	endforeach()

	if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
		message( "Unable to autodeted CPLEX's and Concert's installation directories, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually" )
	else()
		message( "Autodetected CPLEX's include dir   : ${CPLEX_INC_DIR}" )
		message( "Autodetected CPLEX's library dir   : ${CPLEX_LIB_DIR}" )
		message( "Autodetected Concert's include dir : ${CONCERT_INC_DIR}" )
		message( "Autodetected Concert's library dir : ${CONCERT_LIB_DIR}" )
	endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libstdc++")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -g -ggdb -O0 -DDEBUG -DVERBOSE")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_PROFILE "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -O3 -g -ggdb -DNDEBUG")

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wno-unused-private-field -Wno-unused-variable -Wno-gnu-static-float-init" )
add_definitions( -DIL_STD )
configure_file( src/config.h.in ${PROJECT_SOURCE_DIR}/src/config.h )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
endif()

include_directories( ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} "${LIBLEMON_ROOT}/include" "${LIBOGDF_ROOT}/include" "src" "src/externlib" )
link_directories( ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} "${LIBLEMON_ROOT}/lib" "${LIBOGDF_ROOT}/lib" )

set( CommonLibs
	m
	ilocplex
	cplex
	concert
	pthread
	emon
	OGDF
)

add_executable( heinz-mc EXCLUDE_FROM_ALL src/mwcs-mc.cpp ${Heinz_Monte_Carlo_Src} ${Heinz_Monte_Carlo_Hdr} ${Heinz_Hdr} )
target_link_libraries( heinz-mc emon OGDF pthread )

add_executable( heinz src/mwcs.cpp ${Heinz_Src} ${Heinz_Hdr} )
target_link_libraries( heinz ${CommonLibs} )

add_executable( heinz_mwcs_mc src/dimacs/heinz_mwcs_mc.cpp ${Heinz_Monte_Carlo_Src} ${Heinz_Monte_Carlo_Hdr} ${Heinz_Hdr} )
target_link_libraries( heinz_mwcs_mc emon OGDF pthread )

add_executable( heinz_mwcs_dc src/dimacs/heinz_mwcs_dc.cpp ${Heinz_Src} ${Heinz_Hdr} )
target_link_libraries( heinz_mwcs_dc ${CommonLibs} )

add_executable( heinz_mwcs_no_dc src/dimacs/heinz_mwcs_no_dc.cpp ${Heinz_Src} ${Heinz_Hdr} )
target_link_libraries( heinz_mwcs_no_dc ${CommonLibs} )

add_executable( heinz_mwcs_no_pre src/dimacs/heinz_mwcs_no_pre.cpp ${Heinz_Src} ${Heinz_Hdr} )
target_link_libraries( heinz_mwcs_no_pre ${CommonLibs} )

add_executable( heinz_pcst_mc src/dimacs/heinz_pcst_mc.cpp ${Heinz_Monte_Carlo_Src} ${Heinz_Monte_Carlo_Hdr} ${Heinz_Hdr} )
target_link_libraries( heinz_pcst_mc emon OGDF pthread )

add_executable( heinz_pcst_dc src/dimacs/heinz_pcst_dc.cpp ${Heinz_Src} ${Heinz_Hdr} )
target_link_libraries( heinz_pcst_dc ${CommonLibs} )

add_executable( heinz_pcst_no_dc src/dimacs/heinz_pcst_no_dc.cpp ${Heinz_Src} ${Heinz_Hdr} )
target_link_libraries( heinz_pcst_no_dc ${CommonLibs} )

add_executable( heinz_pcst_no_pre src/dimacs/heinz_pcst_no_pre.cpp ${Heinz_Src} ${Heinz_Hdr} )
target_link_libraries( heinz_pcst_no_pre ${CommonLibs} )

add_executable( heinz_rpcst_mc src/dimacs/heinz_rpcst_mc.cpp ${Heinz_Monte_Carlo_Src} ${Heinz_Monte_Carlo_Hdr} ${Heinz_Hdr} )
target_link_libraries( heinz_rpcst_mc emon OGDF pthread )

add_executable( heinz_rpcst_no_dc src/dimacs/heinz_rpcst_no_dc.cpp ${Heinz_Src} ${Heinz_Hdr} )
target_link_libraries( heinz_rpcst_no_dc ${CommonLibs} )

add_executable( heinz_rpcst_no_pre src/dimacs/heinz_rpcst_no_pre.cpp ${Heinz_Src} ${Heinz_Hdr} )
target_link_libraries( heinz_rpcst_no_pre ${CommonLibs} )

add_executable( print EXCLUDE_FROM_ALL ${Heinz_Print_Graph_Src} ${Heinz_Hdr} ${CommonObjs} )
target_link_libraries( print emon OGDF pthread )

add_executable( check_mwcs_solution EXCLUDE_FROM_ALL src/dimacs/check_mwcs_solution.cpp src/utils.cpp )
target_link_libraries( check_mwcs_solution emon OGDF pthread )

add_executable( check_pcst_solution EXCLUDE_FROM_ALL src/dimacs/check_pcst_solution.cpp src/utils.cpp )
target_link_libraries( check_pcst_solution emon OGDF pthread )

enable_testing()
add_test( heinz_mwcs_no_dc ${PROJECT_SOURCE_DIR}/test/run.py ./heinz_mwcs_no_dc ./check_mwcs_solution ${PROJECT_SOURCE_DIR}/test/lymphoma.stp lymphoma.dimacs )
add_test( heinz_mwcs_dc ${PROJECT_SOURCE_DIR}/test/run.py ./heinz_mwcs_dc ./check_mwcs_solution ${PROJECT_SOURCE_DIR}/test/lymphoma.stp lymphoma.dimacs )
add_test( heinz_mwcs_no_pre ${PROJECT_SOURCE_DIR}/test/run.py ./heinz_mwcs_no_pre ./check_mwcs_solution ${PROJECT_SOURCE_DIR}/test/lymphoma.stp lymphoma.dimacs )
add_test( heinz_mwcs_mc ${PROJECT_SOURCE_DIR}/test/run.py ./heinz_mwcs_mc ./check_mwcs_solution ${PROJECT_SOURCE_DIR}/test/lymphoma.stp lymphoma.dimacs )
add_test( heinz_pcst_no_dc ${PROJECT_SOURCE_DIR}/test/run.py ./heinz_pcst_no_dc ./check_pcst_solution ${PROJECT_SOURCE_DIR}/test/i640-003.stp i640-003.dimacs )
add_test( heinz_pcst_dc ${PROJECT_SOURCE_DIR}/test/run.py ./heinz_pcst_dc ./check_pcst_solution ${PROJECT_SOURCE_DIR}/test/i640-003.stp i640-003.dimacs )
add_test( heinz_pcst_no_pre ${PROJECT_SOURCE_DIR}/test/run.py ./heinz_pcst_no_pre ./check_pcst_solution ${PROJECT_SOURCE_DIR}/test/i640-003.stp i640-003.dimacs )
add_test( heinz_pcst_mc ${PROJECT_SOURCE_DIR}/test/run.py ./heinz_pcst_mc ./check_pcst_solution ${PROJECT_SOURCE_DIR}/test/i640-003.stp i640-003.dimacs )
add_test( heinz_rpcst_no_dc ${PROJECT_SOURCE_DIR}/test/run.py ./heinz_rpcst_no_dc ./check_pcst_solution ${PROJECT_SOURCE_DIR}/test/i101M2.stp i101M2.dimacs )
add_test( heinz_rpcst_no_pre ${PROJECT_SOURCE_DIR}/test/run.py ./heinz_rpcst_no_pre ./check_pcst_solution ${PROJECT_SOURCE_DIR}/test/i101M2.stp i101M2.dimacs )
add_test( heinz_rpcst_mc ${PROJECT_SOURCE_DIR}/test/run.py ./heinz_rpcst_mc ./check_pcst_solution ${PROJECT_SOURCE_DIR}/test/i101M2.stp i101M2.dimacs )
add_custom_target( check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS 
  check_mwcs_solution 
  check_pcst_solution 
  heinz_pcst_dc 
  heinz_pcst_mc 
  heinz_pcst_no_dc 
  heinz_pcst_no_pre 
  heinz_rpcst_mc 
  heinz_rpcst_no_dc 
  heinz_rpcst_no_pre 
  heinz_mwcs_dc 
  heinz_mwcs_mc 
  heinz_mwcs_no_dc 
  heinz_mwcs_no_pre 
)
