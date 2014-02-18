cmake_minimum_required( VERSION 2.8 )

project( Heinz )

set( BK_MaxFlow_Hdr
  src/externlib/maxflow-v3.01/graph.h
  src/externlib/maxflow-v3.01/block.h
)

set( Heinz_Hdr
  src/parser/parser.h
  src/parser/mwcsparser.h
  src/parser/mwcsstpparser.h
  src/parser/identityparser.h
  src/preprocessing/mwcspreprocessrootrule.h
  src/preprocessing/mwcspreprocessrule.h
  src/preprocessing/mwcspreprocessrulebase.h
  src/preprocessing/mwcspreprocessrulenegdeg01.h
  src/preprocessing/mwcspreprocessrulenegedge.h
  src/preprocessing/mwcspreprocessruleneghub.h
  src/preprocessing/mwcspreprocessruleposdeg01.h
  src/preprocessing/mwcspreprocessruleposedge.h
  src/verbose.h
  src/mwcs.h
  src/mwcsanalyze.h
  src/mwcsenumerate.h
  src/mwcsenumeratecomp.h
  src/mwcsenumerateroot.h
  src/solver/mwcssolver.h
  src/solver/mwcstreeheuristicsolver.h
  src/solver/mwcstreesolver.h
  src/solver/cplex_cut/nodecutbk.h
  src/solver/cplex_cut/nodecutrootedbk.h
  src/solver/cplex_cut/nodecutrootedbkcallback.h
  src/solver/cplex_cut/nodecutunrootedbk.h
  src/solver/cplex_cut/nodecutunrootedbkcallback.h
  src/solver/cplex_cut/bk_alg.h
  src/solver/mwcscplexsolver.h
  src/solver/mwcssizecutsolver.h
  src/solver/mwcssizetreememsolver.h
  src/solver/mwcssizetreesolver.h
  src/solver/cplex_heuristic/mwcssizecutsolverheuristic.h
  src/solver/cplex_heuristic/mwcscutsolverheuristic.h
  src/solver/mwcscutsolver.h
  src/mwcsgraph.h
  src/mwcspreprocessedgraph.h
  src/utils.h
  src/mwcsgraphparser.h
)

set( Heinz_Src
  src/mwcs.cpp
  src/utils.cpp
  src/verbose.cpp
)

set( Heinz_Monte_Carlo_Src
  src/mwcs-mc.cpp
  src/utils.cpp
  src/verbose.cpp
)

set( Heinz_Size_Src
  src/mwcs-s.cpp
  src/utils.cpp
  src/verbose.cpp
)

set( Heinz_Print_Graph_Src
  src/printgraph.cpp
  src/utils.cpp
  src/verbose.cpp
)

set( GitCommand1 "symbolic-ref HEAD 2> /dev/null | cut -b 12-" )
execute_process( COMMAND git symbolic-ref HEAD  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                 COMMAND tail -c -7 COMMAND tr -d "\\n" OUTPUT_VARIABLE GitBranch )
set( GitCommand2 "log --pretty=format:\\\"%h\\\" -1" )
execute_process( COMMAND git log --pretty=format:%h WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                 COMMAND head -n 1  COMMAND tr -d "\\n" OUTPUT_VARIABLE GitHashtag )
set( HEINZ_VERSION \"${GitBranch}-${GitHashtag}\" )

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
	set( CPLEX_AUTODETECTION_PREFIXES "$ENV{HOME}/ILOG;/opt/ibm/ILOG" CACHE STRING "CPLEX automatic lookup directories" )
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

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -g -ggdb -O0 -DDEBUG -DVERBOSE")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_PROFILE "-Wall -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -O3 -g -ggdb -DNDEBUG")

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings " )
add_definitions( -DIL_STD )
configure_file( src/config.h.in ${PROJECT_SOURCE_DIR}/src/config.h )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
endif()

include_directories( ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} "${LIBLEMON_ROOT}/include" "src" "src/externlib" )
link_directories( ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} "${LIBLEMON_ROOT}/lib" )

set( CommonLibs
	m
	ilocplex
	cplex
	concert
	pthread
	emon
)

add_executable( heinz ${Heinz_Src} ${Heinz_Hdr} ${CommonObjs} )
target_link_libraries( heinz ${CommonLibs} )

add_executable( heinz-mc ${Heinz_Monte_Carlo_Src} ${Heinz_Hdr} ${CommonObjs} )
target_link_libraries( heinz-mc emon )

add_executable( print ${Heinz_Print_Graph_Src} ${Heinz_Hdr} ${CommonObjs} )
target_link_libraries( print emon )
