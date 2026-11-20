# -*- cmake -*-
CMAKE_MINIMUM_REQUIRED(VERSION 3.12)
PROJECT(HighLANDPD)

# ---- create folders ----------
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")
file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/include")
file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/dict")

#-------------------------------------------------------
# ----- Set cmake options ------------------------------
#-------------------------------------------------------

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-g -pipe -Wall -Wnon-virtual-dtor -fpic -Woverloaded-virtual -std=c++17")

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}")

message("CMAKE source dir: ${CMAKE_SOURCE_DIR}")
message("CMAKE binary dir: ${CMAKE_BINARY_DIR}")
message("CMAKE compiler:   ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_FLAGS:  ${CMAKE_CXX_FLAGS}")
message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")

#-------------------------------------------------------
#------- External packages -----------------------------
#-------------------------------------------------------

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

# Look for the required libraries

#------ ROOT -------------------
# Add Minuit for CoherentFit Algorithm
find_package(ROOT REQUIRED COMPONENTS Geom EG Minuit MathMore)
#find_package(ROOT REQUIRED COMPONENTS Geom EG)
include(${ROOT_USE_FILE})

message("ROOT_USE_FILE:            ${ROOT_USE_FILE}")
message("ROOT include directories: ${ROOT_INCLUDE_DIR}")
message("ROOT libraries:           ${ROOT_LIBRARIES}")

#------ HighLAND -------------------
find_package(HIGHLAND REQUIRED)

message("HighLAND_USE_FILE:            ${HIGHLAND_USE_FILE}")
message("HighLAND include directories: ${HIGHLAND_INCLUDE_DIRS}")
message("HighLAND libraries:           ${HIGHLAND_LIBRARIES}")
include_directories(${HIGHLAND_INCLUDE_DIRS})
link_directories(${HIGHLAND_LIBRARY_DIR})

#-------------------------------------------------------
#----- include directories -----------------------------
#-------------------------------------------------------

file(GLOB highlandPD_inc ${CMAKE_SOURCE_DIR}/src/*/src
#                         ${CMAKE_SOURCE_DIR}/src/secondaryKaonAnalysis/src/coherentfit
                         ${CMAKE_BINARY_DIR}/dict)

include_directories(${highlandPD_inc})
message("Include directories: ${INCLUDE_DIRECTORIES}")

#-------------------------------------------------------
#------- ROOT Dictionaries -----------------------------
#-------------------------------------------------------

file(GLOB dict_src ${CMAKE_SOURCE_DIR}/src/*/src/*_LinkDef.h)
#                   ${CMAKE_SOURCE_DIR}/src/secondaryKaonAnalysis/src/coherentfit/*_LinkDef.h)

foreach(file0 ${dict_src})
  string(FIND ${file0} "/"          pos_ini REVERSE)
  string(SUBSTRING ${file0} 0      ${pos_ini}  dir)
  string(REPLACE "${dir}/"     "" file2 ${file0})
  string(REPLACE _LinkDef.h    "" file3 ${file2})  

# remove dict/ from name, new CMake version does not like '/' in names. To be solved.
#  root_generate_dictionary(${file3}_Dict NODEPHEADERS ${file3}.hxx LINKDEF ${file0}  MODULE ${file3} DEPENDENCIES ${ROOT_LIBRARIES})
  root_generate_dictionary(${file3}_Dict ${dir}/${file3}.hxx LINKDEF ${file0}  MODULE ${file3} DEPENDENCIES ${ROOT_LIBRARIES})

# remove dict/ from name, new CMake version does not like '/' in names. To be solved.
  list(APPEND dict_list ${dir}/${file3}.cxx ${file3}_Dict.cxx)
  list(APPEND dict_src  ${CMAKE_BINARY_DIR}/${file3}_Dict.cxx)

  INSTALL(FILES "${CMAKE_BINARY_DIR}/lib${file3}_rdict.pcm" DESTINATION ${CMAKE_SOURCE_DIR}/lib)
  INSTALL(FILES "${CMAKE_BINARY_DIR}/lib${file3}.rootmap"   DESTINATION ${CMAKE_SOURCE_DIR}/lib)
endforeach(file0)


#-------------------------------------------------------
#-------- the highlandPD library -------------------------
#-------------------------------------------------------

file(GLOB highlandPD_src  ${CMAKE_SOURCE_DIR}/src/*/src/*.cxx)
#                          ${CMAKE_SOURCE_DIR}/src/secondaryKaonAnalysis/src/coherentfit/*.cxx)

add_library(highlandPD SHARED ${dict_src} ${highlandPD_src} )
target_include_directories(highlandPD PRIVATE ${HIGHLAND_INCLUDE_DIRS} ${highlandPD_inc})

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

target_link_libraries(highlandPD ${ROOT_LIBRARIES} ${HIGHLAND_LIBRARIES})
INSTALL(TARGETS highlandPD DESTINATION ${CMAKE_SOURCE_DIR}/lib)

set(HIGHLAND_PD_LIBRARIES "-L${CMAKE_BINARY_DIR} -lhighlandPD")

#-------------------------------------------------------
#----- highlandPD includes -------------------------------
#-------------------------------------------------------

file(GLOB highlandPD_inc2 ${CMAKE_SOURCE_DIR}/src/*/src/*.hxx)
#                          ${CMAKE_SOURCE_DIR}/src/secondaryKaonAnalysis/src/coherentfit/*.hxx)

foreach(file0 ${highlandPD_inc2})
   INSTALL(FILES ${file0} DESTINATION ${CMAKE_SOURCE_DIR}/include)
endforeach(file0)

#-------------------------------------------------------
#-------------- executable -----------------------------
#-------------------------------------------------------

# Cada ejecutable se compila a partir de un único Run*.cxx y enlaza con
# libhighlandPD. Solo dependen de la librería, no entre ellos, así que
# pueden compilarse en paralelo.
function(highlandpd_add_executable name source)
  add_executable(${name} ${CMAKE_SOURCE_DIR}/${source})
  target_include_directories(${name} PRIVATE ${HIGHLAND_INCLUDE_DIRS} ${highlandPD_inc})
  target_link_libraries(${name} ${ROOT_LIBRARIES} ${HIGHLAND_LIBRARIES} ${HIGHLAND_PD_LIBRARIES})
  add_dependencies(${name} highlandPD)
  # Mantiene -rdynamic, que CMake < 3.4 añadía siempre (política CMP0065)
  set_target_properties(${name} PROPERTIES ENABLE_EXPORTS ON)
  INSTALL(TARGETS ${name} DESTINATION ${CMAKE_SOURCE_DIR}/bin)
endfunction()

highlandpd_add_executable(secondaryKaon.exe                 src/secondaryKaonAnalysis/app/RunSecondaryKaonAnalysis.cxx)
highlandpd_add_executable(secondaryProton.exe               src/secondaryKaonAnalysis/app/RunSecondaryProtonAnalysis.cxx)
highlandpd_add_executable(secondaryKaonProductionXS.exe     src/secondaryKaonAnalysis/app/RunSecondaryKaonProductionXSAnalysis.cxx)
highlandpd_add_executable(CreateMiniTreeKaon.exe            src/secondaryKaonAnalysis/app/RunCreateMiniTreeKaon.cxx)
highlandpd_add_executable(CreateMiniTreeEventCandidates.exe src/secondaryKaonAnalysis/app/RunCreateMiniTreeEventCandidates.cxx)
highlandpd_add_executable(CreateMiniTreeCandidates.exe      src/secondaryKaonAnalysis/app/RunCreateMiniTreeCandidates.cxx)
highlandpd_add_executable(pdExample.exe                     src/pdExampleAnalysis/app/RunpdExampleAnalysis.cxx)
highlandpd_add_executable(dEdxCalibration.exe               src/dEdxCalibrationAnalysis/app/RundEdxCalibration.cxx)
highlandpd_add_executable(dQdxYZCalibration.exe             src/dEdxCalibrationAnalysis/app/RundQdxYZCalibration.cxx)
highlandpd_add_executable(dQdxXCalibration.exe              src/dEdxCalibrationAnalysis/app/RundQdxXCalibration.cxx)
highlandpd_add_executable(CreateMichelRemovingMiniTree.exe  src/dEdxCalibrationAnalysis/app/RunCreateMichelRemovingMiniTree.cxx)
highlandpd_add_executable(brokenTracks.exe                  src/brokenTracksAnalysis/app/RunBrokenTracksAnalysis.cxx)
highlandpd_add_executable(beamWeights.exe                   src/beamWeightsAnalysis/app/RunBeamWeightsAnalysis.cxx)
highlandpd_add_executable(stoppingProton.exe                src/stoppingProtonAnalysis/app/RunStoppingProtonAnalysis.cxx)

# Desactivados:
#highlandpd_add_executable(CreateMiniTreePD.exe src/pdIO/app/RunCreateMiniTreePD.cxx)
#highlandpd_add_executable(CoherentFit.exe      src/secondaryKaonAnalysis/src/coherentfit/RunCoherentFit.cxx)

#-------------------------------------------------------
#------- custom targets --------------------------------
#-------------------------------------------------------

# target to create dictionaries before libraries
add_custom_target(highlandPD_dicts ALL SOURCES ${dict_list})

# target to change the library name in the rootmap files
add_custom_target(ch_rootmap COMMAND ${CMAKE_SOURCE_DIR}/scripts/chrootmap ${CMAKE_BINARY_DIR})

# target to install header files in the include folder
add_custom_target(include ALL SOURCES ${highlandPD_inc2} ${CMAKE_SOURCE_DIR})

#-------------------------------------------------------
#------- dependencies ----------------------------------
#-------------------------------------------------------

add_dependencies(highlandPD                        ch_rootmap)
add_dependencies(ch_rootmap                        highlandPD_dicts)
add_dependencies(highlandPD_dicts                  include)


#-------------------------------------------------------
#------- other -----------------------------------------
#-------------------------------------------------------

#---Configure and install various files neded later and for clients -----------------------------
include(HighlandPDConfiguration)
