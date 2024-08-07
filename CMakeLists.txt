cmake_minimum_required(VERSION 3.17)
project(actsTest)

# soft installation path
set(SW $ENV{HOME}/mpd/sw/slc7_x86-64)
set(ACTS_SRC ${SW}/ACTS/latest)
set(Boost_DIR ${SW}/Boost/latest/lib/cmake/Boost-1.83.0)
set(Acts_DIR ${SW}/ACTS/latest/lib/cmake/Acts)
set(DD4hep_DIR ${SW}/DD4hep/latest/cmake)
set(Geant4_DIR ${SW}/GEANT4/latest/lib/cmake/Geant4)
set(ROOT_DIR ${SW}/ROOT/latest/cmake)
set(nlohmann_json_DIR ${SW}/nlohmann_json/latest/share/cmake/nlohmann_json)
set(LCIO_DIR ${SW}/LCIO/latest/lib/cmake/LCIO)
set(SIO_DIR ${SW}/LCIO/latest/lib/cmake/SIO)
set(EDM4HEP_DIR ${SW}/EDM4hep/latest/lib/cmake/EDM4HEP)
set(podio_DIR ${SW}/podio/latest/lib/cmake/podio)

find_package(Eigen3 REQUIRED HINTS "${SW}/Eigen3/latest")
find_package(Acts REQUIRED)

include_directories(
        $ENV{PWD}
        ${SW}/TBB/latest/include
        ${SW}/nlohmann_json/latest/include
        ${SW}/ROOT/latest/include
)
link_directories(
        ${SW}/TBB/latest/lib
        ${SW}/ACTS/latest/lib
        ${SW}/ROOT/latest/lib
        )

# main exec
add_executable(${PROJECT_NAME} main.cpp tracker.h tracker_config.h)

target_link_libraries(${PROJECT_NAME} PUBLIC
        ActsCore
        ActsExamplesFramework
        ActsExamplesGenerators
        ActsExamplesDetectorTelescope
        ActsExamplesFatras
        ActsPluginFpeMonitoring
        ActsExamplesMagneticField
        ActsExamplesDigitization
        ActsExamplesIoJson
        ActsExamplesIoRoot
        ActsExamplesIoCsv
        ActsExamplesTrackFinding
        ActsExamplesTruthTracking
        tbb
        Core
        Tree
)

# copy run macros
set(SCRIPTS
        digi-smearing-config.json
        )

foreach (_script ${SCRIPTS})
    configure_file(${_script}
            ${PROJECT_BINARY_DIR}/${_script}
            COPYONLY)
endforeach ()
