add_library(
  ActsExamplesIoRoot SHARED
  src/RootDigitizationWriter.cpp
  src/RootMaterialDecorator.cpp
  src/RootMaterialWriter.cpp
  src/RootMaterialTrackReader.cpp
  src/RootMaterialTrackWriter.cpp
  src/RootPlanarClusterWriter.cpp
  src/RootParticleWriter.cpp
  src/RootPropagationStepsWriter.cpp
  src/RootSimHitWriter.cpp
  src/RootTrackParameterWriter.cpp
  src/RootTrajectoryStatesWriter.cpp
  src/RootTrajectoryParametersWriter.cpp
)
target_include_directories(
  ActsExamplesIoRoot
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(
  ActsExamplesIoRoot
  PUBLIC
    ActsCore ActsPluginDigitization ActsPluginIdentification ActsExamplesDigitization
    ActsExamplesFramework ActsExamplesPropagation ActsExamplesTruthTracking
    Threads::Threads
  PRIVATE ROOT::Core ROOT::Hist ROOT::Tree)

install(
  TARGETS ActsExamplesIoRoot
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
