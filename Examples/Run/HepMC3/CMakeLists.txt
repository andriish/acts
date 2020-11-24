add_executable(
  ActsExampleHepMC3
  HepMC3Example.cpp)
target_link_libraries(
  ActsExampleHepMC3
  PRIVATE
    ActsCore ActsExamplesCommon
    ActsExamplesFramework ActsExamplesIoHepMC3
    Boost::program_options)

add_executable(
  ActsExampleMean
  HepMC3MeanExample.cpp)
target_link_libraries(
  ActsExampleMean
  PRIVATE
    ActsCore ActsExamplesCommon ActsExamplesDetectorDD4hep
    ActsExamplesFramework ActsExamplesIoHepMC3 ActsExamplesHepMC
    Boost::program_options)
    
install(
  TARGETS ActsExampleHepMC3 ActsExampleMean
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  
