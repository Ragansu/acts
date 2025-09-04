add_executable(ActsAnalysisScoreMonitorPlots ScoreMonitorPlots.cpp)
target_link_libraries(
    ActsAnalysisScoreMonitorPlots
    ROOT::Core
    ROOT::Hist
    ROOT::Tree
    ROOT::TreePlayer
    Boost::program_options
)

install(
    TARGETS ActsAnalysisScoreMonitorPlots
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
