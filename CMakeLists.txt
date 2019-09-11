# cuba
project(cuba)

add_library(cuba SHARED
    src/vegas/Vegas.c
    src/common/Data.c)

add_library(cuba::cuba ALIAS cuba)

target_include_directories(cuba PUBLIC include)
target_include_directories(cuba PRIVATE src/common)

# Export interface symbols and hide internal implementation
include(GenerateExportHeader)
generate_export_header(cuba)
target_include_directories(cuba PUBLIC ${PROJECT_BINARY_DIR})
set_target_properties(cuba PROPERTIES CXX_VISIBILITY_PRESET hidden)
