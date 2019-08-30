#cuba

add_library(cuba SHARED 
    src/vegas/Vegas.c	
    src/common/Data.c)

add_library(cuba::cuba ALIAS cuba)

target_include_directories(cuba PUBLIC include)
target_include_directories(cuba PRIVATE src/common)
