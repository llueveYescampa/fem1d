# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
get_filename_component(CurrDir ${CMAKE_CURRENT_SOURCE_DIR} NAME)

#Create a variable called linalg_SOURCES containing all .f90 files:

set(main_Src
    ${myPROJECT_SOURCE_DIR}/${CurrDir}/fem1d.f90
    PARENT_SCOPE
    )

set(fem1d_Inc
    ${myPROJECT_SOURCE_DIR}/fem1d/include
    PARENT_SCOPE
    )
