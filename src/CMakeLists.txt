# Include the directory itself as a path to include directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#Create a variable called linalg_SOURCES containing all .cpp files:

add_subdirectory(fem1d)
add_subdirectory(assemble)
add_subdirectory(boundary)
add_subdirectory(coeffcnt)
add_subdirectory(constrnt)
add_subdirectory(echodata)
add_subdirectory(egnsolvr)
add_subdirectory(eqnsolvr)
add_subdirectory(jacobi)
add_subdirectory(mesh1d)
add_subdirectory(matrxmlt)
add_subdirectory(postproc)
add_subdirectory(reaction)
add_subdirectory(shape1d)
add_subdirectory(timforce)
add_subdirectory(timstres)
add_subdirectory(transfrm)


set(fem1d_Src
    ${main_Src}
    ${assemble_Src}
    ${boundary_Src}
    ${coeffcnt_Src}
    ${constrnt_Src}
    ${echodata_Src}
    ${egnsolvr_Src}
    ${eqnsolvr_Src}
    ${jacobi_Src}
    ${mesh1d_Src}
    ${matrxmlt_Src}
    ${postproc_Src}
    ${reaction_Src}
    ${shape1d_Src}
    ${timforce_Src}
    ${timstres_Src}
    ${transfrm_Src}
    PARENT_SCOPE
    )

set(fem1d_Inc
    ${myPROJECT_SOURCE_DIR}/include
    ${fem1d_Inc}
    PARENT_SCOPE
    )
