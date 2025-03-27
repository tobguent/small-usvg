
macro(enable_ide_support_generic TARGET_NAME SOURCES_LIST SRC_DIRECTORY HEADER_DIRECTORY FOLDER_STRING)
    # Find header files
    file(GLOB HEADERS CONFIGURE_DEPENDS ${HEADER_DIRECTORY}/*.hpp)
    # include header in target
    target_sources(${TARGET_NAME} PRIVATE ${HEADERS})
    # definition for source group
    source_group(TREE "${HEADER_DIRECTORY}" PREFIX "Header Files" FILES ${HEADERS})
    source_group(TREE "${SRC_DIRECTORY}" PREFIX "Source Files" FILES ${SOURCES_LIST})
    # place target into folder
    set_target_properties(${TARGET_NAME} PROPERTIES FOLDER ${FOLDER_STRING})
endmacro()

macro(enable_ide_support TARGET_NAME SOURCES_LIST DIRECTORY FOLDER_STRING)
    enable_ide_support_generic(${TARGET_NAME} "${SOURCES_LIST}" ${DIRECTORY}/src ${DIRECTORY}/include/${TARGET_NAME} ${FOLDER_STRING})
endmacro()

macro(enable_ide_support_vislab TARGET_NAME SOURCES_LIST DIRECTORY)
    enable_ide_support(${TARGET_NAME} "${SOURCES_LIST}" ${DIRECTORY} "vislab")
endmacro()
