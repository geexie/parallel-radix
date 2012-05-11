# find intel compiler on win platform
# INTEL_COMPILER_ARCH
# INTEL_COMPILER_EXE
# INTEL_INCLUDE_DIR
# INTEL_LIBRARY_PATH

macro(find_intel_compiler)

    set(__versions 12 11 10)
    foreach(__each ${__versions})
        set(__icl_path $ENV{ICPP_COMPILER${__each}})
        string(REGEX REPLACE "\\\\" "/" __intel_cxx "${__icl_path}")
        if(NOT "${__intel_cxx}" STREQUAL "")
            set(INTEL_COMPILER_DIR "${__intel_cxx}")
            break()
        endif()
    endforeach()

    if(INTEL_COMPILER_DIR)
        # TODO: select arch
        set(INTEL_COMPILER_ARCH "ia32")
        set(INTEL_COMPILER_EXE "${INTEL_COMPILER_DIR}/bin/${INTEL_COMPILER_ARCH}/icl.exe")
        set(INTEL_INCLUDE_DIR "${INTEL_COMPILER_DIR}compiler/include")
        set(INTEL_LIBRARY_PATH "${INTEL_COMPILER_DIR}compiler/lib/${INTEL_COMPILER_ARCH}")
        message(STATUS "Intel compiler founded in path      ${INTEL_COMPILER_DIR}")
        message(STATUS "Intel compiler architecture         ${INTEL_COMPILER_ARCH}")
        message(STATUS "Intel compiler-stecific includes    ${INTEL_INCLUDE_DIR}")
        message(STATUS "Intel compiler-specific libraries   ${INTEL_LIBRARY_PATH}")
    elseif()
        message("${INTEL_COMPILER_DIR}!!!!!")
        set(INTEL_COMPILER_DIR "icl-NOT_FOUNDED")
    endif()
endmacro()
