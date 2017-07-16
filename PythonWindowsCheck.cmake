# Check whether the found python is the same as ext_python
# HACK: Find python36.dll and compare equality. Probably not the best idea...
# TODO: Check the python version
if(EXISTS "${CMAKE_INSTALL_PREFIX}/python/python36.dll")
    message(STATUS "python36.dll is found in \"${CMAKE_INSTALL_PREFIX}/python/\".")
    file(SHA1 "${CMAKE_INSTALL_PREFIX}/python/python36.dll" _ext_python_dll_sha1)
    get_filename_component(_found_python_dir ${PYTHON_EXECUTABLE} DIRECTORY)
    file(SHA1 "${_found_python_dir}/python36.dll" _found_python_dll_sha1)
    if(NOT ${_ext_python_dll_sha1} STREQUAL ${_found_python_dll_sha1})
        message(FATAL_ERROR "Python versions mismatch.")
    endif()
else()
    message(STATUS "python36.dll is NOT found in \"${CMAKE_INSTALL_PREFIX}/python/\".")
endif()
