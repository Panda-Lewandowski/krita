SET(PREFIX_ext_pyqt "${EXTPREFIX}" )
if (UNIX)
    ExternalProject_Add( ext_pyqt
        DOWNLOAD_DIR ${EXTERNALS_DOWNLOAD_DIR}
        URL http://files.kde.org/krita/build/dependencies/PyQt5_gpl-5.6.tar.gz
        URL_MD5 dbfc885c0548e024ba5260c4f44e0481
        
        CONFIGURE_COMMAND ${PREFIX_ext_pyqt}/bin/python3 <SOURCE_DIR>/configure.py --confirm-license
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        
        BUILD_IN_SOURCE 1

        UPDATE_COMMAND ""
    )
elseif(MINGW)
    list(APPEND PYQT_conf --confirm-license --target-py-version 3.6)
    list(APPEND PYQT_conf --bindir ${PREFIX_ext_pyqt}/bin --qt ${PREFIX_ext_pyqt})
    list(APPEND PYQT_conf --sip ${PREFIX_ext_pyqt}/bin/sip.exe --sip-incdir ${PREFIX_ext_pyqt}/include)
    list(APPEND PYQT_conf --spec win32-g++ --verbose --sipdir ${PREFIX_ext_pyqt}/share/sip)
    list(APPEND PYQT_conf --destdir ${PREFIX_ext_pyqt}/lib/site-packages --no-qml-plugin --no-python-dbus)
    list(APPEND PYQT_conf --no-qsci-api --no-tools --disable QtSql --disable QtTest --disable QtWinExtras)
    ExternalProject_Add( ext_pyqt
        DOWNLOAD_DIR ${EXTERNALS_DOWNLOAD_DIR}
        URL https://sourceforge.net/projects/pyqt/files/PyQt5/PyQt-5.9/PyQt5_gpl-5.9.zip
        URL_MD5 d978884753df265896eda436d8f4e07b
        PATCH_COMMAND  ${PATCH_COMMAND} -p1 -i ${CMAKE_CURRENT_SOURCE_DIR}/pyqt-configure-without-qt-ssl.patch

        CONFIGURE_COMMAND python.exe <SOURCE_DIR>/configure.py ${PYQT_conf}
        BUILD_COMMAND mingw32-make CXXFLAGS=-D_hypot=hypot
        INSTALL_COMMAND mingw32-make install
        
        BUILD_IN_SOURCE 1

        UPDATE_COMMAND ""
    )

endif()

