#!/bin/bash
# Command line for building the CodeLite project files for the Release and Debug cpp0x versions of Gaalet
# for evaluating the library on Linux x64 architectures (a different build script will be needed on other platforms)


#pushd ./cmake/cpp0x/build-Release
#BUILD_ROOT=`pwd`
#cmake -G "CodeLite - Ninja" \
#-DOPENTHREADS_LIBRARY_DEBUG:FILEPATH="/usr/local/lib64/libOpenThreads.so" \
#-DOSGVIEWER_LIBRARY_DEBUG:FILEPATH="/usr/local/lib64/libosgViewer.so" \
#-DOSGVIEWER_LIBRARY:FILEPATH="/usr/local/lib64/libosgViewer.so" \
#-DCMAKE_INSTALL_PREFIX:PATH="/usr/local" \
#-DOSG_LIBRARY:FILEPATH="/usr/local/lib64/libosg.so;/usr/local/lib64/libosgGA.so" \
#-DOSG_LIBRARY_DEBUG:FILEPATH="/usr/local/lib64/libosg.so;/usr/local/lib64/libosgGA.so" \
#-DOPENTHREADS_LIBRARY:FILEPATH="/usr/local/lib64/libOpenThreads.so" \
#-DOSGVIEWER_INCLUDE_DIR:PATH="/usr/local/include" \
#-DCMAKE_MODULE_PATH:PATH="$BUILD_ROOT/../../modules" \
#-DTBB_INCLUDE_DIR="/usr/include/tbb" \
#-DTBB_LIBRARY="/usr/lib/x86_64-linux-gnu" \
#-DCMAKE_BUILD_TYPE:STRING="Release" \
 #../../../cpp0x

#popd

pushd ./cmake/cpp0x/build-Debug
BUILD_ROOT=`pwd`
cmake -G "CodeLite - Ninja" \
-DOPENTHREADS_LIBRARY_DEBUG:FILEPATH="/usr/local/lib64/libOpenThreads.so" \
-DOSGVIEWER_LIBRARY_DEBUG:FILEPATH="/usr/local/lib64/libosgViewer.so" \
-DOSGVIEWER_LIBRARY:FILEPATH="/usr/local/lib64/libosgViewer.so" \
-DCMAKE_INSTALL_PREFIX:PATH="/usr/local" \
-DOSG_LIBRARY:FILEPATH="/usr/local/lib64/libosg.so;/usr/local/lib64/libosgGA.so" \
-DOSG_LIBRARY_DEBUG:FILEPATH="/usr/local/lib64/libosg.so;/usr/local/lib64/libosgGA.so" \
-DOPENTHREADS_LIBRARY:FILEPATH="/usr/local/lib64/libOpenThreads.so" \
-DOSGVIEWER_INCLUDE_DIR:PATH="/usr/local/include" \
-DCMAKE_MODULE_PATH:PATH="$BUILD_ROOT/../../modules" \
-DTBB_INCLUDE_DIR="/usr/include/tbb" \
-DTBB_LIBRARY="/usr/lib/x86_64-linux-gnu" \
-DCMAKE_BUILD_TYPE:STRING="Debug" \
 ../../../cpp0x
popd
