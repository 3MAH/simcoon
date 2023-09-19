:: Simcoon
cd %SRC_DIR%
cmake -G Ninja ^
      -DCMAKE_INSTALL_PREFIX:PATH=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -Wno-dev ^
      -S . -B build

cmake --build build --target simcoon --config Release
cmake --install build
if errorlevel 1 exit 1

:: Python binding
cd %SRC_DIR%\simcoon-python-builder
cmake -G Ninja ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DPython3_ROOT_DIR=%PREFIX% ^
      -Wno-dev ^
      -S . -B build
cat build\CMakeFiles\CMakeOutput.log
dir build\CMakeFiles
cat build\CMakeFiles\CMakeError.log
cmake --build build --target ALL_BUILD --config Release
cmake --install build
if errorlevel 1 exit 1

:: Install simcoon python 
cd %SP_DIR%
mkdir simcoon
cd simcoon
type NUL > __init__.py
xcopy /s /i %SRC_DIR%\simcoon-python-builder\build\lib\Release\simmit.pyd %SP_DIR%\simcoon
if errorlevel 1 exit 1

xcopy /s /i %SRC_DIR%\build\lib\Release\simcoon.dll %PREFIX%\Library\bin
if errorlevel 1 exit 1
