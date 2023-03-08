:: Armadillo
git clone https://gitlab.com/conradsnicta/armadillo-code.git armadillo
cd armadillo
copy /y %SRC_DIR%\armadillo\examples\lib_win64\libopenblas.dll %PREFIX%\Library\bin\
copy /y %SRC_DIR%\armadillo\examples\lib_win64\libopenblas.lib %PREFIX%\Library\lib\
cmake -S . -B build -G "Visual Studio 17 2022" -DBLAS_LIBRARY:FILEPATH="%PREFIX%\Library\lib\libopenblas.lib" -DLAPACK_LIBRARY:FILEPATH="%PREFIX%\Library\lib\libopenblas.lib"
cmake --build build --config RELEASE
cmake --install build --prefix %PREFIX%\Library
if errorlevel 1 exit 1

:: Simcoon
cd %SRC_DIR%
cmake -S . -B build ^
      -G"Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX:PATH=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -Wno-dev
cmake --build build --target simcoon --config Release
cmake --install build
if errorlevel 1 exit 1

:: Arma2numpy
cd arma2numpy-builder
cmake -S . -B build ^
      -G "Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DPython3_ROOT_DIR=%PREFIX% ^
      -Wno-dev 
cmake --build build --target ALL_BUILD --config Release
cmake --install build
if errorlevel 1 exit 1

:: Python binding
cd %SRC_DIR%\simcoon-python-builder
cmake -S . -B build ^
      -G "Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DPython3_ROOT_DIR=%PREFIX% ^
      -Wno-dev
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

REM xcopy /s /i %SRC_DIR%\armadillo\examples\lib_win64\libopenblas.dll %PREFIX%\Library\bin
xcopy /s /i %SRC_DIR%\build\lib\Release\simcoon.dll %PREFIX%\Library\bin
xcopy /s /i %SRC_DIR%\arma2numpy-builder\build\bin\Release\arma2numpy.dll %PREFIX%\Library\bin
if errorlevel 1 exit 1
