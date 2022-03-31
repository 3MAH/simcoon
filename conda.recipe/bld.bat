:: Armadillo
git clone https://gitlab.com/conradsnicta/armadillo-code.git armadillo
cd armadillo
cmake -S . -B build -G "Visual Studio 17 2022"
cmake --build build --config RELEASE
if errorlevel 1 exit 1

cd %SRC_DIR%
cmake -S . -B build ^
      -G"Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX:PATH=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -Wno-dev
cmake --build build --target simcoon --config Release
cmake --install build
if errorlevel 1 exit 1

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

cd %SRC_DIR%\simcoon-python-builder
cmake -S . -B build ^
      -G "Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DPython3_ROOT_DIR=%PREFIX% ^
      -Wno-dev
cat build\CMakeFiles\CMakeOutput.log
cmake --build build --target ALL_BUILD --config Release
cmake --install build
if errorlevel 1 exit 1

cd %SP_DIR%
mkdir simcoon
cd simcoon
type NUL > __init__.py
xcopy /s /i %SRC_DIR%\simcoon-python-builder\build\lib\Release\simmit.pyd %SP_DIR%\simcoon
if errorlevel 1 exit 1

xcopy /s /i %SRC_DIR%\armadillo\examples\lib_win64\libopenblas.dll %PREFIX%\Library\bin
xcopy /s /i %SRC_DIR%\build\lib\Release\simcoon.dll %PREFIX%\Library\bin
xcopy /s /i %SRC_DIR%\arma2numpy-builder\build\bin\Release\arma2numpy.dll %PREFIX%\Library\bin
if errorlevel 1 exit 1
