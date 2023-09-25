:: Simcoon
cd %SRC_DIR%
cmake -G"Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX:PATH=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -Wno-dev ^
      -S . -B build

cmake --build build --target simcoon --config Release
cmake --install build
if errorlevel 1 exit 1

:: Arma2numpy binding
cd %SRC_DIR%\arma2numpy-builder
cmake -G"Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DPython3_ROOT_DIR=%PREFIX% ^
      -Wno-dev ^
      -S . -B build
      
cmake --build build --target arma2numpy --config Release
cmake --install build
if errorlevel 1 exit 1

:: Python binding
cd %SRC_DIR%\simcoon-python-builder
cmake -G"Visual Studio 17 2022" ^
      -DCMAKE_INSTALL_PREFIX=%PREFIX% ^
      -DCMAKE_BUILD_TYPE=Release ^
      -DPython3_ROOT_DIR=%PREFIX% ^
      -Wno-dev ^
      -S . -B build

cmake --build build --target simmit --config Release
cmake --install build
if errorlevel 1 exit 1

:: Install simcoon python 
cd %SRC_DIR%\python-setup
xcopy /s /i %SRC_DIR%\simcoon-python-builder\build\lib\Release\simmit.pyd %SRC_DIR%\python-setup\simcoon
if errorlevel 1 exit 1

%PYTHON% -m pip install .

:: xcopy /s /i %SRC_DIR%\build\lib\Release\simcoon.dll %PREFIX%\Library\bin
:: if errorlevel 1 exit 1
