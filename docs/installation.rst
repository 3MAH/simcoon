
Installation
============

All Platforms (Linux, macOS, Windows)
-------------------------------------

The recommended way to install *simcoon* is with *conda*. You can use the Anaconda GUI or run the following commands in your terminal:

.. code-block:: none

    conda create --name simcoon
    conda activate simcoon
    conda install -c conda-forge -c set3mah simcoon

*simcoon* is now ready to use.


Developer Installation
---------------------

Prerequisites (using conda)
~~~~~~~~~~~~~~~~~~~~~~~~~~

It is recommended to use a dedicated environment for development:

.. code-block:: none
    conda create --name simcoon_build
    conda activate simcoon_build

Install the required dependencies:

.. code-block:: none
    # Compilers and build tools
    conda install -c conda-forge cxx-compiler fortran-compiler cmake ninja uv
    # Libraries
    conda install -c conda-forge armadillo pybind11 numpy gtest carma
    # Python testing and setuptools
    pip install pytest setuptools

For x86 architectures, you may also need MKL:
.. code-block:: none
    conda install -c conda-forge mkl



Prerequisites (without conda)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Alternatively, you can install the dependencies using your system's package manager:

- **Debian/Ubuntu:**

    .. code-block:: none
            sudo apt-get install libarmadillo-dev libgtest-dev ninja-build

- **macOS (Homebrew):**

    .. code-block:: none
            brew install armadillo googletest

- **Windows (vcpkg, PowerShell):**

    .. code-block:: none
            vcpkg install armadillo gtest


Simcoon Installation
~~~~~~~~~~~~~~~~~~~~

Download the Simcoon source code from the GitHub repository:
.. _Simcoon : https://github.com/3MAH/simcoon

.. code-block:: none
    git clone https://github.com/3MAH/simcoon.git
    cd simcoon

To install, you can use the provided script:

.. code-block:: none
    sh Install.sh

Alternatively, you can build manually:

**Linux/macOS:**

.. code-block:: none
    cmake -S . -B build -G Ninja -D CMAKE_BUILD_TYPE=Release
    cmake --build build
    pip install ./build/python-package

**Windows:**

.. code-block:: none
    cmake -S . -B build
    cmake --build build --config Release
    pip install ./build/python-package


Running Tests (All Platforms)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none
    ctest --test-dir build --output-on-failure

The build folder is created automatically in the Simcoon directory. After installation, executables are located in `build/bin`. Python wrappers are available for easier usage.

Additional Information
~~~~~~~~~~~~~~~~~~~~~

- [Armadillo](http://arma.sourceforge.net)

.. image:: _static/Armadillo_logo.png

Note: [FTensor](https://bitbucket.org/wlandry/ftensor) is also used by Simcoon, but it is included in the source for easier installation.


