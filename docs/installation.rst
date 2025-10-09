Installation
============

On All platforms (Linux/MacOS/Windows) platforms
----------------

The easiest way to install *simcoon* is to create a *conda* environnement: You can utilize the Anaconda GUI or type:
(for the installation of an environment called "simcoon")

The simplest way to install simcoon is directly with conda:

.. code-block:: none

    conda create --name simcoon

To activate the environment: 

.. code-block:: none

    conda activate simcoon

Now you can install *simcoon* using 

.. code-block:: none

    conda install -c conda-forge -c set3mah simcoon

That's it! You can start utilizing *simcoon*

Developper installation:
-------------------------

Prerequisites (using conda)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

We still recommand to use a specific environnement :

.. code-block:: none
    conda create --name simcoon_build

To activate the environment: 

.. code-block:: none
    conda activate simcoon_build

Install required dependencies:

.. code-block:: none
    # Compilers and build tools
    conda install -c conda-forge cxx-compiler fortran-compiler cmake ninja uv

    # Libraries
    conda install -c conda-forge armadillo boost pybind11 numpy gtest carma

    # Python testing and setuptools
    pip install pytest setuptools


For x86 architectures, you may also need MKL:
.. code-block:: none
    conda install -c conda-forge mkl


Prerequisites (without conda)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Install required dependencies using your system's package manager.

- On Debian/Ubuntu:

.. code-block:: none
    sudo apt-get install libarmadillo-dev libboost-all-dev libgtest-dev ninja-build

- On macOS with Homebrew:

.. code-block:: none
    brew install armadillo boost googletest

- On Windows with vcpkg:

.. code-block:: none
    #using powershell
    vcpkg install armadillo boost-config boost-dll gtest

Installation of Simcoon
~~~~~~~~~~~~~~~~~~~~~~~~

Next, download the Simcoon sources in the github repository of Simcoon_
.. _Simcoon : https://github.com/3MAH/simcoon or clone or download the repository (recommended):

.. code-block:: none
    git clone https://github.com/3MAH/simcoon.git
    cd simcoon

The last step is to run the installation script:

.. code-block:: none
    sh Install.sh

or:

**Linux/macOS:**

.. code-block:: none
    # Configure
    cmake -S . -B build -G Ninja -D CMAKE_BUILD_TYPE=Release

    # Build
    cmake --build build

    # Install Python package
    pip install ./build/python-package

**Windows:**

.. code-block:: none
    # using powershell
    # Configure
    cmake -S . -B build

    # Build
    cmake --build build --config Release

    # Install Python package
    pip install ./build/python-package

Run tests (All platforms)
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: none
    ctest --test-dir build --output-on-failure

A build folder will be automatically created in the Simcoon folder. At some point you can decide wether you will install or not the Simcoon library. Make sure you have carefully added thje path to your anaconda environnement.
Once the installation is done, the executables can be found in the build/bin folder. The use of python wrappers to those executables are however now easier to handle.

Here are some additional information about the prerequisites and the link to get them and their documentation:


- Boost_ (at least 1.63), including Boost Python
.. _Boost : https://www.boost.org
- Armadillo_ 
.. _Armadillo : http://arma.sourceforge.net

.. image:: _static/boost_logo.png
.. image:: _static/Armadillo_logo.png

Note that FTensor_ .. _FTensor : https://bitbucket.org/wlandry/ftensor
is also utilized by Simcoon but it is integrated to facilitate the installation.


