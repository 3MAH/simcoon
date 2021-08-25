Installation
========

On MacOS/Linux platforms
----------------

The easiest way to install simcoon is to create a *conda* environnement: You can utilize the Anaconda GUI or type:
(for the installation of an environment called "scientific")

.. code-block:: none

    conda create --name scientific

To activate the environment: 

.. code-block:: none

    conda activate scientific

The next step is to install the required packages:

.. code-block:: none

    conda install -c conda-forge armadillo 
    conda install -c conda-forge boost 
    conda install -c conda-forge cgal 
    conda install -c conda-forge numpy

Next, download the Simcoon sources in the github repository of Simcoon_
.. _Simcoon : https://github.com/3MAH/simcoon
Unzip the content in a folder and modify the Install.sh source file to look after you conda environnement path:

anacondaloc=/path/to/anaconda/anaconda3/envs/scientific

The last step is to run the installation script:

.. code-block:: none

sh Install.sh

A build folder will be automatically created in the Simcoon folder. At some point you can decide wether you will install or not the Simcoon library. Make sure you have carefully added thje path to your anaconda environnement.
Once the installation is done, the executables can be found in the build/bin folder. The use of python wrappers to those executables are however now easier to handle.

Note: You shall make sure that you have CMake installed

.. image:: _static/CMake.png

If not installed, for Ubuntu and debian-based systems:

.. code-block:: none

    sudo apt-get install cmake 

And for Mac OS user, you can use brew:

.. code-block:: none

   brew update
   brew install cmake

If you do no want to install Simcoon using a conda environnement, the following dependencies are required to install simcoon: 

- Boost_ (at least 1.63), including Boost Python
.. _Boost : https://www.boost.org
- Armadillo_ 
.. _Armadillo : http://arma.sourceforge.net
- CGAL_
.. _CGAL : https://www.cgal.org

.. image:: _static/boost_logo.png
.. image:: _static/Armadillo_logo.png
.. image:: _static/CGAL_logo.png

Note that FTensor_ is also utilized by Simcoon but it is integrated to facilitate the installation. You can get the sources and docs here
.. _FTensor : https://bitbucket.org/wlandry/ftensor

Make sure that you have access to the folder selected for the installation with Cmake (by default /usr/local on most Unix-based systems).

On Windows platforms
----------------

The following procedure has been tested on Windows 10 64 bits:

The first thing is to download the last version of Visual Studio.

1. Install Anaconda using windows 64bits installer
2. download and execute the CMake Win64 Installer (https://cmake.org/download/). Make sure to set the PATH variable during installation
3. download and install Visual studio (tested with VS 2019). You can get if here_
_here : https://visualstudio.microsoft.com/downloads/
4. Download simcoon from Github : https://github.com/3MAH/simcoon/

*To be completed*


