Installation
========

On MacOS/Linux platforms
----------------

You need to have CMake installed

.. image:: _static/CMake.png

If not installed, for Ubuntu and debian-based systems:

.. code-block:: none

    sudo apt-get install cmake 

And for Mac OS user, you can use brew:

.. code-block:: none

   brew update
   brew install cmake

The following dependencies are required to install simcoon: 

- Boost_ (at least 1.63), including Boost Python
.. _Boost : https://www.boost.org
- Armadillo_ 
.. _Armadillo : http://arma.sourceforge.net
- CGAL_
.. _CGAL : https://www.cgal.org
- FTensor_
.. _FTensor : https://bitbucket.org/wlandry/ftensor

.. image:: _static/boost_logo.png
.. image:: _static/Armadillo_logo.png
.. image:: _static/CGAL_logo.png

Download simcoon from Github : https://github.com/simcoon/simcoon/

Copy the files in a folder. Using a terminal, navitate to such folder and run the installation script:

.. code-block:: none

    sh Install.sh
    
A build folder will be automatically created. At some point you can decide wether you will install or not the Simcoon library (by default in usr/local/)
You can set up, in the CMakeLists.txt file, the location where the library be installed. 
Once the installation is done, the executables can be found in the build/bin folder

Make sure that you have access to /usr/local folder for both MacOS and Linux operating systems.

On Windows platforms
----------------

The following procedure has been tested on Windows 10 64 bits:

The first thing is to download Visual Studio 2019.

1. Install Anaconda using windows 64bits installer
2. download and execute the CMake Win64 Installer (https://cmake.org/download/). Make sure to set the PATH variable during installation
3. download and install Visual studio (tested with VS 2019). You can get if here_
_here : https://visualstudio.microsoft.com/downloads/
4. Download simcoon from Github : https://github.com/simcoon/simcoon/


