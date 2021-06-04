#!/bin/bash

helpFunction()
{
   echo ""
   echo "Usage: $0 -t test -n ncpus"
   echo "\t-t run the script without executing tests"
   echo "\t-n Number of cpus for the compilation"
   exit 1 # Exit script after printing help
}

test=1
ncpus=4
anacondaloc=/Users/ychemisky/opt/anaconda3/envs/test

while getopts "tn:" opt
do
   case $opt in
      t ) test=0 ;;
      n ) ncpus="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$test" ] || [ -z $ncpus ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct

echo "\n-----------------------------"
echo "Start of Simcoon compilation."
echo "-----------------------------\n"

blue=`tput setaf 4`
red=`tput setaf 1`
reset=`tput sgr0`

#Find the current directory
current_dir=$(pwd)

#Get the file of the bash file (possibly necessary)
bashFileName=${0##*/}

# Detect the platform (similar to $OSTYPE)
OS="`uname`"
case $OS in
'Linux')
OS='Linux'
alias finder='find'
alias nprocs='nproc'
;;
'Darwin')
OS='Mac'
alias finder='gfind'
alias nprocs='sysctl -n hw.ncpu'
;;
'WindowsNT')
OS='Windows'
echo 'No equivalent function known to be equivalent to "find or gfind" : Please search for one and replace this "echo" command in ${bashFileName} and replace it by an "alias"'
;;
*) ;;
esac

#Number of procs used to compil (+1 because if nproc=1 => 1/2=0)
#nproc_used=$(( ($(nprocs)+1)/2 ))

if [ ! -d "exec" ]
then
    mkdir ${current_dir}/exec
    echo "exec folder created.\n"
else
    echo "exec directory already exists."
fi

#Test if build exist and if it's necessary to erase it
if [ ! -d "build" ]
then
	mkdir ${current_dir}/build
	echo "Build folder created.\n"
else
	echo "Build directory already exists."
	
	while true; do
		read -p "Do you want to erase old compilation files (Recommended : No) ? " yn
		case $yn in
			[YyOo]* ) rm -r ${current_dir}/build/*; break;;
			[Nn]* ) break;;
			* ) echo "Please answer yes (y) or no (n).";;
		esac
	done
fi

#Ask for installation of the simcoon library
while true; do
	read -p "Do you want to install simcoon library (necessary to use libsimcoon.so and simmit) ? " yn
	case $yn in
		[YyOo]* ) Install_check='OK'; break;;
		[Nn]* ) Install_check='NO'; break;;
		* ) echo "Please answer yes (y) or no (n).";;
	esac
done

#Build SMART+
echo ""
cd ${current_dir}/build
cmake .. -DCMAKE_INCLUDE_PATH=${anacondaloc}/include -DCMAKE_LIBRARY_PATH=${anacondaloc}/lib -DCMAKE_INSTALL_PREFIX=${anacondaloc} -Wno-dev
echo ""
make -j${ncpus}
Install_OK=$?
echo ""
if [ $Install_OK -eq 0 ]
then

	if [ "${Install_check}" = "OK" ]
	then
		make install
		
	fi

    if [ $test -eq 1 ]
    then
        make test
    fi
    Test_OK=$?

	#Create the list of the file to copy after compilation
	executableToCopy="solver identification L_eff Elastic_props ODF PDF Abaqus_apply_inter"
#	objectToCopy="umat_single umat_singleT"

	# Copy all important files (+ final message)
	if [ $Test_OK -eq 0 ]
	then
		echo "\n---------------------------"
		
		#Treatement of object files
#		for object in ${objectToCopy}
#		do
#			#Copy of the "object".o from build/CMakeFiles/umat.dir/software to build/bin
#			if [ -f ${current_dir}/build/CMakeFiles/umat.dir/software/${object}.cpp.o ]
#			then
#				cp ${current_dir}/build/CMakeFiles/umat.dir/software/${object}.cpp.o ${current_dir}/build/bin/${object}.o
#				echo "${blue}${object}.o${reset} copied in ${blue}${current_dir}/build/bin${reset}"
#			fi
#		done
		
		#Treatement of executable files
		for file in ${executableToCopy}
		do
			#if debug exists, copy of the file from build/bin/Debug to build/bin
			if [ -f ${current_dir}/build/bin/Debug/${file} ]
			then
				cp ${current_dir}/build/bin/Debug/${file} ${current_dir}/build/bin
			fi

			#if Release exists, copy of the file from build/bin/Debug to build/bin
			if [ -f ${current_dir}/build/bin/Release/${file} ]
			then
				cp ${current_dir}/build/bin/Release/${file} ${current_dir}/build/bin
			fi
			
			#Copy the file from build/bin to exec
			cp ${current_dir}/build/bin/${file} ${current_dir}/exec/
			echo "${blue}${file}${reset} copied in ${blue}${current_dir}/exec${reset}"
		done

        cp ${current_dir}/build/bin/solver ${current_dir}/examples/elastic-plastic_tension
        cp ${current_dir}/build/bin/solver ${current_dir}/examples/micromechanics
        cp ${current_dir}/build/bin/identification ${current_dir}/examples/multi-layer_identification

		if [ "${Install_check}" = "OK" ]
		then
			echo "${blue}libsimcoon.so${reset} installed in ${blue}${current_dir}/lib${reset}"
		else
			echo "${blue}libsimcoon.so${reset} not installed."
		fi
		
		echo "---------------------------"
		echo "Simcoon compilation done.\n"
	else
		echo "\n---------------------------"
		echo "${red} Simcoon tests failed.\n${reset}"
	fi

else

	echo "\n---------------------------"
	echo "${red} Simcoon compilation failed.\n${reset}"

fi

if [ "${Install_check}" = "OK" ]
then
    cd ${current_dir}/simcoon-python-builder

    #Test if build exist and if it's necessary to erase it
    if [ ! -d "build" ]
    then
    mkdir ${current_dir}/simcoon-python-builder/build
    echo "Folder created.\n"
    else
    echo "Build directory already exists."

    while true; do
    read -p "Do you want to erase old compilation files (Recommended : No) ? " yn
    case $yn in
    [YyOo]* ) rm -r ${current_dir}/simcoon-python-builder/build/*; break;;
    [Nn]* ) break;;
    * ) echo "Please answer yes (y) or no (n).";;
    esac
    done
    fi

    cd ${current_dir}/simcoon-python-builder/build
    cmake .. -DCMAKE_INCLUDE_PATH=${anacondaloc}/include -DCMAKE_LIBRARY_PATH=${anacondaloc}/lib -DCMAKE_INSTALL_PREFIX=${anacondaloc} -Wno-dev
    echo ""
    make
    make install

    cd ..
    cd ..
    cp ${current_dir}/simcoon-python-builder/build/lib/simmit.so ${current_dir}/python-setup/simcoon/simmit.so
    cd ${current_dir}/python-setup

    #Change the current dir and install python library
    #current_dir=$(pwd)

    python setup.py install
    pip install .
    
    cd ${current_dir}/simcoon-python-builder/build
    make test
fi
