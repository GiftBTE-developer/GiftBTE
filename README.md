# GiftBTE

## What is GiftBTE 
GiftBTE is a powerful tool designed to investigate submicron heat conduction. At the submicron scale, the characteristic length of heat conduction becomes comparable to the phonon mean free path, rendering the macroscopic Fourier's law inapplicable. Instead, at this scale, the phonon Boltzmann transport equation accurately describes heat conduction. GiftBTE employs numerical methods to solve this equation and simulates submicron heat conduction. It takes phonon properties from first-principles simulations as input and provides a built-in database for some materials. Additionally, GiftBTE offers an interface with two external unstructured mesh generators and provides several examples. The tool is easy to use, as it only requires an environment with g++, MPI, and Cmake.

### Video Links
https://www.bilibili.com/video/BV1BW4y1d7jn/?share_source=copy_web&vd_source=200fa4770a16926e2993543e680fef05

https://www.bilibili.com/video/BV1sX4y1J7bK/?share_source=copy_web&vd_source=200fa4770a16926e2993543e680fef05

https://www.bilibili.com/video/BV1r14y1Q7je/?share_source=copy_web&vd_source=200fa4770a16926e2993543e680fef05

https://www.bilibili.com/video/BV19u4y1S7ao/?share_source=copy_web&vd_source=200fa4770a16926e2993543e680fef05

## Features
### General
- Deterministic solver of phonon Boltzmann transport equation
- Applicable to any crystalline materials and structures
- Applicable to steady-state and transient problem
- Interface to ShengBTE, ALAMODE, GMSH, COMSOL, and Paraview packages
- Mainly written in C++, parallelized in MPI with CPUs or GPUs

### Application
- Computing thermal conductivity of nanostructures such as nano-porous media, superlattice, nanowires, nano-composite
- Predicting temperature rise in transistors
- Simulating laser heating process
- Other submicron heat conduction problems

## Workflow
![method-figure1](https://github.com/Gift-BTE-developer/Gift-BTE/assets/133620758/86174800-9b5c-4b00-9fab-a4643f14727f)

## Links
- Documentation: https://bte.sjtu.edu.cn/
- Git repository: https://github.com/GiftBTE-developer/GiftBTE.git

## License
This software is distributed under the GNU General Public License v3.0 license. See the LICENSE.txt file for license rights and limitations.

## How to cite GiftBTE
Please cite the following papers when you use GiftBTE:

Yue Hu, Ru Jia, Jiaxuan Xu, Yufei Sheng, Minhua Wen, James Lin, Yongxing Shen, Hua Bao, GiftBTE: An efficient deterministic solver for non-gray phonon Boltzmann transport equation, J. Phys.: Condens. Matter (2023). https://doi.org/10.1088/1361-648X/acfdea.

Yue Hu, Yongxing Shen, Hua Bao, Ultra-efficient and parameter-free computation of submicron thermal transport with phonon Boltzmann transport equation, Fundamental Research (2022). 
https://doi.org/10.1016/j.fmre.2022.06.007.

## Issues & Bug report
- If you find a bug or issue related to GiftBTE, please report it at [Github Issues](https://github.com/Gift-BTE-developer/Gift-BTE/issues)
- Other questions and suggestions can be posted on the Github Disccusions [Github Disccusions](https://github.com/GiftBTE-developer/GiftBTE/discussions)

## Acknowledgment
This project is/was partially supported by the following projects:
 
- National Key R \& D Project from Ministery of Science and Technology of China (Grant No. 2022YFA1203100)
- National Natural Science Foundation of China (Grant No. 52122606)

## Contributors & Contact
Authors: Yue Hu, Ru Jia, Jiaxuan Xu, Yufei Sheng, Minhua Wen, James Lin, Yongxing Shen, Hua Bao*

Mr. Peng Wan (Center for High Performance Computing, Shanghai Jiao Tong University) is the main contributor to the GPU version of GiftBTE

GiftBTE developers in TPEC Lab (Contact: Gift-BTE@outlook.com)

TPEC Lab: https://sites.ji.sjtu.edu.cn/hua-bao/

Global Institute of Future Technology (GIFT)

Shanghai Jiao Tong University

China

# Download
You can download the package from the git repository as: https://github.com/GiftBTE-developer/GiftBTE

Please download the ‘master’ branch.

# Installation
## Requirement
- C++ compiler (gcc is recommended)
- Cmake
- MPI library (openmpi is recommended)
- Linux environment is recommended

## Install and compile
### Step 1. Install all required packages
All packages are widely used packages. The installation of these packages can be easily found in web.

#### Ensure that the C++ compiler is installed 

Type in command line

![屏幕快照 2023-06-12 上午11 00 30](https://github.com/Gift-BTE-developer/Gift-BTE/assets/50352151/57fff5eb-bdf3-4e3e-8284-8b4685ab8950)

when output the version of gcc like:

![屏幕快照 2023-06-12 上午11 01 09](https://github.com/Gift-BTE-developer/Gift-BTE/assets/50352151/17aec67a-ebe0-4d38-9d28-7943da401bd4)

#### Ensure that the MPI library is installed 

Type in command line

![屏幕快照 2023-06-12 上午11 02 57](https://github.com/Gift-BTE-developer/Gift-BTE/assets/50352151/1f6d6c3a-41f2-4709-9e31-225361221bc1)

when output the folder of mpirun like

![屏幕快照 2023-06-12 上午11 03 30](https://github.com/Gift-BTE-developer/Gift-BTE/assets/50352151/dbdc40a0-437e-42f4-9fba-640599be3a3d)

#### Ensure that the Cmake is installed 

Type in command line

![屏幕快照 2023-06-12 上午11 04 59](https://github.com/Gift-BTE-developer/Gift-BTE/assets/50352151/5b961d6b-b7f7-4d22-828a-0318b2e2ef80)

when output the folder of cmake like

![屏幕快照 2023-06-12 上午11 05 03](https://github.com/Gift-BTE-developer/Gift-BTE/assets/50352151/07e2a9a2-0212-483b-9b4d-914a14d3440a)

### Step 2. Download source

You can download the package from the git repository as:
 https://github.com/GiftBTE-developer/GiftBTE
 - $ unzip GiftBTE-master.zip
 
Or directly clone the repository by git (may need installation)

- $ git clone https://github.com/GiftBTE-developer/GiftBTE.git




### Step 3. Build by CMake

Type in the command line under the GiftBTE folder

$ cd GiftBTE-master

$ cmake -B cpu-build -S.  -DCMAKE_BUILD_TYPE=Release

$ cd cpu-build

$ make

Output in command line

![1](https://github.com/GiftBTE-developer/GiftBTE/assets/133620758/0bcdc95d-3da5-449e-bb80-cff89f37a16e)
![2](https://github.com/GiftBTE-developer/GiftBTE/assets/133620758/b469b29f-7c00-4dee-b0f6-1eef223937e2)
![3](https://github.com/GiftBTE-developer/GiftBTE/assets/133620758/4dac43e2-5e82-4e31-a3d8-78d452fb4c62)


Common error: 

- C++ compiler not found: the C++ compiler is not installed 
- MPI not found: the MPI library is not installed
- If the CMake fails, please delete the cpu-build ($ rm -rf cpu-build) folder before retries.

### Step 4. Run example
GiftBTE provides many examples under bin/examples. One can choose an example to test. For example: 

$ cd bin

$ cd examples

$ cd cross-plane

$ cd 1e-6

$ mpirun -np 4 ../../../BTE_CPU

### About 3D cases
The current website version supports steady-state 2D case calculations.
For steady-state and transient-state (TDTR and TTG) 3D cases, you do not need to compile the GiftBTE code yourself — simply use the provided executable BTE_CPU file.

Examples for running 3D cases:

$ chmod -x BTE_CPU

$ cd bin/examples/Finfet/1e-7

$ mpirun -np 4 ../../../../BTE_CPU

The attached executable BTE_CPU file was built with the following environment. It is recommended to use the same or compatible versions:
- OpenMPI: 4.1
- GCC: 7.5
- CMake: 3.14

### The GPU version for transient-state cases (TTG and TDTR)
The compiled executable BTE_GPU file is also attached in the current website, which supports solving transient-state TTG and TDTR cases with GPUs. This executable BTE_GPU file was built with the following environment:
- OpenMPI: 4.1
- GCC: 12.3
- CMake: 3.27
- Cuda: 12.4

Examples for running:

$ chmod -x BTE_GPU

$ cd bin/examples/TDTR

$ export OMPI_MCA_btl=self,vader,tcp

$ mpirun -np 1 ../../../BTE_GPU

# More
More details can be found in [User's Mannual](https://sjtu.feishu.cn/docx/GzB2dXQfaozenFxEtP4cQEJOnKb)
