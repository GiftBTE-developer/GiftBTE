# Gift-BTE

## What is Gift-BTE 
Gift-BTE is a powerful tool designed to investigate submicron heat conduction. At the submicron scale, the characteristic length of heat conduction becomes comparable to the phonon mean free path, rendering the macroscopic Fourier's law inapplicable. Instead, at this scale, the phonon Boltzmann transport equation accurately describes heat conduction. Gift-BTE employs numerical methods to solve this equation and simulates submicron heat conduction. It takes phonon properties from first-principles simulations as input and provides a built-in database for some materials. Additionally, Gift-BTE offers an interface with two external unstructured mesh generators and provides several examples. The tool is easy to use, as it only requires an environment with g++, MPI, and Cmake.

## Features
### General
- Deterministic solver of phonon Boltzmann transport equation
- Applicable to any crystalline materials and structures
- Applicable to steady-state and transient problem
- Interface to ShengBTE, ALAMODE, GMSH, COMSOL, and Paraview packages
- Mainly written in C++, parallelized in MPI

### Application
- Computing thermal conductivity of nanostructures such as nano-porous media, superlattice, nanowires, nano-composite
- Predicting temperature rise in transistors
- Simulating laser heating process
- Other submicron heat conduction problems

## Workflow
![method-figure1](https://github.com/Gift-BTE-developer/Gift-BTE/assets/133620758/86174800-9b5c-4b00-9fab-a4643f14727f)

## Links
- Documentation: https://sjtu.feishu.cn/docx/GzB2dXQfaozenFxEtP4cQEJOnKb
- Git repository: https://github.com/Gift-BTE-developer/Gift-BTE.git

## License
This software is distributed under the GNU General Public License v3.0 license. See the LICENSE.txt file for license rights and limitations.

## How to cite Gift-BTE
Please cite the following article when you use Gift-BTE 
Yue Hu, Yongxing Shen, Hua Bao,
Ultra-efficient and parameter-free computation of submicron thermal transport with phonon Boltzmann transport equation,
Fundamental Research,
2022,
,
ISSN 2667-3258,
https://doi.org/10.1016/j.fmre.2022.06.007.
(https://www.sciencedirect.com/science/article/pii/S2667325822002758)

##Issues & Bug report
- If you find a bug or issue related to Gift-BTE, please report it at [Github Issues](https://github.com/Gift-BTE-developer/Gift-BTE/issues)
- Other questions and suggestions can be posted on the Github Disccusions [Github Disccusions](https://github.com/Gift-BTE-developer/Gift-BTE/discussions)

## Acknowledgment
This project is/was partially supported by the following projects:
 
- National Key R \& D Project from Ministery of Science and Technology of China (Grant No. 2022YFA1203100)
- National Natural Science Foundation of China (Grant No. 52122606)

## Contributors & Contact
Gift-BTE developers in TPEC Lab (Contact: Gift-BTE@outlook.com)

TPEC Lab: https://sites.ji.sjtu.edu.cn/hua-bao/

Global Institute of Future Technology (GIFT)

Shanghai Jiao Tong University

China

# Download
You can download the package from the git repository as:
 https://github.com/Gift-BTE-developer/Gift-BTE.git

If you choose the GitHub version, please use the ‘master’ branch.

# Installation
## Requirement
- C++ compiler (gcc is recommended)
- Cmake
- MPI library (openmpi is recommended)
-  Linux environment is recommended

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
 https://github.com/Gift-BTE-developer/Gift-BTE.git
 
Or directly clone the repository by git (may need installation)

- Git clone https://github.com/Gift-BTE-developer/Gift-BTE.git

- Git Checkout master



### Step 3. Build by CMake

Type in the command line under the Gift-BTE folder

cmake3 -B cpu-build -S. -DCMAKE_BUILD_TYPE=Release

![屏幕快照 2023-06-12 上午11 10 21](https://github.com/Gift-BTE-developer/Gift-BTE/assets/50352151/4f2baef9-7820-46c8-8094-50accf27e9f5)

Output in command line

![屏幕快照 2023-06-12 上午11 15 07](https://github.com/Gift-BTE-developer/Gift-BTE/assets/50352151/46b6867f-b4c5-4416-8481-08bc565b6f4f)


Common error: 

- C++ compiler not found: the C++ compiler is not installed 
- MPI not found: the MPI library is not installed 

### Step 4. Run example
Gift-BTE provides many examples under bin/examples. One can choose an example to test. For example: 
$ cd bin
$ cd examples
$ cd cross-plane
$ cd 1D
$ cd 1e-6
$ mpirun -np 4 ../../../../BTE_CPU

# More details relates to the User mannual can be found in https://sjtu.feishu.cn/docx/GzB2dXQfaozenFxEtP4cQEJOnKb
