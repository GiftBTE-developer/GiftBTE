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
You can download the latest and previous versions of Gift-BTE at
 https://bte.sjtu.edu.cn/
 
You can also download the package from the git repository as:
 https://github.com/Gift-BTE-developer/Gift-BTE.git

If you choose the GitHub version, please use the ‘master’ branch.
