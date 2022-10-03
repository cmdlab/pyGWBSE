# *py*GWBSE
*py*GWBSE is a high-throughtput python workflow package designed to perform automtaed GW-BSE (Bethe-Salpeter Equation) calculation. The *py*GWBSE package leverages well known computational tools: *pymatgen, atomate, fireworks* and first principles software such as: *VASP* and *Wannier90* to perform high-throughput GW-BSE calculations. The workflow creates, computes, analyzes and stores all relevant simulation parameters and results in a queryable MongoDB database that can be accessed through our API. Check out our [documentation](https://cmdlab.github.io/pyGWBSE/) or the [paper](arxiv link needed)!

## Package Description
*py*GWBSE package takes input structure of materials and performs GW-BSE simulations for studying excited state properties. It performs automated convergence calculations required to obtain accurate results from such simulations. *py*GWBSE uses widely used *VASP* software to perform first-principles calculations. It uses *Wannier90* software to obtain quasiparticle (QP) bandstructure (both the one-shot G<sub>0</sub>W<sub>0</sub> and partially self-consistent GW<sub>0</sub> level) using maximally localized wannier functions. The electron and hole effective masses are computed using Sumo package. *py*GWBSE is also capable of solving BSE, which calculates absorption spectra of materials that incorporates excitonic effects and is accurate enough to comapre with experimental spectra. In addition to performing excited state property simulations, *py*GWBSE also computes several key electronic structure properties of materials such as, orbital resolved  density of states (both at DFT and QP level), real and imaginary part of the dielectric function (with and without incorporating electron-hole interaction), the exciton energies, and their corresponding oscillator strengths, band-edges, static dielectric tensors etc.


## Installation Instructions for *py*GWBSE
IMPORTANT NOTE: Atomate and FireWorks do not run on Windows OS. You need a unix-based OS (Mac or Linux) in order for these packages to run. As such, all setup instructions are given for Unix systems. 

1. Download the repo from the green code icon or via github's clone method.
- ``git clone git@github.com:cmdlab/pyGWBSE.git``
2. Install *py*GWBSE in a clean enviromnent using python=3.9. I suggest using Anaconda3 to manange environments. 
- ``conda create --name pygwbse python=3.9``
3. Activate the *py*GWBSE environment and run the line below in the *py*GWBSE directory to install:
- ``pip install -r requirements.txt``
4. After installation, *py*GWBSE needs to be added to your python path. This can be done by running the first line below **OR** by adding the 2nd line listed below to your *.bashrc* file. Only necessary if python cannot find the package or the setup.py failed for some reason.
- ``python setup.py develop`` or ``python setup.py install`` 
- ``export PYTHONPATH="$HOME/path_to_package/pyGWBSE:$PYTHONPATH"``
5. If this is your first time installing the package dependencies listed below, please ensure you have followed the respective setup instructions:
- [atomate](https://atomate.org/)  
- [FireWorks](https://materialsproject.github.io/fireworks/installation.html)
- [pymatgen](https://pymatgen.org/installation.html)
6. To run jupyter notebooks on various resources the ipykernel has to be installed. Sometimes this isn't enough and you need explicitly add the kernel to the list of environments. Via the command line:
- Activate your environment ``conda activate pygwbse``
- ``python -m ipykernel install --user --name pygwbse``

## Setting up dependancies
The *py*GWBSE package dependancies have a lot of documentation to look over. I will highlight the essential documentation to get started as quickly as possible.
1. *atomate* requires the most set up. Mainly, creating a directory scaffold and writing the 5 required files to connect to the database and run jobs. (MongoDB or free Atlas MongoDB is required) 
2. *pymatgen* has a command line tool installed to set up default directory paths called pmg. There are 2 essential commands you have to run to use *py*GWBSE on any system. 
- Reference directory for the VASP POTCARs. You need to have the POTCARs from VASP yourself.
  - `pmg config -p <EXTRACTED_VASP_POTCAR> <MY_PSP>` 
- Default pseudopotential files from VASP 
  - `pmg config --add PMG_DEFAULT_FUNCTIONAL PBE_54`
 3. *VASP*
 4. *Wannier90*
 5. *Sumo*

## Examples
To get started using *py*GWBSE, various tutorials and examples have been created using Jupyter Notebooks. These notebooks demonstrate the basic functionality of *py*GWBSE to enable users to quickly learn how to use the various modules within this package. These can be found under pyGWBSE/examples.

## Issues Installing
1. If you have a new install you likely do not have the gcc compiler that pymatgen requires. In that case run 
- `sudo apt install build-essential`
- `sudo apt-get install manpages-dev`

## How to cite *py*GWBSE
If you use *py*GWBSE in your research, please consider citing the paper!

> Tathagata Biswas, Arunima K. Singh. *pyGWBSE: A high throughput workflow
package for GW-BSE calculations.* Computational Materials Science, 2022, 207, 111238. [doi:10.1016/j.commatsci.2022.111238](https://doi.org/10.1016/j.commatsci.2022.111238)
