# pyGWBSE
python workflow for GW-BSE calculation
The Hetero2d package leverages well known computational tools: pymatgen, MPInterfaces, atomate, fireworks, and custodian to perform high-throughput *ab-initio* calculations. Hetero2d is tailored to addressing scientific questions regarding the stability of 2D-substrate hetero-structured materials using an all-in-one workflow approach to model the hetero-structures formed by arbitrary 2D materials and substrates. The workflow creates, computes, analyzes and stores all relevant simulation parameters and results in a queryable MongoDB database that can be accessed through our API. Check out our [documentation](https://cmdlab.github.io/Hetero2d/) or the [paper](https://doi.org/10.1016/j.commatsci.2022.111238)!

## Package Description
The 2D-substrate hetero-structure workflow takes a given 2D material, 3D phase (bulk) of the 2D material, and a substrate, relaxes the structures using vdW-corrected DFT and creates hetero-structures subject to user constraints. The workflow analyzes and stores energetic stability information of various 2D hetero-structured materials to predict the feasibility of a substrate stabilizing a meta-stable 2D material

Hetero2d provides automated routines for the generation of low-lattice mismatched hetero-structures for arbitrary 2D materials and substrate surfaces, the creation of van der Waals corrected density-functional theory (DFT) input files, the submission and monitoring of simulations on computing resources, the post-processing of the key parameters to compute (a) the interface interaction energy of 2D-substrate hetero-structures, (b) the identification of substrate-induced changes in interfacial structure, and (c) charge doping of the 2D material.

## Installation Instructions for Hetero2d
IMPORTANT NOTE: Atomate and FireWorks do not run on Windows OS. You need a unix-based OS (Mac or Linux) in order for these packages to run. As such, all setup instructions are given for Unix systems. 

1. Download the repo from the green code icon or via github's commandline tool gh. 
- ``gh repo clone cmdlab/Hetero2d`` (gh must be installed)
- ``git clone https://github.com/cmdlab/Hetero2d.git``
2. Install Hetero2d in a clean enviromnent using python=3.9. I suggest using Anaconda3 to manange environments. 
- ``conda create --name hetero2d python=3.9``
3. Activate the Hetero2d environment and run the line below in the Hetero2d directory to install:
- ``pip install -r requirements.txt``
4. After installation, Hetero2d needs to be added to your python path. This can be done by running the first line below **OR** by adding the 2nd line listed below to your *.bashrc* file. Only necessary if python cannot find the package or the setup.py failed for some reason.
- ``python setup.py develop`` or ``python setup.py install`` 
- ``export PYTHONPATH="$HOME/path_to_package/Hetero2d:$PYTHONPATH"``
5. If this is your first time installing the package dependencies listed below, please ensure you have followed the respective setup instructions:
- [atomate](https://atomate.org/)  
- [FireWorks](https://materialsproject.github.io/fireworks/installation.html)
- [pymatgen](https://pymatgen.org/installation.html)
6. To run jupyter notebooks on various resources the ipykernel has to be installed. Sometimes this isn't enough and you need explicitly add the kernel to the list of environments. Via the command line:
- Activate your environment ``conda activate hetero2d``
- ``python -m ipykernel install --user --name hetero2d``

## Setting up dependancies
The Hetero2d package dependancies have a lot of documentation to look over. I will highlight the essential documentation to get started as quickly as possible.
1. *atomate* requires the most set up. Mainly, creating a directory scaffold and writing the 5 required files to connect to the database and run jobs. (MongoDB or free Atlas MongoDB is required) 
2. *pymatgen* has a command line tool installed to set up default directory paths called pmg. There are 2 essential commands you have to run to use Hetero2d on any system. 
- Reference directory for the VASP POTCARs. You need to have the POTCARs from VASP yourself.
  - `pmg config -p <EXTRACTED_VASP_POTCAR> <MY_PSP>` 
- Default pseudopotential files from VASP 
  - `pmg config --add PMG_DEFAULT_FUNCTIONAL PBE_54`

## Examples
To get started using Hetero2d, various tutorials and examples have been created using Jupyter Notebooks. These notebooks demonstrate the basic functionality of Hetero2d to enable users to quickly learn how to use the various modules within this package. These can be found under Hetero2d/examples.

## Issues Installing
1. If you have a new install you likely do not have the gcc compiler that pymatgen requires. In that case run 
- `sudo apt install build-essential`
- `sudo apt-get install manpages-dev`

## How to cite Hetero2d
If you use Hetero2d in your research, please consider citing the paper!

> Tara M. Boland, Arunima K. Singh. *Computational synthesis of 2D materials: A high-throughput approach to materials design.* Computational Materials Science, 2022, 207, 111238. [doi:10.1016/j.commatsci.2022.111238](https://doi.org/10.1016/j.commatsci.2022.111238)
