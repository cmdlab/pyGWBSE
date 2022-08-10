# coding: utf-8

import glob
import gzip
import os
import re
from atomate.common.firetasks.glue_tasks import CopyFiles, get_calc_loc
from atomate.utils.utils import env_chk, get_logger
from fireworks import explicit_serialize, FiretaskBase, FWAction
from pymatgen.io.vasp import Vasprun, Locpot
from pymatgen.io.vasp.inputs import Incar
from pymatgen.core import Structure

from pyGWBSE.inputset import CreateInputs

"""
This module defines tasks that acts as a glue between other vasp Firetasks to allow communication
between different Firetasks and Fireworks. This module also contains tasks that affect the control
flow of the workflow, e.g. tasks to check stability or the gap is within a certain range.
"""

logger = get_logger(__name__)

__author__ = 'Anubhav Jain, Kiran Mathew'
__email__ = 'ajain@lbl.gov, kmathew@lbl.gov'


@explicit_serialize
class CheckBeConv(FiretaskBase):
    """
    Your Comments Here
    """
    required_params = ["niter", "tolerence", "no_conv"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        niter = self["niter"]
        conv = self["no_conv"]
        tol = self["tolerence"]
        filename = str(os.getcwd()) + '/vasprun.xml'
        vasprun = Vasprun(filename)
        gap, cbm, vbm, is_direct = vasprun.eigenvalue_band_properties
        print()
        print('=================================================')
        print('Iteration:', niter)
        print("band gap:", ("%10.4f" % gap))
        print("CBM:", ("%10.4f" % cbm))
        print("VBM:", ("%10.4f" % vbm))
        print("Is the gap direct:", is_direct)
        if niter > 1:
            old_val = fw_spec["conval"]
            new_val = gap
            diff = abs(old_val - new_val)
            print("Difference:", ("%10.4f" % diff))
            print('=================================================')
            print()
            if diff < tol:
                conv = True
        return FWAction(update_spec={"ifconv": conv, "conval": gap})


@explicit_serialize
class MakeWFilesList(FiretaskBase):

    def run_task(self, fw_spec):
        files2copy = glob.glob("W*.tmp")
        return FWAction(update_spec={"wfiles": files2copy})

@explicit_serialize
class SaveNbandsov(FiretaskBase):

    required_params = ["enwinbse"]

    def run_task(self, fw_spec):

        enwinbse = self["enwinbse"]
        filename = str(os.getcwd()) + '/vasprun.xml'
        vasprun = Vasprun(filename)
        qpe=vasprun.eigenvalues
        gap, cbm, vbm, is_direct = vasprun.eigenvalue_band_properties
        nbandso,nbandsv=get_nbandsov(qpe,vbm,cbm,enwinbse)   

        return FWAction(update_spec={"nbandso": nbandso, "nbandsv": nbandsv})

@explicit_serialize
class SaveConvParams(FiretaskBase):

    required_params = ["nomegagw","encutgw","nbands"]

    def run_task(self, fw_spec):
        nomegagw=self["nomegagw"]
        encutgw=self["encutgw"]
        nbands=self["nbands"]
        return FWAction(update_spec={"nomegagw": nomegagw, "encutgw": encutgw, "nbands": nbands})

@explicit_serialize
class StopIfConverged(FiretaskBase):
    """
    Your Comments Here
    """
    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        ifconv = fw_spec["ifconv"]
        if ifconv:
            return FWAction(exit=True)


@explicit_serialize
class PasscalClocsCond(FiretaskBase):
    """
    Your Comments Here
    """
    required_params = ["name"]
    optional_params = ["filesystem", "path"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        calc_locs = list(fw_spec.get("calc_locs", []))
        calc_locs.append({"name": self["name"],
                          "filesystem": env_chk(self.get('filesystem', None), fw_spec),
                          "path": self.get("path", os.getcwd())})
        ifconv = fw_spec["ifconv"]
        if ifconv:
            return FWAction(mod_spec=[{'_push_all': {'calc_locs': calc_locs}}])


@explicit_serialize
class WriteBSEInput(FiretaskBase):
    """
    Your Comments Here
    """
    required_params = ["structure", "reciprocal_density"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        f_incar = str(os.getcwd()) + '/INCAR'
        structure = self["structure"]
        reciprocal_density = self["reciprocal_density"]
        nbandso = fw_spec["nbandso"]
        nbandsv = fw_spec["nbandsv"]
        prev_incar = Incar.from_file(f_incar)
        vis = CreateInputs(structure, mode='BSE', prev_incar=prev_incar, reciprocal_density=reciprocal_density,
                           nbandso=nbandso, nbandsv=nbandsv)
        vis.write_input(".")


@explicit_serialize
class WriteGWInput(FiretaskBase):
    """
    Your Comments Here
    """
    required_params = ["structure", "reciprocal_density", "nbandsgw", "wannier_fw"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        f_incar = str(os.getcwd()) + '/INCAR'
        prev_incar = Incar.from_file(f_incar)
        structure = self["structure"]
        nbandsgw = self["nbandsgw"]
        wannier_fw = self["wannier_fw"]
        reciprocal_density = self["reciprocal_density"]
        encutgw=fw_spec["encutgw"]
        nomegagw=fw_spec["nomegagw"]
        nbands=fw_spec["nbands"]
        vis = CreateInputs(structure, mode='GW', prev_incar=prev_incar, reciprocal_density=reciprocal_density,
                           encutgw=encutgw, nomegagw=nomegagw, nbands=nbands, nbandsgw=nbandsgw, wannier_fw=wannier_fw)
        vis.write_input(".")

def read_emcpyout(fname):                                                       
    f = open(fname)                                                             
    contents = f.readlines()                                                    
    f.close()                                                                   
    ne=0                                                                        
    nh=0                                                                        
    h_res={}                                                                    
    e_res={}                                                                    
    for content in contents:                                                    
        if "m_h:" in content:                                                   
            nh=nh+1                                                             
            lines = str.split(content)                                          
                                                                                
            hmass=eval(lines[1])                                                
            hibnd=eval(lines[4])                                                
            hdir1=lines[9]                                                      
            hdir2=lines[14]                                                     
                                                                                
            h_res['mass'+str(nh)+': '+str(hibnd)+', '+hdir1+', '+hdir2]=hmass            
                                                                                
        if "m_e:" in content:                                                   
            ne=ne+1                                                             
            lines = str.split(content)                                          
                                                                                
            emass=eval(lines[1])                                                
            eibnd=eval(lines[4])                                                
            edir1=lines[9]                                                      
            edir2=lines[14]                                                     
                                                                                
            e_res['mass'+str(ne)+': '+str(eibnd)+', '+edir1+', '+edir2]=emass            
                                                                                
    return h_res, e_res     

def read_epsilon(fname):
    vasprun = Vasprun(fname)
    try:
        epsilon = vasprun.epsilon_static_wolfe
    except:
        epsilon = None
    return epsilon


def calc_delta_evac(x,y,n0):
    dx1=x[n0+1]-x[n0]
    dy1=y[n0+1]-y[n0]
    dx2=x[n0]-x[n0-1]
    dy2=y[n0]-y[n0-1]
    f1=dy1/dx1
    f2=dy2/dx2
    return (f1+f2)/2



def read_vac_level(locpot_fname,fname):
    lpot=Locpot.from_file(locpot_fname)
    vrun=Vasprun(fname)
    vrun_d=vrun.as_dict()
    struct=Structure.from_dict(vrun_d["input"]["crystal"])
    latt_c=struct.lattice.c
    y=lpot.get_average_along_axis(2)
    zarr=lpot.get_axis_grid(2)

    dist_z0=0.
    z0=0.

    for n,z in enumerate(zarr):
        dist=[]
        for site in struct.sites:
            pos_z=site.coords[2]
            if pos_z<0:
                pos_z=pos_z+latt_c
            dist.append(abs(pos_z-z))
            dist.append(abs(pos_z+latt_c-z))
        dist_z=min(dist)
        if dist_z>dist_z0:
            dist_z0=dist_z
            z0=z
            n0=n

    evac=y[n0]
    zvac=zarr[n0]
    delta_evac=calc_delta_evac(zarr,y,n0)

    return zvac,evac,delta_evac

def get_gap_from_dict(qp_energies):

    tol=1e-5
    vbm = -float("inf")
    cbm = float("inf")
    for spin, d in qp_energies.items():
        for k, val in enumerate(d):
            for (eigenval, occu) in val:
                if occu > tol and eigenval > vbm:
                    vbm = eigenval
                elif occu <= tol and eigenval < cbm:
                    cbm = eigenval
    igap=max(cbm-vbm,0)

    gap = []
    filled = 1.0
    tol=1e-5
    vbm=0
    cbm=0
    for spin, d in qp_energies.items():
        for k, val in enumerate(d):
            for (eigenval, occu) in val:
                if abs(occu-filled)<=tol:
                    vbm=eigenval
                else:
                    cbm=eigenval
                    gap.append(cbm-vbm)
                    break
    dgap=min(gap)
    return igap,dgap

def get_nbandsov(qp_energies,vbm,cbm,enwinbse):                                      
                                                                                
    tol=1e-5                                                                    
    nvb=[]                                                                      
    ncb=[]                                                                      
    for spin, d in qp_energies.items():                                         
        for k, val in enumerate(d):                                             
            for ibnd,(eigenval, occ) in enumerate(val):                         
                if abs(occ-1.0)<tol and (vbm-eigenval)<=enwinbse:                  
                    nvb.append(ibnd+1)                                          
                if occ<tol and (eigenval-cbm)<=enwinbse:                           
                    ncb.append(ibnd+1)                                          
                                                                                
    return max(nvb)-min(nvb)+1,max(ncb)-min(ncb)+1                              
                                                                                


@explicit_serialize
class CopyOutputFiles(CopyFiles):
    """
    Copy files from a previous VASP run directory to the current directory.
    By default, copies 'INCAR', 'POSCAR' (default: via 'CONTCAR'), 'KPOINTS', 
    'POTCAR', 'OUTCAR', and 'vasprun.xml'. Additional files, e.g. 'CHGCAR', 
    can also be specified. Automatically handles files that have a ".gz" 
    extension (copies and unzips).

    Note that you must specify either "calc_loc" or "calc_dir" to indicate
    the directory containing the previous VASP run.

    Args:
        none (none): but you must specify either "calc_loc" OR "calc_dir"

    Other Parameters:
        calc_loc (str OR bool): if True will set most recent calc_loc. If str
            search for the most recent calc_loc with the matching name
        calc_dir (str): path to dir that contains VASP output files.
        filesystem (str): remote filesystem. e.g. username@host
        additional_files ([str]): additional files to copy,
            e.g. ["CHGCAR", "WAVECAR"]. Use $ALL if you just want to copy
            everything
        contcar_to_poscar(bool): If True (default), will move CONTCAR to
            POSCAR (original POSCAR is not copied).
    """

    optional_params = ["calc_loc", "calc_dir", "filesystem", "additional_files",
                       "contcar_to_poscar"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        calc_loc = get_calc_loc(self["calc_loc"],
                                fw_spec["calc_locs"]) if self.get(
            "calc_loc") else {}

        # determine what files need to be copied
        files_to_copy = None
        if not "$ALL" in self.get("additional_files", []):
            files_to_copy = ['INCAR', 'POSCAR', 'KPOINTS', 'POTCAR', 'OUTCAR',
                             'vasprun.xml']
            if self.get("additional_files"):
                files_to_copy.extend(self["additional_files"])
            if "wfiles" in fw_spec.keys():
                files2copy = fw_spec["wfiles"]
                print(files2copy)
                files_to_copy.extend(files2copy)

        # decide between poscar and contcar
        contcar_to_poscar = self.get("contcar_to_poscar", True)
        if contcar_to_poscar and "CONTCAR" not in files_to_copy:
            files_to_copy.append("CONTCAR")
            files_to_copy = [f for f in files_to_copy if
                             f != 'POSCAR']  # remove POSCAR

        # setup the copy
        self.setup_copy(self.get("calc_dir", None),
                        filesystem=self.get("filesystem", None),
                        files_to_copy=files_to_copy, from_path_dict=calc_loc)
        # do the copying
        self.copy_files()

    def copy_files(self):
        """
        Your Comments Here
        """
        all_files = self.fileclient.listdir(self.from_dir)
        # start file copy
        for f in self.files_to_copy:
            prev_path_full = os.path.join(self.from_dir, f)
            print(prev_path_full)
            dest_fname = 'POSCAR' if f == 'CONTCAR' and self.get(
                "contcar_to_poscar", True) else f
            dest_path = os.path.join(self.to_dir, dest_fname)

            relax_ext = ""
            relax_paths = sorted(
                self.fileclient.glob(prev_path_full + ".relax*"))
            if relax_paths:
                if len(relax_paths) > 9:
                    raise ValueError(
                        "CopyVaspOutputs doesn't properly handle >9 relaxations!")
                m = re.search('\.relax\d*', relax_paths[-1])
                relax_ext = m.group(0)

            # detect .gz extension if needed - note that monty zpath() did not seem useful here
            gz_ext = ""
            if not (f + relax_ext) in all_files:
                for possible_ext in [".gz", ".GZ"]:
                    if (f + relax_ext + possible_ext) in all_files:
                        gz_ext = possible_ext

            if not (f + relax_ext + gz_ext) in all_files:
                raise ValueError("Cannot find file: {}".format(f))

            # copy the file (minus the relaxation extension)
            self.fileclient.copy(prev_path_full + relax_ext + gz_ext,
                                 dest_path + gz_ext)

            # unzip the .gz if needed
            if gz_ext in ['.gz', ".GZ"]:
                # unzip dest file
                f = gzip.open(dest_path + gz_ext, 'rt')
                file_content = f.read()
                with open(dest_path, 'w') as f_out:
                    f_out.writelines(file_content)
                f.close()
                os.remove(dest_path + gz_ext)
