# coding: utf-8

import glob
import os

from atomate.utils.utils import env_chk
from atomate.vasp.database import VaspCalcDb
from fireworks import explicit_serialize, FiretaskBase, FWAction
from monty.json import jsanitize
from pymatgen.io.vasp.outputs import Vasprun, Outcar

from pyGWBSE.tasks import read_emcpyout, read_epsilon, get_gap_from_dict, read_vac_level
from pyGWBSE.wannier_tasks import read_vbm, read_wannier, read_vasp, read_special_kpts


@explicit_serialize
class gw2db(FiretaskBase):
    """
    Insert quasi-particle energies into the database for a GW calculation.
    """
    required_params = ["structure", "task_label", "db_file", "mat_name"]
    optional_params = ["job_tag"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        # get adddtional tags to parse the directory for
        db_file = env_chk(self.get('db_file'), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        ifconv = fw_spec["ifconv"]
        structure = self["structure"]
        task_label = self["task_label"]
        if "job_tag" in self:
            job_tag = self["job_tag"]
        else:
            job_tag = None
        mat_name = self["mat_name"]
        task_collection = 'QP_Results'
        dir_name = os.getcwd()
        file = glob.glob('vasprun.xml*')[-1]
        vasprun = Vasprun(file)
        file = glob.glob('OUTCAR*')[-1]
        outcar = Outcar(file)
        run_stats=outcar.run_stats
        qp_energies = vasprun.eigenvalues
        en, eps1, eps2 = vasprun.dielectric
        igap, dgap = get_gap_from_dict(qp_energies)
        incar = vasprun.incar
        parameters = vasprun.parameters
        bgap, cbm, vbm, is_band_gap_direct = vasprun.eigenvalue_band_properties
        kpts_dict = vasprun.kpoints.as_dict()
        # dictionary to update the database with
        d = {"structure": structure.as_dict(),
             "formula_pretty": structure.composition.reduced_formula,
             "material_id": mat_name, "run_stats": run_stats,
             "run_directory": dir_name, 'direct_gap': dgap, 'indirect_gap': igap, 
             "qp_energies": qp_energies, "task_label": task_label,
             "frequency": en, "epsilon_1": eps1, "epsilon_2": eps2, 
             "job_tag": job_tag, "ifconv": ifconv, "vbm": vbm, "cbm": cbm, 
             "incar": incar, "parameters": parameters, "kpoints": kpts_dict}
        d = jsanitize(d)
        coll = mmdb.db[task_collection]
        coll.insert_one(d)

        return FWAction(update_spec={"gw_gaps": [igap, dgap]})   


@explicit_serialize
class bse2db(FiretaskBase):
    """
    Insert exciton energies, oscillator strength and dielectric function into the database for a BSE calculation.
    """
    required_params = ["structure", "task_label", "db_file", "mat_name"]
    optional_params = ["job_tag"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        # get adddtional tags to parse the directory for
        db_file = env_chk(self.get('db_file'), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        structure = self["structure"]
        task_label = self["task_label"]
        job_tag = self["job_tag"]
        mat_name = self["mat_name"]
        igap = fw_spec["gw_gaps"][0]
        dgap = fw_spec["gw_gaps"][1]
        task_collection = 'BSE_Results'
        dir_name = os.getcwd()
        filename = glob.glob('vasprun.xml*')[-1]
        if "job_tag" in self:
            job_tag = self["job_tag"]
        else:
            job_tag = None
        with open(filename, "a") as file:
            file.write("</modeling>")
        vasprun = Vasprun(filename)
        incar = vasprun.incar
        parameters = vasprun.parameters
        optical_transition = vasprun.optical_transition
        en, eps1, eps2 = vasprun.dielectric
        kpts_dict = vasprun.kpoints.as_dict()
        file = glob.glob('OUTCAR*')[-1]                                         
        outcar = Outcar(file)                                                   
        run_stats=outcar.run_stats                                              
        d = {"structure": structure.as_dict(),
             "formula_pretty": structure.composition.reduced_formula,
             "material_id": mat_name, "run_directory": dir_name,
             "frequency": en, "epsilon_1": eps1, "epsilon_2": eps2,
             'direct_gap': dgap, 'indirect_gap': igap, "run_stats": run_stats,
             "optical_transition": optical_transition, "task_label": task_label, "job_tag": job_tag,
             "incar": incar, "parameters": parameters, "kpoints": kpts_dict}
        d = jsanitize(d)
        coll = mmdb.db[task_collection]
        coll.insert_one(d)

@explicit_serialize
class rpa2db(FiretaskBase):
    """
    Insert exciton energies, oscillator strength and dielectric function into the database for a BSE calculation.
    """
    required_params = ["structure", "task_label", "db_file", "mat_name"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        # get adddtional tags to parse the directory for
        db_file = env_chk(self.get('db_file'), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        structure = self["structure"]
        task_label = self["task_label"]
        mat_name = self["mat_name"]
        task_collection = 'RPA_Results'
        dir_name = os.getcwd()
        filename = glob.glob('vasprun.xml*')[-1]
        vasprun = Vasprun(filename)
        incar = vasprun.incar
        parameters = vasprun.parameters
        dft_energies = vasprun.eigenvalues
        en, eps1, eps2 = vasprun.dielectric
        kpts_dict = vasprun.kpoints.as_dict()
        file = glob.glob('OUTCAR*')[-1]                                         
        outcar = Outcar(file)                                                   
        run_stats=outcar.run_stats                                              
        d = {"structure": structure.as_dict(),
             "formula_pretty": structure.composition.reduced_formula,
             "material_id": mat_name, "run_directory": dir_name,
             "dft_energies": dft_energies, "run_stats": run_stats,
             "frequency": en, "epsilon_1": eps1, "epsilon_2": eps2, "task_label": task_label,
             "incar": incar, "parameters": parameters, "kpoints": kpts_dict}
        d = jsanitize(d)
        coll = mmdb.db[task_collection]
        coll.insert_one(d)

@explicit_serialize
class emc2db(FiretaskBase):
    """
    Insert effective masses for a SUMO-BANDSTATS calculation.
    """
    required_params = ["structure", "db_file", "mat_name"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        # get adddtional tags to parse the directory for
        db_file = env_chk(self.get('db_file'), fw_spec)
        dir_name = os.getcwd()
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        structure = self["structure"]
        mat_name = self["mat_name"]
        task_collection = 'EMC_Results'
        filename = glob.glob('sumo-bandstats.log*')[-1]
        hmass, emass = read_emcpyout(filename)
        # dictionary to update the database with
        d = {"structure": structure.as_dict(),
            "formula_pretty": structure.composition.reduced_formula,
            "material_id": mat_name, "run_directory": dir_name,
            "hole_effective_mass": hmass, "electron_effective_mass": emass}
        d = jsanitize(d)
        coll = mmdb.db[task_collection]
        coll.insert_one(d)

@explicit_serialize
class eps2db(FiretaskBase):
    """
    Insert macroscopic dielectric constants for LEPSILON=TRUE calculation.
    """
    required_params = ["structure", "db_file", "mat_name"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        # get additional tags to parse the directory for
        db_file = env_chk(self.get('db_file'), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        dir_name = os.getcwd()
        structure = self["structure"]
        mat_name = self["mat_name"]
        task_collection = 'EPS_Results'
        filename = glob.glob('vasprun.xml*')[-1]
        try:
            locpot_fname = glob.glob('LOCPOT*')[-1]
            zvac, evac, delta_evac=read_vac_level(locpot_fname, filename)
            ifvac=True
        except:
            ifvac=False
        epsilon = read_epsilon(filename)
        vrun = Vasprun(filename, parse_projected_eigen=True)
        ks_energies = vrun.eigenvalues
        proj_eigs=vrun.projected_eigenvalues
        kwgs=vrun.actual_kpoints_weights
        file = glob.glob('OUTCAR*')[-1]                                         
        outcar = Outcar(file)                                                   
        run_stats=outcar.run_stats                                              
        igap, dgap = get_gap_from_dict(ks_energies)
        bgap, cbm, vbm, is_band_gap_direct = vrun.eigenvalue_band_properties
        # dictionary to update the database with
        d = {"structure": structure.as_dict(),
                "formula_pretty": structure.composition.reduced_formula,
                "material_id": mat_name, "dielectric constant": epsilon, 
                "run_stats": run_stats, "run_directory": dir_name, 
                'direct_gap': dgap, 'indirect_gap': igap,
                "kpoint_weights": kwgs, "cbm": cbm, "vbm": vbm,
                "projected_eigs": proj_eigs, "ks_energies": ks_energies}
        if ifvac:
            d.update({"z_vacuum": zvac, "e_vacuum": evac, "delta_e_vacuum": delta_evac})
        d = jsanitize(d)
        coll = mmdb.db[task_collection]
        coll.insert_one(d)



@explicit_serialize
class Wannier2DB(FiretaskBase):
    """
    Insert a optical transition and quasi-particle energies
    into the database for a GW and BSE calculation.

    Args:
        db_file (str): path to file containing the database credentials.
            Supports env_chk. Default: write data to JSON file.
        task_collection (str): The name of the task collection you
            want to push the data to.
        wf_name (str): The name of the workflow that this analysis is part of.
    """

    required_params = ["structure", "task_label", "db_file", "compare_vasp", "mat_name"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        # get adddtional tags to parse the directory for
        db_file = env_chk(self.get('db_file'), fw_spec)
        mmdb = VaspCalcDb.from_db_file(db_file, admin=True)
        structure = self["structure"]
        task_label = self["task_label"]
        compare_vasp = self["compare_vasp"]
        mat_name = self["mat_name"]
        task_collection = 'WANNIER_Results'
        dir_name = os.getcwd()
        fname_band = glob.glob('wannier90_band.dat*')[-1]
        fname_kpt = glob.glob('wannier90_band.kpt*')[-1]
        fname_gnu = glob.glob('wannier90_band.gnu*')[-1]
        fname_vasp = glob.glob('vasprun.xml*')[-1]
        gap, vbm = read_vbm(fname_vasp)
        kpts, eigs_wann = read_wannier(fname_band, fname_kpt, vbm)
        spkptl, spkptc = read_special_kpts(fname_gnu)
        if compare_vasp:
            kpt_vasp, eigs_vasp = read_vasp(fname_vasp, vbm)
            d = {"structure": structure.as_dict(),
                 "formula_pretty": structure.composition.reduced_formula,
                 "material_id": mat_name,
                 "run_directory": dir_name,
                 "wannier_kpoints": kpts, "wannier_eigenvalues": eigs_wann,
                 "actual_kpoints": kpt_vasp, "actual_eigenvalues": eigs_vasp,
                 "special_kpoint_labels": spkptl,
                 "special_kpoint_coordinates": spkptc,
                 "task_label": task_label}
        else:
            d = {"structure": structure.as_dict(),
                 "formula_pretty": structure.composition.reduced_formula,
                 "material_id": mat_name,
                 "run_directory": dir_name,
                 "wannier_kpoints": kpts, "wannier_eigenvalues": eigs_wann,
                 "special_kpoint_labels": spkptl,
                 "special_kpoint_coordinates": spkptc,
                 "task_label": task_label}
        d = jsanitize(d)
        coll = mmdb.db[task_collection]
        coll.insert_one(d)
