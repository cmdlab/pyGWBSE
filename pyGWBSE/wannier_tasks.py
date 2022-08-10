# coding: utf-8

"""
This module defines tasks for WANNIER90 calculation 
"""

import os

import numpy as np
from atomate.utils.utils import get_logger
from fireworks import explicit_serialize, FiretaskBase
from pymatgen import Structure
from pymatgen.io.vasp import Vasprun
from pymatgen.io.vasp.inputs import Incar, Potcar, PotcarSingle
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath

from pyGWBSE.inputset import CreateInputs

logger = get_logger(__name__)

__author__ = 'Tathagata Biswas'
__email__ = 'tbiswas3@asu.edu'


@explicit_serialize
class WriteWannierInputForGW(FiretaskBase):
    """
    Your Comments Here
    """
    required_params = ["structure", "reciprocal_density", "nbandsgw"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        f_incar = str(os.getcwd()) + '/INCAR'
        structure = self["structure"]
        nbandsgw = self["nbandsgw"]
        reciprocal_density = self["reciprocal_density"]
        prev_incar = Incar.from_file(f_incar)
        wann_inp = str(os.getcwd()) + '/wannier90.win'
        vasprunfile = str(os.getcwd()) + '/vasprun.xml'
        poscarfile = str(os.getcwd()) + '/POSCAR'
        potcarfile = str(os.getcwd()) + '/POTCAR'
        vasprun = Vasprun(vasprunfile)
        incar = vasprun.incar
        nbands = incar["NBANDS"]
        elements = read_potcar(potcarfile, poscarfile)
        labels, kpts = kpath_finder(poscarfile)
        numwan = 0
        for element in elements:
            numwan = numwan + element[5]
        if numwan > nbandsgw:
            nbandsgw = numwan
        encutgw=fw_spec["encutgw"]
        nomegagw=fw_spec["nomegagw"]
        nbands=fw_spec["nbands"]
        vis = CreateInputs(structure, mode='GW', reciprocal_density=reciprocal_density,
                           prev_incar=prev_incar, encutgw=encutgw, nomegagw=nomegagw, nbands=nbands, nbandsgw=nbandsgw)
        vis.write_input(".")
        labels, kpts = kpath_finder(poscarfile)
        write_wannier_input(numwan, nbands, labels, kpts, wann_inp, elements,False)


@explicit_serialize
class WriteWannierInputForDFT(FiretaskBase):
    """
    Your Comments Here
    """
    required_params = ["structure", "reciprocal_density", "ppn", "write_hr"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        f_incar = str(os.getcwd()) + '/INCAR'
        structure = self["structure"]
        reciprocal_density = self["reciprocal_density"]
        ppn = self["ppn"]
        write_hr = self["write_hr"]
        prev_incar = Incar.from_file(f_incar)
        wann_inp = str(os.getcwd()) + '/wannier90.win'
        vasprunfile = str(os.getcwd()) + '/vasprun.xml'
        poscarfile = str(os.getcwd()) + '/POSCAR'
        potcarfile = str(os.getcwd()) + '/POTCAR'
        vasprun = Vasprun(vasprunfile)
        incar = vasprun.incar
        nbands = incar["NBANDS"]
        elements = read_potcar(potcarfile, poscarfile)
        numwan = 0
        for element in elements:
            numwan = numwan + element[5]
        print(numwan, nbands)
        if numwan > nbands:
            nbands = (int(numwan / ppn) + 1) * ppn
        vis = CreateInputs(structure, mode='STATIC', prev_incar=prev_incar, reciprocal_density=reciprocal_density,
                           nbands=nbands)
        vis.write_input(".")
        labels, kpts = kpath_finder(poscarfile)
        write_wannier_input(numwan, nbands, labels, kpts, wann_inp, elements,write_hr)


@explicit_serialize
class CopyKptsWan2vasp(FiretaskBase):
    """
    Your Comments Here
    """

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        f_wannkpt = str(os.getcwd()) + '/wannier90_band.kpt'
        f_vaspkpt = str(os.getcwd()) + '/KPOINTS'
        f = open(f_wannkpt)
        contents = f.readlines()
        f.close()
        lines = str.split(contents[0])
        nkpts = eval(lines[0])
        f = open(f_vaspkpt, 'w')
        f.write('kpoints file generated from wannier90_band.kpt' + '\n')
        f.write(str(nkpts) + '\n')
        f.write('Reciprocal' + '\n')
        for ik in range(nkpts):
            f.write(contents[ik + 1])
        f.close()


def write_wannier_input(numwan, nbands, labels, kpts, wann_inp, elements,write_hr):
    """
    Your Comments Here
    """
    if write_hr:
        f = open(wann_inp, 'a')
        f.write("write_hr = true" + "\n")
        f.close()
    else:
        f = open(wann_inp, 'w')
        f.write("num_wann = " + str(numwan) + "\n")
        if numwan < nbands:
            f.write("exclude_bands " + str(numwan + 1) + "-" + str(nbands) + "\n")
        if numwan == nbands - 1:
            f.write("exclude_bands " + str(numwan + 1) + "\n")
        f.write("bands_plot = true" + "\n")
        f.write("begin kpoint_path" + "\n")
        for i in range(len(labels)):
            for j in range(len(labels[i])-1):
                label1=labels[i][j]
                label2=labels[i][j+1]
                string = write_kpath_lines(kpts[label1], label1, kpts[label2], label2)
                f.write(string)
        f.write("end kpoint_path" + "\n")
        f.write("bands_num_points 40" + "\n")
        f.write("bands_plot_format gnuplot" + "\n")
        f.write("Begin Projections" + "\n")
        f.write("random" + "\n")
        for element in elements:
            f.write(element[0] + ':')
            for l in range(4):
                if element[l + 1] > 0:
                    if l == 0:
                        f.write('l=' + str(l))
                    else:
                        f.write(';' + 'l=' + str(l))
            f.write('\n')
        f.write("End Projections" + "\n")
        f.write("num_iter = 500" + "\n")
        f.close()


def write_kpath_lines(kpt1, label1, kpt2, label2):
    """
    Your Comments Here
    """
    string = label1 + "  " + "%10.5f" % kpt1[0] + "%10.5f" % kpt1[1] + "%10.5f" % kpt1[
        2] + "  " + label2 + "  " + "%10.5f" % kpt2[0] + "%10.5f" % kpt2[1] + "%10.5f" % kpt2[2] + "\n"
    return string


def kpath_finder(filename):
    """
    Your Comments Here
    """
    struct = Structure.from_file(filename)
    spcg = SpacegroupAnalyzer(struct)
    pstruct = spcg.get_primitive_standard_structure(international_monoclinic=False)
    hskp = HighSymmKpath(pstruct)
    kpts = hskp.kpath["kpoints"]
    labels = hskp.kpath["path"]
    return labels, kpts


def read_potcar(potcarfile, poscarfile):
    """
    Your Comments Here
    """
    struct = Structure.from_file(poscarfile)
    potcar = Potcar.from_file(potcarfile)
    sites = struct.species
    specs = potcar.spec
    wan_info = []
    for spec in specs:
        single = PotcarSingle.from_symbol_and_functional(symbol=spec['symbol'], functional='PBE_54')
        orbitals = single.electron_configuration
        element = single.element
        nat = str(sites).count(single.element)
        spdf = []
        for orbital in orbitals:
            spdf.append(orbital[1])
        scnt = spdf.count('s')
        pcnt = spdf.count('p')
        dcnt = spdf.count('d')
        fcnt = spdf.count('f')
        numwan = nat * (scnt * 1 + pcnt * 3 + dcnt * 5 + fcnt * 7)
        wan_info.append([element, scnt, pcnt, dcnt, fcnt, numwan])
    return wan_info


def read_vbm(fname):
    """
    Your Comments Here
    """
    vasprun = Vasprun(fname, parse_potcar_file=False)
    gap, cbm, vbm, is_direct = vasprun.eigenvalue_band_properties
    return gap, vbm


def read_wannier(fname_band, fname_kpt, vbm):
    """
    Your Comments Here
    """
    f = open(fname_band)
    contents = f.readlines()
    f.close()

    kvals = []
    evals = []

    for content in contents:
        lines = str.split(content)
        if len(lines) != 0:
            kvals.append(eval(lines[0]))
            evals.append(eval(lines[1]))
    with open(fname_kpt) as f:
        first_line = f.readline()
    nkpt = eval(str.split(first_line)[0])
    nband = int(len(evals) / nkpt)
    kpts_flat = np.array(kvals)
    eigs_flat = np.array(evals) - vbm
    kpts = kpts_flat.reshape(nband, nkpt)
    eigs = eigs_flat.reshape(nband, nkpt)

    return kpts, eigs


def read_vasp(fname, vbm):
    """
    Your Comments Here
    """
    vasprun = Vasprun(fname)
    eigs = vasprun.eigenvalues
    nkpt = len(vasprun.actual_kpoints)
    kptvasp = vasprun.actual_kpoints
    out = []
    for spin, d in eigs.items():
        for k, val in enumerate(d):
            for n, eigenval in enumerate(val):
                out.append(eigenval[0] - vbm)
    nbands = int(len(out) / nkpt)
    eigvasp = np.reshape(out, (nkpt, nbands))

    return kptvasp, eigvasp


def read_special_kpts(fname):
    """                                                                         
    Your Comments Here                                                          
    """
    f = open(fname)
    contents = f.readlines()
    f.close()
    spkptl = []
    spkptc = []
    for content in contents:
        if 'set xtics' in content:
            lines = str.split(content,sep='"')
            for n, words in enumerate(lines):
                if n!=0 and n%2!=0:
                    spkptl.append(words)
                elif n!=0:
                    spkptc.append(eval(words[:9]))

    return spkptl,spkptc

