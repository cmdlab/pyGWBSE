# coding: utf-8

"""
This module defines tasks that support running vasp in various ways.
"""

import os
import subprocess

from atomate.utils.utils import env_chk, get_logger
from fireworks import explicit_serialize, FiretaskBase

__author__ = 'Anubhav Jain <ajain@lbl.gov>'
__credits__ = 'Shyue Ping Ong <ong.sp>'

logger = get_logger(__name__)


@explicit_serialize
class Run_Vasp(FiretaskBase):
    """
    Execute a command directly (no custodian).

    Args:
        cmd (str): the name of the full executable to run. Supports env_chk.

    Other Parameters:
        expand_vars (str): Set to true to expand variable names in the cmd.
    """

    required_params = ["vasp_cmd"]
    optional_params = ["expand_vars"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        cmd = env_chk(self["vasp_cmd"], fw_spec)
        if self.get("expand_vars", False):
            cmd = os.path.expandvars(cmd)

        logger.info("Running command: {}".format(cmd))
        return_code = subprocess.call(cmd, shell=True)
        logger.info("Command {} finished running with returncode: {}".format(cmd, return_code))


@explicit_serialize
class Run_Sumo(FiretaskBase):
    """
    Your Comments Here
    """
    required_params = ["sumo_cmd"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        cmd = env_chk(self["sumo_cmd"], fw_spec)
        logger.info("Running command: {}".format(cmd))
        return_code = subprocess.call(cmd, shell=True)
        logger.info("Command {} finished running with returncode: {}".format(cmd, return_code))


@explicit_serialize
class Run_Wannier(FiretaskBase):
    """
    Your Comments Here
    """
    required_params = ["wannier_cmd"]

    def run_task(self, fw_spec):
        """
        Your Comments Here
        """
        cmd = env_chk(self["wannier_cmd"], fw_spec)
        logger.info("Running command: {}".format(cmd))
        return_code = subprocess.call(cmd, shell=True)
        logger.info("Command {} finished running with returncode: {}".format(cmd, return_code))
