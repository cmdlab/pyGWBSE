# coding: utf-8

from fireworks import FiretaskBase, FWAction, explicit_serialize

from pymongo import MongoClient


def connect_database(username, password, databaseName):
    """
    Your Comments Here
    """
    # MongoDB connection info
    hostname = '209.147.162.7'
    port = 27017
    client = MongoClient(hostname, port)
    db = client[databaseName]
    db.authenticate(username, password)

    return db


def get_nbands(mid):
    """
    Your Comments Here
    """
    db = connect_database('tatha', 'tatha', 'tatha_db')
    found = False
    eqpcollection = db.get_collection('QP_Results')
    for x in eqpcollection.find({"material_id": mid, "task_label": {"$regex": 'Conv'}, 'ifconv': True}):
        name = x["formula_pretty"]
        nbands = x["incar"]["NBANDS"]
        found = True
    if found == True:
        return nbands


@explicit_serialize
class getfromdb(FiretaskBase):
    """
    Insert quasi-particle energies into the database for a GW calculation.
    """
    required_params = ["mat_name"]

    def run_task(self, fw_spec):
        mat_name = self['mat_name']
        nbands = get_nbands(mat_name)
        return FWAction(update_spec={"conv_nbands": nbands})
