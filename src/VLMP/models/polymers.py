import sys, os

from VLMP.models import modelBase

import numpy as np

class wlc(modelBase):
    """
    This model produces a WLC with N particles.

    Parameters:
        N: Number of particles
        mass: Mass of the particles
        b: Distance between two consecutive particles
        Kb: Spring constant for bonds
        Ka: Spring constant for angles

    Selectors:
        particleIndex: Select a particle by its index
    """

    availableParameters = ["N","mass","b","Kb","Ka"]
    availableSelectors  = ["particleIndex"]

    def __init__(self,**kwargs):
        super().__init__()
        self.logger.debug("[WLC] Using model Worm Like Chain (WLC)")

        ############################################################
        ####################  Model Parameters  ####################
        ############################################################

        # Check if all parameters given by kwargs are available
        for key in kwargs:
            if key not in self.availableParameters:
                self.logger.error(f"[WLC] Parameter \"{key}\" not available")
                sys.exit(1)

        #Read the number of particles
        self.N = kwargs.get("N",0)

        #Read the particle mass and the distance between two consecutive particles
        self.mass = kwargs.get("mass",1.0)
        self.b    = kwargs.get("b",1.0)

        #Read values for the spring constant for bonds and angles
        self.Kb = kwargs.get("Kb",1.0)
        self.Ka = kwargs.get("Ka",1.0)

        self.logger.debug(f"[WLC] Generating WLC with \"N\":{self.N} particles")
        self.logger.debug(f"[WLC] Parameters \"mass\":{self.mass}, \"b\":{self.b}, \"Kb\":{self.Kb}, \"Ka\":{self.Ka}")

        ########################################################
        ############### Generate the WLC model #################
        ########################################################

        #Generate coordinates, a line along the z axis
        self.coords = np.zeros((self.N,3))
        for i in range(self.N):
            self.coords[i,2] = i*self.b

        #Generate topology
        self.topology = {}
        #Generate particle types
        self.topology["particleTypes"] = {}
        self.topology["particleTypes"]["labels"] = ["name", "mass", "radius", "charge"]
        self.topology["particleTypes"]["data"] = [["A", self.mass, 0.5*self.b, 0.0]]

        #Generate structure
        self.topology["structure"] = {}
        self.topology["structure"]["labels"] = ["id","type","modelId"]

        self.topology["structure"]["data"] = []
        for i in range(self.N):
            self.topology["structure"]["data"].append([i,"A",0])

        #Generate force field
        self.topology["forceField"] = {}

        #Generate bonds
        self.topology["forceField"]["bonds"] = {}
        self.topology["forceField"]["bonds"]["type"]       = ["Bond2","HarmonicConst_K_r0"]
        self.topology["forceField"]["bonds"]["parameters"] = {"K":self.Kb, "r0":self.b}
        self.topology["forceField"]["bonds"]["labels"]     = ["id_i","id_j"]
        self.topology["forceField"]["bonds"]["data"]       = []

        for i in range(self.N-1):
            self.topology["forceField"]["bonds"]["data"].append([i,i+1])

        #Generate angles
        self.topology["forceField"]["angles"] = {}
        self.topology["forceField"]["angles"]["type"]       = ["Bond3","KratkyPorodConst_K"]
        self.topology["forceField"]["angles"]["parameters"] = {"K":self.Ka}
        self.topology["forceField"]["angles"]["labels"]     = ["id_i","id_j","id_k"]
        self.topology["forceField"]["angles"]["data"]       = []

        for i in range(self.N-2):
            self.topology["forceField"]["angles"]["data"].append([i,i+1,i+2])

    def selection(self,**kwargs):
        # Check if all selectors given by kwargs are available
        for key in kwargs:
            if key not in self.availableSelectors:
                self.logger.error(f"[WLC] Selector \"{key}\" not available")
                sys.exit(1)

        return kwargs.get("particleIndex",[])

    def write(self, filePath:str):
        pass
