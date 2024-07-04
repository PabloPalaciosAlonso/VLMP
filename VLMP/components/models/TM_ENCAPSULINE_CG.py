from VLMP.components.models import modelBase

import os
import json
import numpy as np
from scipy.spatial import cKDTree

class TM_ENCAPSULIN_CG(modelBase):
    """
    Component name: TM_ENCAPSULIN_CG
    Component type: model

    Author: Pablo Palacios-Alonso
    Date: 01/07/2024

    Model coarse grained of termothoga marithima encapsulin, with one bead per monomer
    """


    def __applyPBC(self, positions, box_size):
        """
        Create periodic images of the system.
        
        Parameters:
        positions (np.ndarray): Array of positions of shape (N, D) where N is the number of points and D is the dimensionality.
        box_size (float): The size of the periodic box.
        
        Returns:
        np.ndarray: Array of replicated positions considering periodic boundary conditions.
        """
        shifts = np.array([-1, 0, 1])
        offsets = np.array(np.meshgrid(shifts, shifts, shifts)).T.reshape(-1, 3)
        
        replicated_positions = []
        for offset in offsets:
            replicated_positions.append(positions + offset * box_size)
            
        return np.vstack(replicated_positions)

    def __computeNewPosition(self, sphPositions,X,Y,Z, radius,
                             heightMean, heightStd, heightReference,
                             ntries):
    
        newPosition = []
        count = 0
        while count<ntries:
            if heightStd > 0.0:
                height = np.random.normal(heightMean,heightStd)
            else:
                height = heightMean
            height += heightReference
            if height > Z - radius or height < -Z + radius:
                continue

            x = np.random.uniform(-X,X)
            y = np.random.uniform(-Y,Y)
            center = [x,y,height]
        
            if len(sphPositions)>0:
                sphPositionsPBC = self.__applyPBC(sphPositions, 2*X)
                tree = cKDTree(sphPositionsPBC)
                minDst, minDstIndex = tree.query(center, k=1)
            else:
                minDst = np.inf

            if minDst > 2.0*radius*1.05:            
                newPosition = center
                break
            count+=1
            
        return newPosition


    def __init__(self,name,**params):
        super().__init__(_type = self.__class__.__name__,
                         _name= name,
                         availableParameters = {"particleName",
                                                "blobMass",
                                                "blobCharge",
                                                "blobRadius",
                                                "numberOfEncapsulins",
                                                "KInterPentamerFirst",
                                                "KIntraPentamerFirst",
                                                "KInterPentamerSecond",
                                                "KIntraPentamerSecond",
                                                "heightMean","heightStd",
                                                "heightReference",
                                                "maxTries",
                                                "encapsulinDistribution",
                                                "encapsulinCenters"},
                         requiredParameters  = {"KInterPentamerFirst",
                                                "KIntraPentamerFirst",
                                                "blobRadius",
                                                "heightMean",
                                                "numberOfEncapsulins"},
                         definedSelections   = {"particleId"},
                         **params)

        ############################################################

        encapsulinDistribution = params.get("encapsulinDistribution", "random")
        if encapsulinDistribution not in ["random", "fixed"]:
            self.logger.error("Encapsulin distribution type not recognized")
            raise Exception("Encapsulin distribution type not recognized")

        blobName = params.get("particleName","Encapsulin")

        numberOfEncapsulins = params["numberOfEncapsulins"]
        blobRadius          = params["blobRadius"]
        blobMass            = params.get("blobMass",0)
        blobCharge          = params.get("blobCharge",0)
        heightMean          = params["heightMean"]
        heightStd           = params.get("heightStd",0.0)        
        heightReference     = params.get("heightReference",0.0)

        KInterPentamerFirst = params["KInterPentamerFirst"]
        KIntraPentamerFirst = params["KIntraPentamerFirst"]
        
        KInterPentamerSecond = params.get("KInterPentamerSecond",0)
        KIntraPentamerSecond = params.get("KIntraPentamerSecond",0)

        self.maxTries    = params.get("maxTries",100)
        
        box = self.getEnsemble().getEnsembleComponent("box")

        X = box[0]/2.0
        Y = box[1]/2.0
        Z = box[2]/2.0
        
        #Check box height and generation height
        if heightReference > Z or heightReference < -Z:
            self.logger.error(f"Height reference is out of box {heightReference} > {Z} or {heightReference} < {-Z}")
            raise ValueError("Height reference is out of box")

        if heightMean + heightReference > Z or heightMean + heightReference < -Z:
            self.logger.error(f"Generation height is out of box {heightMean + heightReference} > {Z} or "
                              f"{heightMean + heightReference} < {-Z}")
            raise ValueError("Generation height is out of box")

        ############################################################

        
        dataPath = "data/TM_ENCAPSULIN.json"
        dataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),dataPath)
        
        with open(dataPath) as f:
            encapsulinInfo = json.load(f)

        encapsulinRadius = encapsulinInfo["meanRadius"]
        if encapsulinDistribution == "fixed":
            encapCenters = params["encapsulinCenters"]
            encapCenters = [np.array(center) for center in encapCenters]
        else:
            encapCenters = []        
            i              = 0
            tries          = 0
            while i < numberOfEncapsulins:
                newPosition = self.__computeNewPosition(encapCenters,X,Y,Z, encapsulinRadius,
                                                        heightMean, heightStd, heightReference,
                                                        self.maxTries*100*(i+1))
            
                if len(newPosition)>0:
                    encapCenters.append(newPosition)
                    i+=1
                else:
                    encapCenters = []
                    i            = 0
                    tries       += 1
    
                if tries >= self.maxTries:
                    print("Unable to find a correct configuration")
                    raise ValueError("The number of spheres is too high for the box size")

        # Generate encapsulins
        state = {}
        state["labels"]=["id","position"]
        state["data"]  =[]

        idCounter = 0
        encap2ids = []
        encap2pos = []
        
        # Extract positions from the data
        positions = [entry[1] for entry in encapsulinInfo["coordinates"]["data"]]

        # Convert to a numpy array
        positions = np.array(positions)
        
        for center in encapCenters:
            pos = center+positions
            ids = []
            for p in pos:
                state["data"].append([idCounter,list(p)])
                ids.append(idCounter)
                idCounter += 1
            encap2ids.append(ids.copy())
            encap2pos.append(pos.copy())


        
        types = self.getTypes()
        types.addType(name = blobName,
                      mass = blobMass,
                      radius = blobRadius,
                      charge = blobCharge)

        
        structure = {}
        structure["labels"] = ["id", "type", "modelId"]
        structure["data"]   = []
        
        for i in range(numberOfEncapsulins):
            for j in encap2ids[i]:
                structure["data"].append([j,blobName,i])


        bondsInterPentamerFirst = encapsulinInfo["bondsFirstNeighborsInterPentamer"]["data"]
        r0_interFirst           = encapsulinInfo["bondsFirstNeighborsInterPentamer"]["parameters"]["r0"]

        bondsIntraPentamerFirst = encapsulinInfo["bondsFirstNeighborsIntraPentamer"]["data"]
        r0_intraFirst           = encapsulinInfo["bondsFirstNeighborsIntraPentamer"]["parameters"]["r0"]

        bondsInterPentamerSecond = encapsulinInfo["bondsSecondNeighborsInterPentamer"]["data"]
        r0_interSecond           = encapsulinInfo["bondsSecondNeighborsInterPentamer"]["parameters"]["r0"]

        bondsIntraPentamerSecond = encapsulinInfo["bondsSecondNeighborsIntraPentamer"]["data"]
        r0_intraSecond           = encapsulinInfo["bondsSecondNeighborsIntraPentamer"]["parameters"]["r0"]
        
        forceField = {}
        forceField["BondPair"] = {}
        forceField["BondPair"]["parameters"] = {}
        forceField["BondPair"]["type"] = ["Bond2", "Harmonic"]
        forceField["BondPair"]["labels"] = ["id_i", "id_j", "K", "r0"]
        forceField["BondPair"]["data"] = []

        bonds = []
        # for i in range(numberOfEncapsulins):
        #     for j in range(59):
        #         for l in range(j+1, 60):
        #             pj = positions[j]
        #             pl = positions[l]
        #             dist = np.linalg.norm(pj-pl)
        #             bonds.append([encap2ids[i][j],
        #                           encap2ids[i][l],
        #                           KInterPentamer,
        #                           dist])
                

        for i in range(numberOfEncapsulins):
            for bond in bondsInterPentamerFirst:
                bonds.append([encap2ids[i][bond[0]],
                              encap2ids[i][bond[1]],
                              KInterPentamerFirst,
                              r0_interFirst])
            
            for bond in bondsIntraPentamerFirst:
                bonds.append([encap2ids[i][bond[0]],
                              encap2ids[i][bond[1]],
                              KIntraPentamerFirst,
                              r0_intraFirst])

            if KInterPentamerSecond > 0:
                for bond in bondsInterPentamerSecond:
                    bonds.append([encap2ids[i][bond[0]],
                                  encap2ids[i][bond[1]],
                                  KInterPentamerSecond,
                                  r0_interSecond])

            if KIntraPentamerSecond > 0:
                for bond in bondsIntraPentamerSecond:
                    bonds.append([encap2ids[i][bond[0]],
                                  encap2ids[i][bond[1]],
                                  KIntraPentamerSecond,
                                  r0_intraSecond])
        
        for i,j,k,r0 in bonds:
            forceField["BondPair"]["data"].append([i,j,k,r0])

            
        self.setState(state)
        self.setStructure(structure)
        self.setForceField(forceField)


    def processSelection(self,**params):

        sel = []
        if "particleId" in params:
            sel += params["particleId"]
        return sel

