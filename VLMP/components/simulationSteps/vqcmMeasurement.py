import sys, os

import logging

from . import simulationStepBase

class vqcmMeasurement(simulationStepBase):
    """
    Component name: vqcmMeasurement
    Component type: simulationStep

    Author: Pablo Palacios-Alonso and Pablo Ibáñez-Freire
    Date: 2/11/2023

    """


    def __init__(self,name,**params):
        super().__init__(_type = self.__class__.__name__,
                         _name = name,
                         availableParameters = {"outputFilePath",
                                                "resonatorImpedance",
                                                "kernel",
                                                "tolerance",
                                                "f0",
                                                "overtone",
                                                "hydrodynamicRadius",
                                                "viscosity",
                                                "vwall",
                                                "fluidDensity",
                                                "maxNIterations",
                                                "toleranceConvergence",
                                                "maxNIterations",
                                                "memory",
                                                "damping",
                                                "notAcceleratedInterval",
                                                "h",
                                                "printSteps",
                                                "tetherInteractorNames"},
                         requiredParameters  = {"outputFilePath",
                                                "f0","overtone",
                                                "hydrodynamicRadius",
                                                "viscosity",
                                                "vwall",
                                                "fluidDensity"},
                         availableSelections = set(),
                         requiredSelections  = set(),
                         **params)

        ############################################################
        ############################################################
        ############################################################

        parameters = {}

        parameters["outputFilePath"] = params["outputFilePath"]

        parameters["kernel"]    = params.get("kernel", "Peskin3p")
        if parameters["kernel"] != "Gaussian" and parameters["kernel"] != "Peskin3p":
            self.logger.error("[vqcmMeasurement] The selected kernel is not valid. Valid kernels: Gaussian, Peskin3p")

        if "tolerance" in params and parameters["kernel"] == "Peskin3p":
            self.logger.error("[vqcmMeasurement] tolerance parameter is not available for Peskin3p kernel")
            
        if parameters["kernel"] == "Gaussian":
            parameters["tolerance"] = params.get("tolerance", 1e-5)
            
        parameters["f0"]                 = params["f0"]
        parameters["overtone"]           = params["overtone"]
        parameters["hydrodynamicRadius"] = params["hydrodynamicRadius"]
        parameters["viscosity"]          = params["viscosity"]
        parameters["vwall"]              = params["vwall"]
        parameters["fluidDensity"]       = params["fluidDensity"]
    
        parameters["maxNIterations"]         = params.get("maxNIterations", 10000)
        parameters["toleranceConvergence"]   = params.get("toleranceConvergence", 1e-4)
        parameters["memory"]                 = params.get("memory", 5)
        parameters["damping"]                = params.get("damping", 1e-5)
        parameters["notAcceleratedInterval"] = params.get("notAcceleratedInterval", 2)
        parameters["h"]                      = params.get("h", 0.0) # 0.0 means it is computed by UAMMD
        parameters["resonatorImpedance"]     = params.get("resonatorImpedance", -1.0) #-1.0 means UAMMD ignores it
        parameters["printSteps"]             = params.get("printSteps", 0)
        tetherInteractorNames                = params.get("tetherInteractorNames",[])
        if len(tetherInteractorNames)>0:
            parameters["tetherInteractorNames"] = tetherInteractorNames
        
        simulationStep = {
            name:{
              "type":["OscillatingFluidMeasure","VQCMMeasure"],
              "parameters":{**parameters}
            }
        }

        ############################################################

        self.setSimulationStep(simulationStep)



