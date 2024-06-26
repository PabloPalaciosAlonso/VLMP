import sys, os

import logging

from . import simulationStepBase

class meanMagnetizationMeasurement(simulationStepBase):
    """
    Component name: meanMagnetizationMeasurement
    Component type: simulationStep

    Author: P. Palacios Alonso
    Date: 18/10/2023

    This component writes the mean magnetization of the system to a file.

    :param outputFilePath: Path to the output file
    :type outputFilePath: str

    """

    def __init__(self,name,**params):
        super().__init__(_type = self.__class__.__name__,
                         _name = name,
                         availableParameters = {"outputFilePath",
                                                "startStep"},
                         requiredParameters  = {"outputFilePath"},
                         availableSelections = {"selection"},
                         requiredSelections  = set(),
                         **params)

        ############################################################
        ############################################################
        ############################################################

        parameters = {}

        parameters["outputFilePath"] = params["outputFilePath"]
        parameters["startStep"]      = params.get("startStep", 0)
        
        simulationStep = {
            name:{
              "type":["MagneticMeasure","MeasureMeanMagnetization"],
              "parameters":{**parameters}
            }
        }

        ############################################################

        self.setGroup("selection")
        self.setSimulationStep(simulationStep)



