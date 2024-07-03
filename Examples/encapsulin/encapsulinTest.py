import VLMP

Lx = 120
Ly = 120
Lz = 120

z0    = 40
z0Std = 0

KInter = 60
KIntra = 60

temperature = 1
blobRadius  = 3
blobMass    = 1

dt                  = 0.001
numberOfEncapsulins = 1

nsteps = 5000000
nsave  = 5000


simulationPool = [{"system":[{"type":"simulationName","parameters":{"simulationName":"testEncapsulin"}},
                             {"type":"backup","parameters":{"backupIntervalStep":100000}}],
                   "units":[{"type":"none"}],
                   "types":[{"type":"basic"}],
                   "ensemble":[{"type":"NVT","parameters":{"box":[Lx,Ly,Lz],"temperature":temperature}}],
                   "integrators":[{"type":"BBK","parameters":{"timeStep":dt,
                                                              "frictionConstant":1,
                                                              "integrationSteps":nsteps}}],
                   "models":[{"type":"TM_ENCAPSULIN_CG",
                              "parameters":{"numberOfEncapsulins":numberOfEncapsulins,
                                            "blobRadius":blobRadius,
                                            "blobMass":blobMass,
                                            "KInterPentamerFirst":KInter,
                                            "KIntraPentamerFirst":KIntra,
                                            "KInterPentamerSecond":KInter,
                                            "KIntraPentamerSecond":KIntra,
                                            "heightMean":z0,
                                            "heightStd":z0Std,
                                            "heightReference":-Lz/2.0}}],
                   "simulationSteps":[{"type":"saveState","parameters":{"intervalStep":nsave,
                                                                        "outputFilePath":"test",
                                                                        "outputFormat":"sp"}},
                                      {"type":"info","parameters":{"intervalStep":nsave*10}}]

                           }]


vlmp = VLMP.VLMP()

vlmp.loadSimulationPool(simulationPool)
vlmp.distributeSimulationPool("one")
vlmp.setUpSimulation("TEST")

