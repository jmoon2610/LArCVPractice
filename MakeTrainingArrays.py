import ROOT
from ROOT import *
import numpy
from sys import argv
import matplotlib.pyplot as plt

def IsContained(x,y,z,edge_x,edge_y,edge_z):

    # --------------------------------------------------------------------------------------- #                                                       
    # Function uses detector active volume coordinates in conjuction with an input parameter                                                            
    # for distance from edge of active volume to define a fiducial volume and return 0/1 if                                                             
    # the given coordinates are contained in the so defined fiducial volume                                                                          
    # --------------------------------------------------------------------------------------- #                                                          

    #--Defines active volume coordinate boundaries--#                                                                                              
    xmax =  256.25
    xmin =  0
    ymax =  116.5
    ymin = -116.5
    zmax =  1036.8
    zmin =  0
    #-----------------------------------------------#                                                                                                      

    if x < (xmax - edge_x) and x > (xmin + edge_x) and y < (ymax - edge_y) and y > (ymin + edge_y) and z < (zmax - edge_z) and z > (zmin + edge_z):
        return 1
    else:
        return 0

def SetupLarliteManager():
    
    global manager

    manager = larlite.storage_manager()
    manager.set_io_mode(manager.kREAD)
    for infile in argv[1:]:
        manager.add_in_filename(infile)
    manager.set_in_rootdir("")
    manager.open()

#--------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------#
#--------------------------------------------------------------------------------------------------------------#

MC         = True
CCQE       = False
Contained  = False
nu_E_cut   = False

nu_E_max        = 0.5       #GeV
nu_E_min        = 0   
WindowLength    = 130   
CoincWindow     = 6
WinStart        = 220 if MC == True else 210       #220 for MC , 210 for Data
WinEnd          = 350 if MC == True else 340       #350 for MC , 340 for Data
TickSize        = 0.015625                         #15.625 ns per tick
MaxEvents       = 10

SetupLarliteManager()

count = 0
while manager.next_event():

    if count > MaxEvents: continue

    # ---------- Grab vectors containing particle info --------- #
    ophits   = manager.get_data(larlite.data.kOpHit,"ophit")
    if MC == True:
        mcdata   = manager.get_data(larlite.data.kMCTrack,"mcreco")                                                                                                
        mc       = manager.get_data(larlite.data.kMCTruth,"generator")                                                                                             
    #shower   = manager.get_data(larlite.data.kMCShower,"mcreco")
    # ---------------------------------------------------------- #

    if mc == True:

        if mcdata.size() == 0: continue
    
        nu_trajectory     = mc[0].GetParticles()[0].Trajectory()
        nu_E              = nu_trajectory[0].E()
        nu_xf             = nu_trajectory[nu_trajectory.size()-1].X()
        nu_yf             = nu_trajectory[nu_trajectory.size()-1].Y()
        nu_zf             = nu_trajectory[nu_trajectory.size()-1].Z()
        interaction_code  = mc[0].GetNeutrino().InteractionType()
        scndry_PDG        = mcdata[0].PdgCode() 
        scndry_x0         = mcdata[0].Start().X()
        scndry_y0         = mcdata[0].Start().Y()
        scndry_z0         = mcdata[0].Start().Z()
        scndry_xf         = mcdata[0].End().X()
        scndry_yf         = mcdata[0].End().Y()
        scndry_zf         = mcdata[0].End().Z()

        # ----------- Performs some data filters on CCQE events, containment, initial neutrino energy ------------- #
        if scndry_PDG not in [11,-11,13,-13]: continue

        if CCQE == True and interaction_code != 1001:
            continue

        if Contained == True and (IsContained(nu_xf,nu_yf,nu_zf,5,5,5) == 0 or IsContained(scndry_x0,scndry_y0,scndry_z0,5,5,5) == 0 or IsContained(scndry_xf,scndry_yf,scndry_zf,0,0,0) == 0):
            continue

        if nu_E_cut  == True and (nu_E > nu_E_max or nu_E < nu_E_min):
            continue
        # --------------------------------------------------------------------------------------------------------- #

    clusters = [[0]*(WindowLength/CoincWindow + 1) for _ in xrange(32)]
    
    for ophit in ophits:
        if ophit.PeakTime()/TickSize > WinEnd or ophit.PeakTime()/TickSize < WinStart:
            continue
        clusters[ophit.OpChannel()][int((ophit.PeakTime()/TickSize - WinStart)/CoincWindow)]+=ophit.PE()

    plt.clf()
    heatmap = plt.pcolor(numpy.asarray(clusters))
    pltname = "heatmap%i"%(count)
    plt.savefig(pltname)

    count+=1
