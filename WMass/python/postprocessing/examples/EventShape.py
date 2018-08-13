import ROOT
import os, array
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

#sci-kit for thrust calculation
from sklearn.decomposition import PCA

from CMGTools.WMass.postprocessing.framework.datamodel import Collection 
from CMGTools.WMass.postprocessing.framework.eventloop import Module
#include the event shape from TOP analysis
from CMGTools.WMass.postprocessing.framework.eventShapeTools import EventShapeTool

from math import *

class EventShapeProducer(Module):
    def __init__(self):
        self.vars = ("pt","eta","phi","mass","pdgid")
        self.genwvars = ("charge","pt","eta","phi","mass","mt","y","costcs","cost2d","phics","costcm","decayId")

        #added by me
        if "EventShapes_cc.so" not in ROOT.gSystem.GetLibraries():
            print "Compile Event shape script"
            ROOT.gSystem.Load("/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/fastjet/3.1.0/lib/libfastjet.so")
            ROOT.gSystem.AddIncludePath("-I/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/fastjet/3.1.0/include")
            ROOT.gSystem.Load("/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/fastjet-contrib/1.020/lib/libfastjetcontribfragile.so")
            ROOT.gSystem.AddIncludePath("-I/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/fastjet-contrib/1.020/include")
            ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/postprocessing/helpers/EventShapes.cc+" % os.environ['CMSSW_BASE'])
        else:
            print "EventShapes_cc.so found in ROOT libraries"
        #

        self._worker = ROOT.EventShapes()
        self._workertool = EventShapeTool()

    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree) # initReaders must be called in beginFile
        self.out = wrappedOutputTree
        #include your branches of the Ntuples
        self.out.branch("tau1", "F")
        self.out.branch("tau2", "F")
        self.out.branch("tau3", "F")
        self.out.branch("tau4", "F")

        self.out.branch("sphericity", "F")
        self.out.branch("aplanarity", "F")
        self.out.branch("C", "F")
        self.out.branch("D", "F")
        self.out.branch("detST", "F")

        #thrust
        self.out.branch("thrustMinor", "F")
        self.out.branch("thrustMajor", "F")
        self.out.branch("thrust", "F")
        self.out.branch("oblateness", "F")

        self.out.branch("thrustTransverse", "F")
        self.out.branch("thrustTransverseMinor", "F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def initReaders(self,tree): # this function gets the pointers to Value and ArrayReaders and sets them in the C++ worker class
        try:
            self.nChgPart = tree.valueReader("nChgPart")
            for B in ("pt","eta","phi","mass","pdgid") : setattr(self,"ChgPart_"+B, tree.arrayReader("ChgPart_"+B))

            self._worker.setChgParticles(self.nChgPart,self.ChgPart_pt,self.ChgPart_eta,self.ChgPart_phi,self.ChgPart_mass, self.ChgPart_pdgid) 

        except:
            print '[EventShapeProducer][Warning] Unable to attach to charged particles'


        try:
            for B in ("pt","eta","phi") : setattr(self,"LepGood_"+B, tree.arrayReader("LepGood_"+B))

            self._worker.setGoodLepton(self.LepGood_pt,self.LepGood_eta,self.LepGood_phi)

        except:
            print '[EventShapeProducer][Warning] Unable to attach to good leptons'

        self._ttreereaderversion = tree._ttreereaderversion # self._ttreereaderversion must be set AFTER all calls to tree.valueReader or tree.arrayReader

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        if event._tree._ttreereaderversion > self._ttreereaderversion: # do this check at every event, as other modules might have read further branches
            self.initReaders(event._tree)

        # do NOT access other branches in python between the check/call to initReaders and the call to C++ worker code
        ## Subjettiness
        self._worker.runNjettiness()

        #define the charged particles
        chgparticlesForEvShape = self._worker.buildchgparticles()
        
        #these can be used also for the PCA for the thrust
        #first define the three vectors needed for the thrust
        Vect3chg = []

        for i in range(0,chgparticlesForEvShape.size()):
            Vect3chg.append(chgparticlesForEvShape[i].Vect())

        #then, you can define and fit the PCA
        pca3d = PCA(n_components=3)
        pca3d.fit(Vect3chg)
        
        #This gives you already the axes along which the thrust can be calculated, through pca.components
        thrust_1=0
        thrust_2=0
        thrust_3=0
        momentum_sum = 0

        for i in range(0,chgparticlesForEvShape.size()):
            thrust_1 += abs(np.dot(Vect3chg[i],pca3d.components_[0]))
            thrust_2 += abs(np.dot(Vect3chg[i],pca3d.components_[1]))
            thrust_3 += abs(np.dot(Vect3chg[i],pca3d.components_[2]))
            momentum_sum += sqrt(pow(Vect3chg[i].Px(),2)+pow(Vect3chg[i].Py(),2)+pow(Vect3chg[i].Pz(),2))

        thrust_1/=momentum_sum
        thrust_2/=momentum_sum
        thrust_3/=momentum_sum

        Thrust3D = sorted([thrust_1,thrust_2,thrust_3])

        #Calculating the transverse thrust
        Vect2chg = []

        for i in range(0,chgparticlesForEvShape.size()):
            SingleChg2d = [chgparticlesForEvShape[i].Px(),chgparticlesForEvShape[i].Py()]
            Vect2chg.append(SingleChg2d)

        #then, you can define and fit the PCA
        pca2d = PCA(n_components=2)
        pca2d.fit(Vect2chg)
        
        #This gives you already the axes along which the thrust can be calculated, through pca.components
        thrust_1=0
        thrust_2=0
        momentum_sum = 0

        for i in range(0,chgparticlesForEvShape.size()):
            thrust_1 += abs(np.dot(Vect2chg[i],pca2d.components_[0]))
            thrust_2 += abs(np.dot(Vect2chg[i],pca2d.components_[1]))
            momentum_sum += sqrt(pow(Vect2chg[i][0],2)+pow(Vect2chg[i][1],2))
        
        thrust_1/=momentum_sum
        thrust_2/=momentum_sum

        Thrust2D = sorted([thrust_1,thrust_2])

        #other event shapes

        self._workertool.analyseNewEvent(chgparticlesForEvShape)

        #from event shape tool
        self._workertool.computeEventShapes()
        sphericity = self._workertool.sphericity
        aplanarity = self._workertool.aplanarity
        C = self._workertool.C
        D = self._workertool.D
        detST = self._workertool.detST

        #this is for the njettiness coming from the helper
        tau1 = self._worker.tau1()
        tau2 = self._worker.tau2()
        tau3 = self._worker.tau3()
        tau4 = self._worker.tau4()
        
        self.out.fillBranch("tau1", tau1)
        self.out.fillBranch("tau2", tau2)
        self.out.fillBranch("tau3", tau3)
        self.out.fillBranch("tau4", tau4)

        #event shapes
        self.out.fillBranch("sphericity", sphericity)
        self.out.fillBranch("aplanarity", aplanarity)
        self.out.fillBranch("C", C)
        self.out.fillBranch("D", D)
        self.out.fillBranch("detST", detST)

        self.out.fillBranch("thrust", Thrust3D[2])
        self.out.fillBranch("thrustMajor", Thrust3D[1])
        self.out.fillBranch("thrustMinor", Thrust3D[0])
        self.out.fillBranch("oblateness", Thrust3D[1]-Thrust3D[0])

        self.out.fillBranch("thrustTransverse", Thrust2D[1])
        self.out.fillBranch("thrustTransverseMinor", Thrust2D[0])

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

EventShape = lambda : EventShapeProducer()
