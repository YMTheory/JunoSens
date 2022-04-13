import numpy as np
from JunoDetector import JunoDetector
from JunoReactorFlux import JunoReactorFlux
import scipy.integrate as integrate

class JunoIBDSignal(object):

    def __init__(self) -> None :

        # binning methods
        self.Ee_edge = []
        self.Ee_edge.append(0.80)
        self.Ee_edge.append(0.94)
        self.Ee_edge.append(7.44)
        self.Ee_edge.append(7.80)
        self.Ee_edge.append(8.20)
        self.Ee_edge.append(12.0)
        self.nbin = []
        self.nbin.append(1)
        self.nbin.append(325)
        self.nbin.append(9)
        self.nbin.append(4)
        self.nbin.append(1)
        self.width = []
        self.width.append( (self.Ee_edge[1] - self.Ee_edge[0]) / self.nbin[0] )
        self.width.append( (self.Ee_edge[2] - self.Ee_edge[1]) / self.nbin[1] )
        self.width.append( (self.Ee_edge[3] - self.Ee_edge[2]) / self.nbin[2] )
        self.width.append( (self.Ee_edge[4] - self.Ee_edge[3]) / self.nbin[3] )
        self.width.append( (self.Ee_edge[5] - self.Ee_edge[4]) / self.nbin[4] )

        self.bin_edge = []
        for i in range(len(self.width)):
            for j in range(self.nbin[i]):
                self.bin_edge.append(self.Ee_edge[i] + self.width[i]*j)
        self.bin_edge.append(self.Ee_edge[-1]) 


        ### Initialization and dataloading...
        self.reactor = JunoReactorFlux()
        self.reactor.CommonInput()

        self.det = JunoDetector()
        self.det.CommonInput()



    def BinnedNe(self, i):

        eff = self.det.GetEfficiency()
        Np  = self.det.GetNproton()

        Ep_low, Ep_high = self.bin_edge[i], self.bin_edge[i+1]
        print("Bin Range [%.5f, %.5f]" %(Ep_low, Ep_high) )
        
        f = lambda Ep, ct, Enu : self.reactor.arrivedFlux(Enu) * self.det.IBDdiffXsec(Enu, ct) * Np * eff * self.det.detector_response(Enu, ct, Ep)

        Ti = integrate.tplquad(f, 1.8, 15, lambda Enu: -1, lambda Enu: 1, lambda Enu, ct: Ep_low, lambda Enu, ct: Ep_high)
        print("Predicted Bin Count = ", Ti)
        return Ti



    def PredPromptSignalSpec(self):
        T = []
        for i in range(len(self.bin_edge)-1) :
            T.append(self.BinnedNe(i))

        print(T)
        return T

















