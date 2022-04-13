import numpy as np
import ROOT


class JunoReactorFlux(object):

    def __init__(self) -> None :
        self.baseline  = 53000    # m
        self.power     = 2.9e9    # W
        self.time      = 6 * 365 * 24 * 60 * 60     # seconds in 6 years

        
        self.MO = 1

        # with DayaBay defect
        self.fU235  = 0.58
        self.fU238  = 0.07
        self.fPu239 = 0.3
        self.fPu241 = 0.05

        self.EU235  = 202.36  # MeV +- 0.26
        self.EU238  = 205.99  # MeV +- 0.52
        self.EPu239 = 211.12  # MeV +- 0.34
        self.EPu241 = 214.26  # MeV +- 0.33

        self.hU235_x,  self.hU235_y  = [], []
        self.hU238_x,  self.hU238_y  = [], []
        self.hPu239_x, self.hPu239_y = [], []
        self.hPu241_x, self.hPu241_y = [], []
        self.sin2theta12 = 0.307      # +- 0.013
        self.Deltam212 = 7.53e-5      # +- 0.18  eV2
        self.Deltam322IO = -2.546e-3  # +0.034 -0.040 eV2
        self.Deltam322NO = 2.453e-3   # +- 0.034 eV2
        self.sin2theta13 = 2.18e-2    # +- 0.07



    def SetBaseline(self, bl):
        self.baseline = bl

    def SetPower(self, p):
        self.power = p

    def SetMO(self, mo):
        self.MO = mo


    def CommonInput(self):
        ff = ROOT.TFile("JUNOInputs2022_01_06.root", "read")
        hU235  = ff.Get("HuberMuellerFlux_U235")
        nbin = hU235.GetNbinsX()
        for i in range(nbin):
            self.hU235_x.append(hU235.GetBinCenter(i+1))
            self.hU235_y.append(hU235.GetBinContent(i+1))
        hU238  = ff.Get("HuberMuellerFlux_U238")
        nbin = hU238.GetNbinsX()
        for i in range(nbin):
            self.hU238_x.append(hU238.GetBinCenter(i+1))
            self.hU238_y.append(hU238.GetBinContent(i+1))
        hPu239 = ff.Get("HuberMuellerFlux_Pu239")
        nbin = hPu239.GetNbinsX()
        for i in range(nbin):
            self.hPu239_x.append(hPu239.GetBinCenter(i+1))
            self.hPu239_y.append(hPu239.GetBinContent(i+1))
        hPu241 = ff.Get("HuberMuellerFlux_Pu241")
        nbin = hPu241.GetNbinsX()
        for i in range(nbin):
            self.hPu241_x.append(hPu241.GetBinCenter(i+1))
            self.hPu241_y.append(hPu241.GetBinContent(i+1))


        del hU235, hU238, hPu239, hPu241, ff

    
    def reactorFlux(self, Enu):
        JtoMeV = 6.242e12  # J to MeV
        phi = self.power * JtoMeV * self.time / (self.fU235*self.EU235 + self.fU238 * self.EU238 + self.fPu239 * self.EPu239 + self.fPu241 * self.EPu241) * (self.fU235 * np.interp(Enu, self.hU235_x, self.hU235_y) + self.fU238 * np.interp(Enu, self.hU238_x, self.hU238_y) + self.fPu239 * np.interp(Enu, self.hPu239_x, self.hPu239_y) + self.fPu241 * np.interp(Enu, self.hPu241_x, self.hPu241_y) )

        return phi      # unit: MeV^{-1}
        


    def survivalProb(self, Enu):
        if self.MO == 1:   ## Normal Ordering
            fast = 4*self.sin2theta13*(1-self.sin2theta13)*((1-self.sin2theta12)*np.sin((self.Deltam212 + self.Deltam322NO) * 1.27 * self.baseline /Enu)**2 + self.sin2theta12 * np.sin(1.27*self.Deltam322NO*self.baseline/Enu)**2 )
            slow = (1 - self.sin2theta13)**2 * 4*self.sin2theta12*(1-self.sin2theta12) * np.sin(1.27*self.Deltam212*self.baseline/Enu)**2

        elif self.MO == 2:   ## Inverted Ordering
            fast = 4*self.sin2theta13*(1-self.sin2theta13)*((1-self.sin2theta12)*np.sin((self.Deltam212 + self.Deltam322IO) * 1.27 * self.baseline /Enu)**2 + self.sin2theta12 * np.sin(1.27*self.Deltam322IO*self.baseline/Enu)**2 )
            slow = (1 - self.sin2theta13)**2 * 4*self.sin2theta12*(1-self.sin2theta12) * np.sin(1.27*self.Deltam212*self.baseline/Enu)**2

        else:
            print("Wrong Mass Ordering Input!!!!!!!!")

        return 1 - fast - slow


    
    def arrivedFlux(self, Enu):
        mtocm = 100
        return self.survivalProb(Enu) / (4*np.pi*self.baseline**2*mtocm**2) * self.reactorFlux(Enu)    # unit: MeV^-1 * cm^-2















    





