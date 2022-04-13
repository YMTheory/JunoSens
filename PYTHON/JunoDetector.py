import numpy as np
import ROOT



class JunoDetector(object):

    def __init__(self) -> None:


        # IBD differential cross section
        self.fIBDdiffXsec = ROOT.TF2()
        self.fEpositron   = ROOT.TF2()

        self.IBD_x = []
        self.IBD_y = []

        # detector configuration
        self.Nproton = 1.43512e33
        self.efficiency = 0.822
        self.hNL_x = []
        self.hNL_y = []
        self.a = 0.0261
        self.b = 0.0082
        self.c = 0.0123

    
    def GetNproton(self):
        return self.Nproton


    def GetEfficiency(self):
        return self.efficiency



    def CommonInput(self):
        ff = ROOT.TFile("JUNOInputs2022_01_06.root", "read")
        self.fIBDdiffXsec = ff.Get("dsigma_dcos_Enu_cos_DYB")
        self.fEpositron = ff.Get("Epositron_Enu_cos_DYB")
        hNL = ff.Get("positronScintNL")
        for i in range(hNL.GetNbinsX()):
            self.hNL_x.append(hNL.GetBinCenter(i+1))
            self.hNL_y.append(hNL.GetBinContent(i+1))
        hIBDXsec = ff.Get("IBDXsec_VogelBeacom_DYB")
        for i in range(hIBDXsec.GetNbinsX()):
            self.IBD_x.append(hIBDXsec.GetBinCenter(i+1))
            self.IBD_y.append(hIBDXsec.GetBinContent(i+1))


    
    def IBDdiffXsec(self, Enu, costheta):
        return self.fIBDdiffXsec.Eval(Enu, costheta)


    def IBDtotXsec(self, Enu):
        return np.interp(Enu, self.IBD_x, self.IBD_y)

    
    def nonlinearity(self, Ee):
        return np.interp(Ee+1.022, self.hNL_x, self.hNL_y) * Ee
    


    def resolution(self, Evis):
        return np.sqrt(self.a**2/Evis + self.b**2 + self.c**2/Evis**2) * Evis


    
    def detector_response(self, Enu, ct, Ep):
        Ee = self.fEpositron(Enu, ct) 
        Edep = Ee + 1.022
        Evis = self.nonlinearity(Edep)
        sigma = self.resolution(Evis)
        #print("neutrino energy", Enu)
        #print("costheta", ct)
        #print("positron KE", Ee)
        #print("positron Edep",  Edep)
        #print("positron Evis", Evis)
        #print("positron Esigma", sigma)
        #print("calc prompt E", Ep)

        return 1 / (np.sqrt(2*np.pi)*sigma) * np.exp(-(Ep - Evis)**2/(2*sigma**2)) 
    






