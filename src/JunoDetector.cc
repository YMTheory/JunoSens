#include "JunoDetector.hh"
#include "JunoPullTerms.hh"

#include <TMath.h>
#include <TFile.h>

JunoDetector::JunoDetector()
{
    m_Nproton = 1.43512e33;
    m_effciency = 0.82;

    // Calibration Paper 
    m_a = 0.0262;
    m_b = 0.0082;
    m_c = 0.0123;

}


JunoDetector::~JunoDetector()
{
    delete hNLnominal;
    delete hNL0;
    delete hNL1;
    delete hNL2;
    delete hNL3;
    delete hIBDtotXsec;
    delete fIBDdiffXsec;
    delete fEpositron;
}


void JunoDetector::LoadCommonInputs()
{
    TFile* ff = new TFile("JUNOInputs2022_01_06.root", "read");
    fIBDdiffXsec = (TF2*)ff->Get("dsigma_dcos_Enu_cos_DYB");
    fEpositron   = (TF2*)ff->Get("Epositron_Enu_cos_DYB");
    hNLnominal   = (TH1D*)ff->Get("positronScintNL");
    hNL0         = (TH1D*)ff->Get("positronScintNLpull0");
    hNL1         = (TH1D*)ff->Get("positronScintNLpull1");
    hNL2         = (TH1D*)ff->Get("positronScintNLpull2");
    hNL3         = (TH1D*)ff->Get("positronScintNLpull3");
    hIBDtotXsec  = (TH1D*)ff->Get("IBDXsec_VogelBeacom_DYB");
}


double JunoDetector::IBDdiffXsec(double Enu, double costheta)
{
    return fIBDdiffXsec->Eval(Enu, costheta);
}


double JunoDetector::PositronEnergy(double Enu, double costheta)
{
    return fEpositron->Eval(Enu, costheta);
}

double JunoDetector::IBDtotXsec(double Enu)
{
    return hIBDtotXsec->Interpolate(Enu);
}



double JunoDetector::Nonlinearity(double Edep)
{
    //return Edep * hNLnominal->Interpolate(Edep);
    return Edep * ( hNLnominal->Interpolate(Edep) + JunoPullTerms::alpha_l0*(hNL0->Interpolate(Edep) - hNLnominal->Interpolate(Edep)) + JunoPullTerms::alpha_l1*(hNL1->Interpolate(Edep) - hNLnominal->Interpolate(Edep)) + JunoPullTerms::alpha_l2*(hNL2->Interpolate(Edep) - hNLnominal->Interpolate(Edep)) + JunoPullTerms::alpha_l3*(hNL3->Interpolate(Edep) - hNLnominal->Interpolate(Edep)) );
}


double JunoDetector::Resolution(double Evis)
{
    double a = m_a * (1 + JunoPullTerms::alpha_ea);
    double b = m_b * (1 + JunoPullTerms::alpha_eb);
    double c = m_c * (1 + JunoPullTerms::alpha_ec);

    return TMath::Sqrt(a*a/Evis + b*b + c*c/Evis/Evis) * Evis;
}


double JunoDetector::DetectorResponse(double Enu, double costheta, double Ep)
{
    double Te    = PositronEnergy(Enu, costheta) - 0.511;
    double Edep  = Te + 1.022;
    double Evis  = Nonlinearity(Edep);
    double sigma = Resolution(Evis);
    
    return 1/(sigma * TMath::Sqrt(2*TMath::Pi())) * TMath::Exp(-(Ep - Evis) * (Ep - Evis) /2/sigma/sigma);
}













