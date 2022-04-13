#include "JunoDetector.hh"

#include <TMath.h>
#include <TFile.h>

JunoDetector::JunoDetector()
{
    m_Nproton = 1.43512e33;
    m_effciency = 0.822;

    m_a = 0.0261;
    m_b = 0.0082;
    m_c = 0.0123;
}


JunoDetector::~JunoDetector()
{}


void JunoDetector::LoadCommonInputs()
{
    TFile* ff = new TFile("JUNOInputs2021_05_28.root", "read");
    fIBDdiffXsec = (TF2*)ff->Get("dsigma_dcos_Enu_cos_DYB");
    fEpositron   = (TF2*)ff->Get("Epositron_Enu_cos_DYB");
    hNL          = (TH1D*)ff->Get("positronScintNL");
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
    return Edep * hNL->Interpolate(Edep);
}


double JunoDetector::Resolution(double Evis)
{
    return TMath::Sqrt(m_a*m_a/Evis + m_b*m_b + m_c*m_c/Evis/Evis) * Evis;
}


double JunoDetector::DetectorResponse(double Enu, double costheta, double Ep)
{
    double Te    = PositronEnergy(Enu, costheta) - 0.511;
    double Edep  = Te + 1.022;
    double Evis  = Nonlinearity(Edep);
    double sigma = Resolution(Evis);
    
    return 1/(sigma * TMath::Sqrt(2*TMath::Pi())) * TMath::Exp(-(Ep - Evis) * (Ep - Evis) /2/sigma/sigma);
}













