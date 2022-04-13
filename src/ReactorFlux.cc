#include "ReactorFlux.hh"

#include <iostream>
#include <TFile.h>
#include <TMath.h>

ReactorFlux::ReactorFlux()
{
    m_name        = "Reactor";              // name
    m_baseline    = 53000;                  // m
    m_power       = 2.9e9;                  // GW
    m_time        = 6 * 365 * 24 * 3600;    //s
    m_MO          = 1;


    fU235  = 0.58;
    fU238  = 0.07;
    fPu239 = 0.3;
    fPu241 = 0.05;

    eU235 = 202.36;
    eU238 = 205.99;
    ePu239 = 211.12;
    ePu241 = 214.26;

    sin2theta12 = 0.307;
    sin2theta13 = 2.18e-2;
    Deltam212 = 7.53e-5;
    Deltam322NO = 2.453e-3;
    Deltam322IO = -2.546e-3;

}


ReactorFlux::~ReactorFlux()
{
    delete hU235;
    delete hU238;
    delete hPu239;
    delete hPu241;
}


void ReactorFlux::LoadCommonInputs()
{
    TFile* ff = new TFile("JUNOInputs2021_05_28.root", "read");
    hU235 = (TH1D*)ff->Get("HuberMuellerFlux_U235");
    hU238 = (TH1D*)ff->Get("HuberMuellerFlux_U238");
    hPu239 = (TH1D*)ff->Get("HuberMuellerFlux_Pu239");
    hPu241 = (TH1D*)ff->Get("HuberMuellerFlux_Pu241");
}

double ReactorFlux::InitialReactorFlux(double Enu)
{
    double JtoMeV = 6.242e12; 
    double phi = m_power * m_time * JtoMeV / (fU235*eU235 + fU238*eU238 + fPu239*ePu239 + fPu241*ePu241) * (fU235*hU235->Interpolate(Enu) + fU238*hU238->Interpolate(Enu) + fPu239*hPu239->Interpolate(Enu) + fPu241*hPu241->Interpolate(Enu));
    return phi;
}

double ReactorFlux::SurvivalProbability(double Enu)
{
    double fast, slow;
    if (m_MO == 1) { // Normal Ordering
        fast = 4*sin2theta13*(1-sin2theta13)*((1-sin2theta12)*TMath::Sin((Deltam212 + Deltam322NO) * 1.27 * m_baseline /Enu)*TMath::Sin((Deltam212 + Deltam322NO) * 1.27 * m_baseline /Enu) + sin2theta12 * TMath::Sin(1.27*Deltam322NO*m_baseline/Enu)*TMath::Sin(1.27*Deltam322NO*m_baseline/Enu) );
        slow = (1 - sin2theta13)*(1-sin2theta13) * 4*sin2theta12*(1-sin2theta12) * TMath::Sin(1.27*Deltam212*m_baseline/Enu)*TMath::Sin(1.27*Deltam212*m_baseline/Enu);
    }
    else if (m_MO == 2) {   // Inverted Ordering
        fast = 4*sin2theta13*(1-sin2theta13)*((1-sin2theta12)*TMath::Sin((Deltam212 + Deltam322IO) * 1.27 * m_baseline /Enu)*TMath::Sin((Deltam212 + Deltam322IO) * 1.27 * m_baseline /Enu) + sin2theta12 * TMath::Sin(1.27*Deltam322IO*m_baseline/Enu)*TMath::Sin(1.27*Deltam322IO*m_baseline/Enu) );
        slow = (1 - sin2theta13)*(1-sin2theta13) * 4*sin2theta12*(1-sin2theta12) * TMath::Sin(1.27*Deltam212*m_baseline/Enu)*TMath::Sin(1.27*Deltam212*m_baseline/Enu);
    }
    else{
        cout << "Wrong Mass Ordering Input !!!!!!" << m_MO << endl;
    }

    return 1 - fast - slow;
}


double ReactorFlux::ArrivedReactorFlux(double Enu)
{
    double mtocm = 100;
    return SurvivalProbability(Enu) / (4*TMath::Pi()*m_baseline*m_baseline*mtocm*mtocm) * InitialReactorFlux(Enu);
}






