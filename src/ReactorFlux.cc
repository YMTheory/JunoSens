#include "ReactorFlux.hh"
#include "JunoPullTerms.hh"

#include <iostream>
#include <TFile.h>
#include <TMath.h>

// Flux arrived the detector for a single reactor core

ReactorFlux::ReactorFlux()
{
    m_name        = "Reactor";              // name
    m_baseline    = 53000;                  // m
    m_power       = 2.9e9;                  // GW
    m_time        = 6 * 365.25 * 24 * 3600;    //s
    m_duty_cycle  = 1; //11./12;
    m_MO          = 1;

    m_rho = 2.45;  // +- 0.15 g/cm3

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
    delete hSNF;
    delete hNonEq;
    delete hDYBratio;
}


void ReactorFlux::LoadCommonInputs()
{
    TFile* ff = new TFile("JUNOInputs2022_01_06.root", "read");
    hU235 = (TH1D*)ff->Get("HuberMuellerFlux_U235");
    hU238 = (TH1D*)ff->Get("HuberMuellerFlux_U238");
    hPu239 = (TH1D*)ff->Get("HuberMuellerFlux_Pu239");
    hPu241 = (TH1D*)ff->Get("HuberMuellerFlux_Pu241");

    hSNF = (TH1D*)ff->Get("SNF_FluxRatio");
    hNonEq = (TH1D*)ff->Get("NonEq_FluxRatio");
    hDYBratio = (TH1D*)ff->Get("DYBFluxBump_ratio");

}

double ReactorFlux::InitialReactorFlux(double Enu)
{
    double JtoMeV = 6.242e12; 
    double phi = m_power * m_time * m_duty_cycle * JtoMeV / (fU235*eU235 + fU238*eU238 + fPu239*ePu239 + fPu241*ePu241) * (fU235*hU235->Interpolate(Enu) + fU238*hU238->Interpolate(Enu) + fPu239*hPu239->Interpolate(Enu) + fPu241*hPu241->Interpolate(Enu)); // unit : MeV-1
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

double ReactorFlux::SurvivalProbability_matter(double Enu) 
{
    double GF = 1.166e-5; // GeV-2
    double NA = 6.022e23;
    double Ye = 0.5;
    double a = -1.52e-4 * Ye * m_rho * (1+JunoPullTerms::alpha_rho) * 1e-3;   // unit : eV2, Enu should be in unit of MeV
    double fast, slow;

    if (m_MO == 1) {
        double dmee2 = (1 - sin2theta12)* (Deltam322NO + Deltam212) + sin2theta12 * Deltam322NO;
        double dmee2_bar = dmee2 * TMath::Sqrt(((1-2*sin2theta13) - a / dmee2)*((1-2*sin2theta13) - a / dmee2) + 4 * sin2theta13 * (1-sin2theta13));
        double costwotheta13_bar = (dmee2*(1-2*sin2theta13) - a) / (dmee2_bar);
        double a_bar = (a + dmee2 - dmee2_bar) / 2;
        double cos2dtheta = (dmee2_bar + dmee2 - a * (1-2*sin2theta13)) / (2 * dmee2_bar);
        double Deltam212_bar = Deltam212 * TMath::Sqrt(((1-2*sin2theta12)-a_bar/Deltam212) * ((1-2*sin2theta12)-a_bar/Deltam212) + cos2dtheta * 4 * sin2theta12 * (1-sin2theta12));
        double costwotheta12_bar = (Deltam212 * (1-2*sin2theta12) - a_bar) / Deltam212_bar;
        double Deltam312_bar = dmee2_bar + (1 - costwotheta12_bar)/2 * Deltam212_bar;
        double Deltam322_bar = Deltam312_bar - Deltam212_bar;

        double sin2theta13_bar = (1 - costwotheta13_bar) / 2;
        double sin2theta12_bar = (1 - costwotheta12_bar) / 2;

        //cout << "Deltam212 "    << Deltam212 << " " << Deltam212_bar << "\n"
        //     << "Deltam322 "    << Deltam322NO << " " << Deltam322_bar << "\n"
        //     << "sin2theta12 "  << sin2theta12 << " " << ((1 - costwotheta12_bar)/2. )<< "\n"
        //     << "sin2theta13 "  << sin2theta13 << " " << ((1 - costwotheta13_bar)/2. );
        fast = 4*sin2theta13_bar*(1-sin2theta13_bar)*((1-sin2theta12_bar)*TMath::Sin((Deltam212_bar + Deltam322_bar) * 1.27 * m_baseline /Enu)*TMath::Sin((Deltam212_bar + Deltam322_bar) * 1.27 * m_baseline /Enu) + sin2theta12_bar* TMath::Sin(1.27*Deltam322_bar*m_baseline/Enu)*TMath::Sin(1.27*Deltam322_bar*m_baseline/Enu) );
        slow = (1 - sin2theta13_bar)*(1-sin2theta13_bar) * 4*sin2theta12_bar*(1-sin2theta12_bar) * TMath::Sin(1.27*Deltam212_bar*m_baseline/Enu)*TMath::Sin(1.27*Deltam212_bar*m_baseline/Enu);
    }

    if (m_MO == 2) {
        double dmee2 = (1 - sin2theta12)* (Deltam322IO + Deltam212) + sin2theta12 * Deltam322IO;
        double dmee2_bar = dmee2 * TMath::Sqrt(((1-2*sin2theta13) - a / dmee2)*((1-2*sin2theta13) - a / dmee2) + 4 * sin2theta13 * (1-sin2theta13));
        double costwotheta13_bar = (dmee2*(1-2*sin2theta13) - a) / (dmee2_bar);
        double a_bar = (a + dmee2 - dmee2_bar) / 2;
        double cos2dtheta = (dmee2_bar + dmee2 - a * (1-2*sin2theta13)) / (2 * dmee2_bar);
        double Deltam212_bar = Deltam212 * TMath::Sqrt(((1-2*sin2theta12)-a_bar/Deltam212) * ((1-2*sin2theta12)-a_bar/Deltam212) + cos2dtheta * 4 * sin2theta12 * (1-sin2theta12));
        double costwotheta12_bar = (Deltam212 * (1-2*sin2theta12) - a_bar) / Deltam212_bar;
        double Deltam312_bar = dmee2_bar + (1 - costwotheta12_bar)/2 * Deltam212_bar;
        double Deltam322_bar = Deltam312_bar - Deltam212_bar;

        double sin2theta13_bar = (1 - costwotheta13_bar) / 2;
        double sin2theta12_bar = (1 - costwotheta12_bar) / 2;

        fast = 4*sin2theta13_bar*(1-sin2theta13_bar)*((1-sin2theta12_bar)*TMath::Sin((Deltam212_bar + Deltam322_bar) * 1.27 * m_baseline /Enu)*TMath::Sin((Deltam212_bar + Deltam322_bar) * 1.27 * m_baseline /Enu) + sin2theta12_bar* TMath::Sin(1.27*Deltam322_bar*m_baseline/Enu)*TMath::Sin(1.27*Deltam322_bar*m_baseline/Enu) );
        slow = (1 - sin2theta13_bar)*(1-sin2theta13_bar) * 4*sin2theta12_bar*(1-sin2theta12_bar) * TMath::Sin(1.27*Deltam212_bar*m_baseline/Enu)*TMath::Sin(1.27*Deltam212_bar*m_baseline/Enu);
    }

    return 1 - fast - slow;

}


double ReactorFlux::ArrivedReactorFlux(double Enu)
{
    double mtocm = 100;
    double phir = SurvivalProbability_matter(Enu) / (4*TMath::Pi()*m_baseline*m_baseline*mtocm*mtocm) * InitialReactorFlux(Enu);   // unit: cm-2 * MeV-1
    //double phir = SurvivalProbability(Enu) / (4*TMath::Pi()*m_baseline*m_baseline*mtocm*mtocm) * InitialReactorFlux(Enu);   // unit: cm-2 * MeV-1
    return phir * ( 1 + hSNF->Interpolate(Enu) * (1+JunoPullTerms::alpha_SNF) + hNonEq->Interpolate(Enu) * (1 + JunoPullTerms::alpha_NonEq) ) * hDYBratio->Interpolate(Enu) ;
}






