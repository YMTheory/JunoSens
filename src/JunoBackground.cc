#include "JunoBackground.hh"

#include <TFile.h>

TH1D* JunoBackground::hAccidentalBkg;
TH1D* JunoBackground::hLi9He8Bkg;
TH1D* JunoBackground::hGeoNeutrino;
TH1D* JunoBackground::hFastNeutron;
TH1D* JunoBackground::hAlphaNBkg;
TH1D* JunoBackground::hGlobalReactorBkg;
TH1D* JunoBackground::hAtmNuBkg;

double JunoBackground::m_time;
double JunoBackground::m_GeoNu_perday;
double JunoBackground::m_Acc_perday;
double JunoBackground::m_FN_perday;
double JunoBackground::m_Li9He8_perday;
double JunoBackground::m_AlphaN_perday;
double JunoBackground::m_AtmNu_perday;
double JunoBackground::m_GlobalReactor_perday;

JunoBackground::JunoBackground()
{
}

JunoBackground::~JunoBackground()
{}


void JunoBackground::LoadCommonInputs()
{
    TFile* ff = new TFile("JUNOInputs2022_01_06.root", "read");
    hAccidentalBkg = (TH1D*)ff->Get("AccBkgHistogramAD");
    hLi9He8Bkg = (TH1D*)ff->Get("Li9BkgHistogramAD");
    hGeoNeutrino = (TH1D*)ff->Get("GeoNuHistogramAD");
    hAlphaNBkg = (TH1D*)ff->Get("AlphaNBkgHistogramAD");
    hFastNeutron = (TH1D*)ff->Get("FnBkgHistogramAD");
    hGlobalReactorBkg = (TH1D*)ff->Get("OtherReactorSpectrum_L300km");
    hAtmNuBkg = (TH1D*)ff->Get("AtmosphericNeutrinoModelNuWroN4");

}


void JunoBackground::CalculateBackground()
{
    m_time = 6 * 365 ;    // day
    m_Acc_perday = 0.8;
    m_FN_perday = 0.1;
    m_GeoNu_perday = 1.2;
    m_Li9He8_perday = 0.8;
    m_AlphaN_perday = 0.05;
    m_AtmNu_perday = 0.16;
    m_GlobalReactor_perday = 1;

    hAccidentalBkg->Scale(m_Acc_perday*m_time/hAccidentalBkg->Integral());
    hLi9He8Bkg->Scale(m_Li9He8_perday*m_time/hLi9He8Bkg->Integral());
    hGeoNeutrino->Scale(m_GeoNu_perday*m_time/hGeoNeutrino->Integral()); 
    hAlphaNBkg->Scale(m_AlphaN_perday*m_time/hAlphaNBkg->Integral());
    hFastNeutron->Scale(m_FN_perday*m_time/hFastNeutron->Integral());
    hGlobalReactorBkg->Scale(m_GlobalReactor_perday*m_time/hGlobalReactorBkg->Integral());
    hAtmNuBkg->Scale(m_AtmNu_perday*m_time/hAtmNuBkg->Integral());

}




