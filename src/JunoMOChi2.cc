#include "JunoMOChi2.hh"

#include <TH1D.h>
#include <TFile.h>

JunoIBDSignal* JunoMOChi2::junoIBD;


double JunoMOChi2::m_chi2;
double JunoMOChi2::m_chi2Min;
int    JunoMOChi2::m_nParameter;


JunoMOChi2::JunoMOChi2()
{
    junoIBD = new JunoIBDSignal(1);

}

JunoMOChi2::~JunoMOChi2()
{
    delete junoIBD;
}


void JunoMOChi2::LoadData()
{
    JunoBackground::LoadCommonInputs();
    JunoBackground::CalculateBackground();
}



void JunoMOChi2::Plot()
{
    TH1D* hIBD = junoIBD->PredictedVisibleEnergySpectrum();

    //TH1D* hAccidentalBkg = JunoBackground::GetAccidentalBkg();
    //TH1D* hLi9He8Bkg     = JunoBackground::GetLi9H38Bkg();
    //TH1D* hGeoNeutrino   = JunoBackground::GetGeoNeutrino();
    //TH1D* hFastNeutron   = JunoBackground::GetFastNeutron();
    //TH1D* hAlphaNBkg     = JunoBackground::GetAlphaNBkg();
    TH1D* hAccidentalBkg    = JunoBackground::hAccidentalBkg;
    TH1D* hLi9He8Bkg        = JunoBackground::hLi9He8Bkg;
    TH1D* hGeoNeutrino      = JunoBackground::hGeoNeutrino;
    TH1D* hFastNeutron      = JunoBackground::hFastNeutron;
    TH1D* hAlphaNBkg        = JunoBackground::hAlphaNBkg;
    TH1D* hGlobalReactorBkg = JunoBackground::hGlobalReactorBkg;
    TH1D* hAtmNuBkg         = JunoBackground::hAtmNuBkg;

    TFile* ff = new TFile("PredTotSpec.root", "recreate");
    hIBD->Write();
    hAccidentalBkg->Write();
    hLi9He8Bkg->Write();
    hFastNeutron->Write();
    hGeoNeutrino->Write();
    hAlphaNBkg->Write();
    hGlobalReactorBkg->Write();
    hAtmNuBkg->Write();

    ff->Close();
}





