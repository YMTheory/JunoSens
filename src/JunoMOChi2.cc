#include "JunoMOChi2.hh"

#include <TH1D.h>
#include <TFile.h>


double JunoMOChi2::m_chi2;
double JunoMOChi2::m_chi2Min;
int    JunoMOChi2::m_nParameter;

JunoSpectrum* JunoMOChi2::junoSpec;


JunoMOChi2::JunoMOChi2()
{
    junoSpec = new JunoSpectrum();
    junoSpec->MeasuredSpectrum();
}

JunoMOChi2::~JunoMOChi2()
{
    delete junoSpec;
}


double JunoMOChi2::GetChi2(double maxChi2)
{
    double chi2 = junoSpec->GetChi2();
    if (chi2 > maxChi2)
        return maxChi2;
    else
        return chi2;
}



void JunoMOChi2::ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag)
{
    SetParameters(par);
    fval = GetChi2();
}


void JunoMOChi2::SetParameters(double* par){}




void JunoMOChi2::Plot()
{
    //TH1D* hIBD = junoIBD->PredictedVisibleEnergySpectrum();

    //TH1D* hAccidentalBkg    = JunoBackground::hAccidentalBkg;
    //TH1D* hLi9He8Bkg        = JunoBackground::hLi9He8Bkg;
    //TH1D* hGeoNeutrino      = JunoBackground::hGeoNeutrino;
    //TH1D* hFastNeutron      = JunoBackground::hFastNeutron;
    //TH1D* hAlphaNBkg        = JunoBackground::hAlphaNBkg;
    //TH1D* hGlobalReactorBkg = JunoBackground::hGlobalReactorBkg;
    //TH1D* hAtmNuBkg         = JunoBackground::hAtmNuBkg;

    //TFile* ff = new TFile("PredTotSpec.root", "recreate");
    //hIBD->Write();
    //hAccidentalBkg->Write();
    //hLi9He8Bkg->Write();
    //hFastNeutron->Write();
    //hGeoNeutrino->Write();
    //hAlphaNBkg->Write();
    //hGlobalReactorBkg->Write();
    //hAtmNuBkg->Write();

    //ff->Close();
}





