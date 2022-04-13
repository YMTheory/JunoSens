#include "JunoBackground.hh"

#include <TFile.h>

TH1D* JunoBackground::hAccidentalBkg;
TH1D* JunoBackground::hLi9He8Bkg;
TH1D* JunoBackground::hGeoNeutrino;
TH1D* JunoBackground::hFastNeutron;
TH1D* JunoBackground::hAlphaNBkg;


JunoBackground::JunoBackground()
{}

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

}







