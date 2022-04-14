#include "JunoSpectrum.hh"

#include <TFile.h>

JunoSpectrum::JunoSpectrum()
{
    junoIBD = new JunoIBDSignal(1);
    JunoBackground::LoadCommonInputs();
    JunoBackground::CalculateBackground();
}

JunoSpectrum::~JunoSpectrum()
{
    delete junoIBD;
}


TH1D* JunoSpectrum::PredictedSpectrum()
{
    hPred = junoIBD->PredictedVisibleEnergySpectrum();
    return hPred;
}


TH1D* JunoSpectrum::MeasuredSpectrum()
{
    TFile* ff = new TFile("./ToyData/OnlySignalNoStat.root", "read");
    hMea = (TH1D*)ff->Get("OnlySignalStatIO");
    return hMea;
}



double JunoSpectrum::GetChi2()
{
    double chi2 = 0;
    for (int i=0; i<hMea->GetNbinsX(); i++) {
        if (hPred->GetBinContent(i+1))
            chi2 += (hMea->GetBinContent(i+1) - hPred->GetBinContent(i+1)) * (hMea->GetBinContent(i+1) - hPred->GetBinContent(i+1)) / hPred->GetBinContent(i+1);
        else
            continue;
    }

    return chi2;
}






