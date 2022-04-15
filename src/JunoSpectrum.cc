#include "JunoSpectrum.hh"
#include "JunoPullTerms.hh"

#include <TFile.h>

JunoSpectrum::JunoSpectrum()
{
    junoIBD = new JunoIBDSignal(2);
    JunoBackground::LoadCommonInputs();
    JunoBackground::CalculateBackground();

    JunoPullTerms::LoadCommonInputs();
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
    //TFile* ff = new TFile("./ToyData/OnlySignalNoStat.root", "read");
    TFile* ff = new TFile("./ToyData/OnlySignalNoStat.root", "read");
    hMea = (TH1D*)ff->Get("OnlySignalStatNO");
    return hMea;
}



double JunoSpectrum::GetChi2()
{
    PredictedSpectrum(); // calculate predicted spectrum each time

    double chi2 = 0;
    for (int i=0; i<hMea->GetNbinsX(); i++) {

        double Evis = hMea->GetBinCenter(i+1);
        
        double Mi = hMea->GetBinContent(i+1);
        double Ti = hPred->GetBinContent(i+1);

        double wa_r = 0;
        for (int iReactor=0; iReactor<10; iReactor++) {
            wa_r += junoIBD->GetWeight(i+1, iReactor+1) * JunoPullTerms::alpha_r;
        }

        double sig1 = (Mi - Ti) * (1 + JunoPullTerms::alpha_C + wa_r * JunoPullTerms::alpha_D);
        double sig2 = Ti + (Ti * JunoPullTerms::GetBin2BinError(Evis)) * (Ti * JunoPullTerms::GetBin2BinError(Evis));
        
        chi2 += sig1*sig1/sig2;

        // Pull Terms ...

        // reactor correlated uncertainty :
        chi2 += (JunoPullTerms::alpha_C / JunoPullTerms::sigma_C) * (JunoPullTerms::alpha_C / JunoPullTerms::sigma_C);
        
        // reactor uncorrelated uncertainty :
        for (int iReactor=0; iReactor<10; iReactor++) {
            chi2 += (JunoPullTerms::alpha_r / JunoPullTerms::sigma_r) * (JunoPullTerms::alpha_r / JunoPullTerms::sigma_r);
        }

        // detector correlated uncertainty :
        chi2 += (JunoPullTerms::alpha_D / JunoPullTerms::sigma_D) * (JunoPullTerms::alpha_D / JunoPullTerms::sigma_D);

        // nonlinearity uncertainty :
        chi2 += (JunoPullTerms::alpha_l0/JunoPullTerms::sigma_l0) * (JunoPullTerms::alpha_l0/JunoPullTerms::sigma_l0);
        chi2 += (JunoPullTerms::alpha_l1/JunoPullTerms::sigma_l1) * (JunoPullTerms::alpha_l1/JunoPullTerms::sigma_l1);
        chi2 += (JunoPullTerms::alpha_l2/JunoPullTerms::sigma_l2) * (JunoPullTerms::alpha_l2/JunoPullTerms::sigma_l2);
        chi2 += (JunoPullTerms::alpha_l3/JunoPullTerms::sigma_l3) * (JunoPullTerms::alpha_l3/JunoPullTerms::sigma_l3);


        // resolution uncertainty:
        chi2 += (JunoPullTerms::alpha_ea/JunoPullTerms::sigma_ea) * (JunoPullTerms::alpha_ea/JunoPullTerms::sigma_ea);
        chi2 += (JunoPullTerms::alpha_eb/JunoPullTerms::sigma_eb) * (JunoPullTerms::alpha_eb/JunoPullTerms::sigma_eb);
        chi2 += (JunoPullTerms::alpha_ec/JunoPullTerms::sigma_ec) * (JunoPullTerms::alpha_ec/JunoPullTerms::sigma_ec);

        // SNF :
        chi2 += (JunoPullTerms::alpha_SNF/JunoPullTerms::sigma_SNF) * (JunoPullTerms::alpha_SNF/JunoPullTerms::sigma_SNF);

        // non-eq
        chi2 += (JunoPullTerms::alpha_NonEq/JunoPullTerms::sigma_NonEq) * (JunoPullTerms::alpha_NonEq/JunoPullTerms::sigma_NonEq);

    }

    return chi2;
}






