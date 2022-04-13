#include "JunoIBDSignal.hh"

#include <TFile.h>

JunoIBDSignal::JunoIBDSignal(int MO)
{
    m_MO = MO;

    reactor = new ReactorFlux();
    reactor->LoadCommonInputs();
    det = new JunoDetector();
    det->LoadCommonInputs();

    JunoConvCore::Initialize();
    //JunoConvCore::SetMO(m_MO);
    fVisibleEnergySpectrum = new TF3("fVisibleEnergySpectrum", JunoConvCore::fVisibleSpectrum, 1.8, 15, -1, 1, 0, 12);


    // binning strategy :
    int N = 0;
    for (int i=0; i<5; i++) {
        double bin_width = (zone_edge[i+1] - zone_edge[i]) / nbins[i];
        for (int j=0; j<nbins[i]; j++){
            bin_edge[N] = zone_edge[i] + j*bin_width;
            N++;
        }
    }
    bin_edge[340] = zone_edge[5];

    hPredEvisSpec = new TH1D("hPredEvisSpec", "Predicted visible energy spectrum", 340, bin_edge);
}


JunoIBDSignal::~JunoIBDSignal()
{
    delete reactor;
    delete det;
}


double JunoIBDSignal::BinnedNeutrinoEnergySpectrum(double Enu)
{
    double arrivedNu = reactor->ArrivedReactorFlux(Enu);
    double nproton = det->GetNproton();
    double eff = det->GetEfficiency();
    double IBDXsec = det->IBDtotXsec(Enu);
    
    return arrivedNu * IBDXsec * nproton * eff;

}


double JunoIBDSignal::BinnedVisibleEnergySpectrum(double Epmin, double Epmax)
{
    double sigma = det->Resolution(Epmax);
    double Enumin = Epmin - 5*sigma + 0.8;
    double Enumax = Epmax + 5*sigma + 0.8;
    return fVisibleEnergySpectrum->Integral(Enumin, Enumax, -1, 1, Epmin, Epmax, 1.e-6);

}


TH1D* JunoIBDSignal::PredictedVisibleEnergySpectrum()
{
    for(int ibin=0; ibin<340; ibin++) {
        cout << ibin << " " << bin_edge[ibin] << " " << bin_edge[ibin+1] << " " << BinnedVisibleEnergySpectrum(bin_edge[ibin], bin_edge[ibin+1]) << endl;
        hPredEvisSpec->SetBinContent(ibin+1, BinnedVisibleEnergySpectrum(bin_edge[ibin], bin_edge[ibin+1]) / (bin_edge[ibin+1] - bin_edge[ibin]) * 0.02);
    }
    return hPredEvisSpec;
}



void JunoIBDSignal::Plot()
{
    TFile* ff = new TFile("PredEvisSpec.root", "recreate");
    hPredEvisSpec->Write();
    ff->Close();
}







