#include "JunoIBDSignal.hh"

#include <TFile.h>

JunoIBDSignal::JunoIBDSignal(int MO)
{
    m_MO = MO;

    cout << endl;
    cout << "===============================================" << endl;
    cout << "Prediction Mass Ordering is " ;
    if (MO == 1)
        cout << "Normal Ordering !" << endl;
    if (MO == 2)
        cout << "Inverted Ordering !" << endl;
    cout << "===============================================" << endl;
    cout << endl;

    //reactor = new ReactorFlux();
    //reactor->LoadCommonInputs();
    det = new JunoDetector();
    det->LoadCommonInputs();

    JunoConvCore::Initialize(m_MO);
    fVisibleEnergySpectrum = new TF3("fVisibleEnergySpectrum", JunoConvCore::fVisibleSpectrum, 1.8, 15, -1, 1, 0, 12, 1);


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

    hWeight = new TH2D("hWeight", "ith reactor ratio in jth bin", 340, 0, 340, 10, 0, 10);

    cout << "\n";
}


JunoIBDSignal::~JunoIBDSignal()
{
    //delete reactor;
    delete det;
    delete fVisibleEnergySpectrum;
}


//double JunoIBDSignal::BinnedNeutrinoEnergySpectrum(double Enu)
//{
//    double arrivedNu = reactor->ArrivedReactorFlux(Enu);
//    double nproton = det->GetNproton();
//    double eff = det->GetEfficiency();
//    double IBDXsec = det->IBDtotXsec(Enu);
//    
//    return arrivedNu * IBDXsec * nproton * eff;
//
//}


void JunoIBDSignal::CalculateReactorBinRatio()
{
    double reactorFlux[10];
    for(int ibin=0; ibin<340; ibin++) {
        double totFlux = 0;
        for( int i=0; i<10; i++) {
            reactorFlux[i] = BinnedVisibleEnergySpectrum(bin_edge[ibin], bin_edge[ibin+1], i);
            totFlux += reactorFlux[i];
        } 
        for( int i=0; i<10; i++) {
            hWeight->SetBinContent(ibin+1, i+1, reactorFlux[i]/totFlux);
        }
    }
}



double JunoIBDSignal::BinnedVisibleEnergySpectrum(double Epmin, double Epmax, int reactorNo)
{
    double sigma = det->Resolution(Epmax);
    double Enumin = Epmin - 5*sigma + 0.8;
    double Enumax = Epmax + 5*sigma + 0.8;
    fVisibleEnergySpectrum->SetParameter(0, reactorNo);
    return fVisibleEnergySpectrum->Integral(Enumin, Enumax, -1, 1, Epmin, Epmax, 1.e-3);

}




TH1D* JunoIBDSignal::PredictedVisibleEnergySpectrum()
{
    CalculateReactorBinRatio();

    for(int ibin=0; ibin<340; ibin++) {
        //cout << ibin << " " << bin_edge[ibin] << " " << bin_edge[ibin+1] << " " << BinnedVisibleEnergySpectrum(bin_edge[ibin], bin_edge[ibin+1]) << endl;
        hPredEvisSpec->SetBinContent(ibin+1, BinnedVisibleEnergySpectrum(bin_edge[ibin], bin_edge[ibin+1], 10) / (bin_edge[ibin+1] - bin_edge[ibin]) * 0.02);
    }
    return hPredEvisSpec;
}



void JunoIBDSignal::Plot()
{
    TFile* ff = new TFile("PredEvisSpec.root", "recreate");
    hPredEvisSpec->Write();
    ff->Close();
}







