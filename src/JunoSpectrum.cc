#include "JunoSpectrum.hh"
#include "JunoPullTerms.hh"

#include <TFile.h>

JunoSpectrum::JunoSpectrum(int dataMO, int predMO)
{
    m_data_MO = dataMO;
    m_pred_MO = predMO;

    junoIBD = new JunoIBDSignal(predMO);
    //JunoBackground::LoadCommonInputs();
    //JunoBackground::CalculateBackground();

    //JunoPullTerms::LoadCommonInputs();
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
    cout << endl;
    cout << "===============================================" << endl;
    cout << "Data Mass Ordering is ";
    if (m_data_MO == 1) 
        cout << "Normal Ordering !" << endl;
    if (m_data_MO == 2)
        cout << "Inverted Ordering !" << endl;
    cout << "===============================================" << endl;
    cout << endl;


    //TFile* ff = new TFile("/junofs/users/miaoyu/JunoSens/ToyData/OnlySignalNoStat.root", "read");
    if (m_data_MO == 1) {
        TFile* ff = new TFile("/junofs/users/miaoyu/JunoSens/ToyData/Asimov_NO.root", "read");
        hMea = (TH1D*)ff->Get("hPredEvisSpec");
    }
    if (m_data_MO == 2){
        TFile* ff = new TFile("/junofs/users/miaoyu/JunoSens/ToyData/Asimov_IO.root", "read");
        hMea = (TH1D*)ff->Get("hPredEvisSpec");
    }
    
    return hMea;
}


TH1F* JunoSpectrum::Loadb2bUncertaintyTAO()
{
    cout << endl;
    cout << "===============================================" << endl;
    TFile* ff = new TFile("/junofs/users/miaoyu/JunoSens/JUNOInputs2022_01_06.root", "read");
    hb2bTAO = (TH1F*)ff->Get("TAOUncertainty");

    return hb2bTAO;

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
        wa_r += junoIBD->GetWeight(i, 0) * JunoPullTerms::alpha_r0;
        wa_r += junoIBD->GetWeight(i, 1) * JunoPullTerms::alpha_r1;
        wa_r += junoIBD->GetWeight(i, 2) * JunoPullTerms::alpha_r2;
        wa_r += junoIBD->GetWeight(i, 3) * JunoPullTerms::alpha_r3;
        wa_r += junoIBD->GetWeight(i, 4) * JunoPullTerms::alpha_r4;
        wa_r += junoIBD->GetWeight(i, 5) * JunoPullTerms::alpha_r5;
        wa_r += junoIBD->GetWeight(i, 6) * JunoPullTerms::alpha_r6;
        wa_r += junoIBD->GetWeight(i, 7) * JunoPullTerms::alpha_r7;
        wa_r += junoIBD->GetWeight(i, 8) * JunoPullTerms::alpha_r8;
        wa_r += junoIBD->GetWeight(i, 9) * JunoPullTerms::alpha_r9;

        double sig1 = (Mi - Ti * (1 + JunoPullTerms::alpha_C + wa_r + JunoPullTerms::alpha_D) );
        //double sig2 = Ti + (Ti * JunoPullTerms::GetBin2BinError(Evis)) * (Ti * JunoPullTerms::GetBin2BinError(Evis));
        //double sig2 = Ti ;
        double sig2 = Ti + (Ti * hb2bTAO->Interpolate(Evis))*(Ti * hb2bTAO->Interpolate(Evis));

        chi2 += sig1*sig1/sig2;
        //cout << i << " " << Evis << " " 
        //     << Mi << " " << Ti   << " " ;

        //for(int iReactor=0; iReactor<10; iReactor++) {
        //    cout << junoIBD->GetWeight(i+1, iReactor+1) <<  " ";
        //}
        //     cout << sig1 << " " << sig2 << " " << sig1*sig1/sig2 << " " << chi2 << endl;

    }

    // Pull Terms ...
    // ************************************************************* //
    // reactor correlated uncertainty :
    chi2 += (JunoPullTerms::alpha_C / JunoPullTerms::sigma_C) * (JunoPullTerms::alpha_C / JunoPullTerms::sigma_C);
    
    // reactor uncorrelated uncertainty :
    chi2 += (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r0) * (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r0);
    chi2 += (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r1) * (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r1);
    chi2 += (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r2) * (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r2);
    chi2 += (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r3) * (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r3);
    chi2 += (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r4) * (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r4);
    chi2 += (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r5) * (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r5);
    chi2 += (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r6) * (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r6);
    chi2 += (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r7) * (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r7);
    chi2 += (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r8) * (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r8);
    chi2 += (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r9) * (JunoPullTerms::alpha_r0 / JunoPullTerms::sigma_r9);

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

    // matter density:
    chi2 += (JunoPullTerms::alpha_rho/JunoPullTerms::sigma_rho) * (JunoPullTerms::alpha_rho/JunoPullTerms::sigma_rho);

    // Backgrounds...


    return chi2;
}






