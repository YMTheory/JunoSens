#include "JunoMOChi2.hh"
#include "JunoPullTerms.hh"

#include <TH1D.h>
#include <TFile.h>


double JunoMOChi2::m_chi2;
double JunoMOChi2::m_chi2Min;
int    JunoMOChi2::m_nParameter;
double JunoMOChi2::m_bestFit[20];
double JunoMOChi2::m_bestFitError[20];

JunoSpectrum* JunoMOChi2::junoSpec;


JunoMOChi2::JunoMOChi2()
{
    junoSpec = new JunoSpectrum(1, 2);
    junoSpec->MeasuredSpectrum();
    junoSpec->Loadb2bUncertaintyTAO();
}

JunoMOChi2::~JunoMOChi2()
{
    delete junoSpec;
}


double JunoMOChi2::GetChi2(double maxChi2)
{
    double chi2 = junoSpec->GetChi2();
    cout << ">>>>>>>>>>>>>>>> Current chi2 = " << chi2 << endl;
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


void JunoMOChi2::SetParameters(double* par){
    JunoPullTerms::alpha_C = par[0];
    JunoPullTerms::alpha_r0 = par[1];
    JunoPullTerms::alpha_r1 = par[2];
    JunoPullTerms::alpha_r2 = par[3];
    JunoPullTerms::alpha_r3 = par[4];
    JunoPullTerms::alpha_r4 = par[5];
    JunoPullTerms::alpha_r5 = par[6];
    JunoPullTerms::alpha_r6 = par[7];
    JunoPullTerms::alpha_r7 = par[8];
    JunoPullTerms::alpha_r8 = par[9];
    JunoPullTerms::alpha_r9 = par[10];
    JunoPullTerms::alpha_D = par[11];
    JunoPullTerms::alpha_l0 = par[12];
    JunoPullTerms::alpha_l1 = par[13];
    JunoPullTerms::alpha_l2 = par[14];
    JunoPullTerms::alpha_l3 = par[15];
    JunoPullTerms::alpha_ea = par[16];
    JunoPullTerms::alpha_eb = par[17];
    JunoPullTerms::alpha_ec = par[18];
    JunoPullTerms::alpha_SNF = par[19];
    JunoPullTerms::alpha_NonEq = par[20];
    JunoPullTerms::alpha_rho = par[21];
}


double JunoMOChi2::GetChiSquare(double maxChi2)
{
    junoMinuit = new TMinuit();
    junoMinuit->SetFCN(ChisqFCN);
    junoMinuit->SetPrintLevel(1);

    double arglist[10];
    int ierrflag = 0;

    int iPar = 0;
    junoMinuit->mnexcm("CLEAR", arglist, 0, ierrflag);

    junoMinuit->mnparm(iPar, "alpha_C",     0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_r0",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_r1",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_r2",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_r3",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_r4",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_r5",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_r6",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_r7",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_r8",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_r9",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_D",     0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_l0",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_l1",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_l2",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_l3",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_ea",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_eb",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_ec",    0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_SNF",   0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_NonEq", 0, 0.001, -1, 1, ierrflag); iPar++;
    junoMinuit->mnparm(iPar, "alpha_rho",   0, 0.001, -1, 1, ierrflag); iPar++;


    // Minimization strategy
    junoMinuit->SetErrorDef(1);
    arglist[0]=2;
    junoMinuit->mnexcm("SET STR",arglist,1,ierrflag);

    arglist[0] = 5000; //maxCalls
    arglist[1] = 0.01; // tolerance
    junoMinuit->mnexcm("MIGrad", arglist, 1, ierrflag);

    junoMinuit->fCstatu.Data();

    double min, edm, errdef;
    int nvpar, nparx, icstat;
    junoMinuit->mnstat(min, edm, errdef, nvpar, nparx, icstat);

    m_nParameter = junoMinuit->GetNumPars();
	for(int i=0; i<m_nParameter; i++)
	{
	    junoMinuit->GetParameter(i, m_bestFit[i], m_bestFitError[i]);
	}


    cout << " ====================== " << endl;
    cout << "    minChi2: " << min << " with nDataBin = " << 340 << " and nPar = " << m_nParameter << endl;
    cout << " ====================== " << endl;

    for(int i=0; i<m_nParameter; i++) {
        cout << m_bestFit[i]  << " " ;
    }
    cout << endl;
    for(int i=0; i<m_nParameter; i++) {
        cout << m_bestFitError[i]  << " " ;
    }
    cout << endl;
    
    delete junoMinuit;
    return min;




}






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





