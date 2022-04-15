#include "JunoPullTerms.hh"

#include <TFile.h>

double JunoPullTerms::alpha_C = 0;
double JunoPullTerms::alpha_D = 0;
double JunoPullTerms::alpha_r = 0;
double JunoPullTerms::alpha_ME = 0;
double JunoPullTerms::alpha_ea = 0;
double JunoPullTerms::alpha_eb = 0;
double JunoPullTerms::alpha_ec = 0;
double JunoPullTerms::alpha_l0 = 0;
double JunoPullTerms::alpha_l1 = 0;
double JunoPullTerms::alpha_l2 = 0;
double JunoPullTerms::alpha_l3 = 0;
double JunoPullTerms::alpha_SNF = 0;
double JunoPullTerms::alpha_NonEq = 0;
double JunoPullTerms::alpha_Acc = 0;
double JunoPullTerms::alpha_Li9He8 = 0;
double JunoPullTerms::alpha_Atm = 0;
double JunoPullTerms::alpha_other = 0;
double JunoPullTerms::alpha_Geo = 0;
double JunoPullTerms::alpha_AlphaN = 0;
double JunoPullTerms::alpha_FN = 0;

double JunoPullTerms::sigma_C = 0.02;
double JunoPullTerms::sigma_D = 0.01;
double JunoPullTerms::sigma_r = 0.008;
double JunoPullTerms::sigma_ME = 0.06;
double JunoPullTerms::sigma_ea = 0.008;
double JunoPullTerms::sigma_eb = 0.012;
double JunoPullTerms::sigma_ec = 0.033;
double JunoPullTerms::sigma_l0 = 0.005;
double JunoPullTerms::sigma_l1 = 0.005;
double JunoPullTerms::sigma_l2 = 0.005;
double JunoPullTerms::sigma_l3 = 0.005;
double JunoPullTerms::sigma_SNF = 0.3;
double JunoPullTerms::sigma_NonEq = 0.3;
double JunoPullTerms::sigma_Acc = 0.01;
double JunoPullTerms::sigma_Li9He8 = 0.2;
double JunoPullTerms::sigma_Atm = 0.5;
double JunoPullTerms::sigma_other = 0.02;
double JunoPullTerms::sigma_Geo = 0.3;
double JunoPullTerms::sigma_AlphaN = 0.5;
double JunoPullTerms::sigma_FN = 1;

TH1D* JunoPullTerms::hBin2BinError;

JunoPullTerms::JunoPullTerms()
{}

JunoPullTerms::~JunoPullTerms()
{;}

void JunoPullTerms::LoadCommonInputs()
{
    TFile* ff = new TFile("JUNOInputs2022_01_06.root", "read");
    hBin2BinError = (TH1D*)ff->Get("YBUncertainty");
    
}
