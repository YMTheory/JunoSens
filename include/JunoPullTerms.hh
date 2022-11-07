#ifndef JunoPullTerms_h
#define JunoPullTerms_h

#include <TH1D.h>

class JunoPullTerms
{
    public:
        JunoPullTerms();
        ~JunoPullTerms();

    public:
        static double alpha_C;
        static double alpha_D;
        static double alpha_r0;
        static double alpha_r1;
        static double alpha_r2;
        static double alpha_r3;
        static double alpha_r4;
        static double alpha_r5;
        static double alpha_r6;
        static double alpha_r7;
        static double alpha_r8;
        static double alpha_r9;
        static double alpha_ME;
        static double alpha_ea;
        static double alpha_eb;
        static double alpha_ec;
        static double alpha_l0;
        static double alpha_l1;
        static double alpha_l2;
        static double alpha_l3;
        static double alpha_SNF;
        static double alpha_NonEq;
        static double alpha_rho;
        static double alpha_Acc;
        static double alpha_Li9He8;
        static double alpha_Atm;
        static double alpha_other;
        static double alpha_Geo;
        static double alpha_AlphaN;
        static double alpha_FN;

    public:
        static double sigma_C;
        static double sigma_D;
        static double sigma_r0;
        static double sigma_r1;
        static double sigma_r2;
        static double sigma_r3;
        static double sigma_r4;
        static double sigma_r5;
        static double sigma_r6;
        static double sigma_r7;
        static double sigma_r8;
        static double sigma_r9;
        static double sigma_ME;
        static double sigma_ea;
        static double sigma_eb;
        static double sigma_ec;
        static double sigma_l0;
        static double sigma_l1;
        static double sigma_l2;
        static double sigma_l3;
        static double sigma_SNF;
        static double sigma_NonEq;
        static double sigma_rho;
        static double sigma_Acc;
        static double sigma_Li9He8;
        static double sigma_Atm;
        static double sigma_other;
        static double sigma_Geo;
        static double sigma_AlphaN;
        static double sigma_FN;


    private:
        static TH1D* hBin2BinError;

    public:
        static void LoadCommonInputs();
        static double GetBin2BinError(double Evis)   {return hBin2BinError->Interpolate(Evis);}

};

#endif
