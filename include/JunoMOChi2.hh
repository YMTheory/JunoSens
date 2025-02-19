#ifndef JunoMOChi2_h
#define JunoMOChi2_h

#include "JunoSpectrum.hh"

#include <TMinuit.h>

class JunoMOChi2
{
    public:
        JunoMOChi2();
        ~JunoMOChi2();

    public:
        double GetChiSquare( double maxChi2 = 100000 );
        static void SetParameters( double *par );
        static double GetChi2( double maxChi2 = 100000);

        static void Plot();


    private:

        static void ChisqFCN(Int_t &npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag);

        TMinuit* junoMinuit;

        static JunoSpectrum* junoSpec;

        static double m_chi2;
        static double m_chi2Min;
        static int m_nParameter;
        static double m_bestFit[20];
        static double m_bestFitError[20];

};

#endif
