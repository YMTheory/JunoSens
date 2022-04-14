#ifndef JunoSpectrum_h
#define JunoSpectrum_h

#include "JunoIBDSignal.hh"
#include "JunoBackground.hh"

#include <TH1D.h>

class JunoSpectrum
{
    public:
        JunoSpectrum();
        ~JunoSpectrum();

    
    public:
        TH1D* PredictedSpectrum();
        TH1D* MeasuredSpectrum();

        double GetChi2();

    private:
        JunoIBDSignal* junoIBD;


    private:
        TH1D* hPred;
        TH1D* hMea;


};

#endif
