#ifndef JunoSpectrum_h
#define JunoSpectrum_h

#include "JunoIBDSignal.hh"
#include "JunoBackground.hh"

#include <TH1D.h>

class JunoSpectrum
{
    public:
        JunoSpectrum(int dataMO, int predMO);
        ~JunoSpectrum();

    
    public:
        TH1D* PredictedSpectrum();
        TH1D* MeasuredSpectrum();

        double GetChi2();

    public:
        TH1F* Loadb2bUncertaintyTAO();

    private:
        JunoIBDSignal* junoIBD;

        int m_data_MO;
        int m_pred_MO;

    public:
        int GetDataMO()         { return m_data_MO;}
        void SetDataMO(int val) { m_data_MO = val; }
        int GetPredMO()         { return m_pred_MO;}
        void SetPredMO(int val) { m_pred_MO = val; }


    private:
        TH1D* hPred;
        TH1D* hMea;

        TH1F* hb2bTAO;
        TH1D* hb2bModel;
        TH1D* hb2bDYB;
        TH1D* hb2bYB;


};

#endif
