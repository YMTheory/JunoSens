#ifndef JunoIBDSignal_h
#define JunoIBDSignal_h

#include "ReactorFlux.hh"
#include "JunoDetector.hh"
#include "JunoConvCore.hh"

#include <TF3.h>
#include <TH1D.h>
#include <TH2D.h>

class JunoIBDSignal
{
    public:
        JunoIBDSignal(int MO);
        ~JunoIBDSignal();

    public:
        void SetMO(int MO)   {m_MO = MO;}

        //double BinnedNeutrinoEnergySpectrum(double Enu);
        double BinnedVisibleEnergySpectrum(double Epmin, double Epmax, int reactorNo);
        void CalculateReactorBinRatio();
        TH1D* PredictedVisibleEnergySpectrum();

        void Plot();

    private:
        //ReactorFlux* reactor;
        JunoDetector* det;

    private:
        TF3* fVisibleEnergySpectrum;

    private:

        double zone_edge[6] = {0.8, 0.94, 7.44, 7.8, 8.2, 12};
        int nbins[5] = {1, 325, 9, 4, 1};
        double bin_edge[341];

        int m_MO;

    private:
        TH1D* hPredEvisSpec;
        TH2D* hWeight;


};



#endif
