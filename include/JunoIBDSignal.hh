#ifndef JunoIBDSignal_h
#define JunoIBDSignal_h

#include "ReactorFlux.hh"
#include "JunoDetector.hh"
#include "JunoConvCore.hh"

#include <TF3.h>
#include <TH1D.h>

class JunoIBDSignal
{
    public:
        JunoIBDSignal();
        ~JunoIBDSignal();

    public:
        double BinnedNeutrinoEnergySpectrum(double Enu);
        double BinnedVisibleEnergySpectrum(double Epmin, double Epmax);
        void PredictedVisibleEnergySpectrum();

        void Plot();

    private:
        ReactorFlux* reactor;
        JunoDetector* det;

    private:
        TF3* fVisibleEnergySpectrum;

    private:

        double zone_edge[6] = {0.8, 0.94, 7.44, 7.8, 8.2, 12};
        int nbins[5] = {1, 325, 9, 4, 1};
        double bin_edge[341];

    private:
        TH1D* hPredEvisSpec;


};



#endif
