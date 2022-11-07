#include "junosens.hh"
#include "JunoSpectrum.hh"
#include "JunoIBDSignal.hh"
#include "ReactorFlux.hh"

int main()
{
    JunoIBDSignal* junoIBD = new JunoIBDSignal(2);
    junoIBD->PredictedVisibleEnergySpectrum();
    junoIBD->Plot();

    //JunoSpectrum* junoSpec = new JunoSpectrum(1);
    //junoSpec->MeasuredSpectrum();
    //junoSpec->PredictedSpectrum();
    //cout << junoSpec->GetChi2();

    //JunoMOChi2* junoMO = new JunoMOChi2();
    //junoMO->GetChiSquare();
    //junoMO->Plot();

}

