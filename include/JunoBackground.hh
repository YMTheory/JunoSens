#ifndef JunoBackground_h
#define JunoBackground_h

#include <TH1D.h>

class JunoBackground
{
    public:
        JunoBackground();
        ~JunoBackground();


    public:
        static void LoadCommonInputs();
        static TH1D* GetAccidentalBkg()     { return hAccidentalBkg;}
        static TH1D* GetLi9H38Bkg()         { return hLi9He8Bkg;}
        static TH1D* GetGeoNeutrino()       { return hGeoNeutrino;}
        static TH1D* GetFastNeutron()       { return hFastNeutron;}
        static TH1D* GetAlphaNBkg()         { return hAlphaNBkg;}

    
    private:
        static TH1D* hAccidentalBkg;
        static TH1D* hLi9He8Bkg;
        static TH1D* hGeoNeutrino;
        static TH1D* hFastNeutron;
        static TH1D* hAlphaNBkg;
};


#endif
