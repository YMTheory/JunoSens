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
        static void CalculateBackground();
        
        static TH1D* GetAccidentalBkg()     { return hAccidentalBkg;}
        static TH1D* GetLi9H38Bkg()         { return hLi9He8Bkg;}
        static TH1D* GetGeoNeutrino()       { return hGeoNeutrino;}
        static TH1D* GetFastNeutron()       { return hFastNeutron;}
        static TH1D* GetAlphaNBkg()         { return hAlphaNBkg;}


        static void SetTime(double time)    { m_time = time; }
        static double GetTime()             { return m_time; }
    

    //private:
    public:
        static double m_time;
        
        static TH1D* hAccidentalBkg;
        static TH1D* hLi9He8Bkg;
        static TH1D* hGeoNeutrino;
        static TH1D* hFastNeutron;
        static TH1D* hAlphaNBkg;
        static TH1D* hGlobalReactorBkg;
        static TH1D* hAtmNuBkg;

        static double m_GeoNu_perday;
        static double m_Acc_perday;
        static double m_FN_perday;
        static double m_Li9He8_perday;
        static double m_AlphaN_perday;
        static double m_GlobalReactor_perday;
        static double m_AtmNu_perday;

};


#endif
