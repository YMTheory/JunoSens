#ifndef JunoDetector_h
#define JunoDetector_h

#include <TF2.h>
#include <TH1D.h>

using namespace std;

class JunoDetector
{
    public:
        JunoDetector();
        ~JunoDetector();

    public:
        void LoadCommonInputs();
        
        // IBD interaction
        double IBDdiffXsec(double Enu, double costheta);
        double PositronEnergy(double Enu, double costheta);
        double IBDtotXsec(double Enu);
        
        // detector response
        double Nonlinearity(double Edep);
        double Resolution(double Evis);
        double DetectorResponse(double Enu, double costheta, double Ep);

    public:
        void SetNproton(double nproton)             {m_Nproton = nproton;}
        double GetNproton()                         {return m_Nproton;}
        void SetEfficiency(double eff)              {m_effciency = eff;}
        double GetEfficiency()                      {return m_effciency;}

        void SetEres_a(double a)                    {m_a = a;}
        double GetEres_a()                          {return m_a;}
        void SetEres_b(double b)                    {m_b = b;}
        double GetEres_b()                          {return m_b;}
        void SetEres_c(double c)                    {m_c = c;}
        double GetEres_c()                          {return m_c;}


            
    private:
        double m_Nproton;
        double m_effciency;
        double m_a;
        double m_b;
        double m_c;

        TF2* fIBDdiffXsec;
        TF2* fEpositron;

        TH1D* hNLnominal;
        TH1D* hNL0;
        TH1D* hNL1;
        TH1D* hNL2;
        TH1D* hNL3;
        TH1D* hIBDtotXsec;
};

#endif
