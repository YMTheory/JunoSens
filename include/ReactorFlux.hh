#ifndef ReactorFlux_h
#define ReactorFlux_h

#include <TH1D.h>
#include <TString.h>
#include <vector>

using namespace std;

class ReactorFlux
{
    public :
        ReactorFlux();
        ~ReactorFlux();

    public:
        void SetName(TString name)          {m_name = name;}
        TString GetName()                   {return m_name;}
        void SetBaseline(double baseline)   {m_baseline = baseline;}
        double GetBaseline()                {return m_baseline;} 
        void SetPower(double power)         {m_power = power;}
        double GetPower()                   {return m_power;}
        void SetTime(double time)           {m_time = time;} 
        double GetTime()                    {return m_time;}
        void SetDutyCycle(double cycle)     {m_duty_cycle = cycle;}
        double GetDutyCycle()               {return m_duty_cycle;}
        void SetDensity(double rho)         {m_rho = rho;}
        double GetDensity()                 { return m_rho;}
        void SetMO(int mo)                  {m_MO = mo;}
        double GetMO()                      {return m_MO;}

        void   LoadCommonInputs();
        double InitialReactorFlux(double Enu);
        double SurvivalProbability(double Enu);
        double SurvivalProbability_matter(double Enu);
        double ArrivedReactorFlux(double Enu);

    
    private:
        TString m_name;
        double m_baseline;
        double m_power;
        double m_time;
        double m_duty_cycle;
        double m_MO;

    private:
        double fU235;
        double fU238;
        double fPu239;
        double fPu241;
        double eU235;
        double eU238;
        double ePu239;
        double ePu241;


    private:
        double m_rho;


        TH1D* hU235;
        TH1D* hU238;
        TH1D* hPu239;
        TH1D* hPu241;

        double fSNF;
        double fNonEq;

        TH1D* hSNF;
        TH1D* hNonEq;

        TH1D* hDYBratio;


    private:
        double sin2theta12;
        double sin2theta13;
        double Deltam212;
        double Deltam322NO;
        double Deltam322IO;
    

};

#endif
