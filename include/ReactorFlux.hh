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
        void SetName(TString name)         {m_name = name;}
        TString GetName()                  {return m_name;}
        void SetBaseline(double baseline)  {m_baseline = baseline;}
        double GetBasline()                {return m_baseline;} 
        void SetPower(double power)        {m_power = power;}
        double GetPower()                  {return m_power;}
        void SetTime(double time)          {m_time = time;} 
        double GetTime()                   {return m_time;}
        void SetMO(int mo)                 {m_MO = mo;}
        double GetMO()                     {return m_MO;}

        void SetAlphaSNF(double snf)       {alpha_SNF = snf;}
        double GetAlphaSNF()               {return alpha_SNF;}
        void SetAlphaNonEq(double noneq)   {alpha_NonEq = noneq;}
        double GetAlphaNonEq()             {return alpha_NonEq;}
    
        void   LoadCommonInputs();
        double InitialReactorFlux(double Enu);
        double SurvivalProbability(double Enu);
        double ArrivedReactorFlux(double Enu);

    
    private:
        TString m_name;
        double m_baseline;
        double m_power;
        double m_time;
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

        TH1D* hU235;
        TH1D* hU238;
        TH1D* hPu239;
        TH1D* hPu241;

        double fSNF;
        double alpha_SNF;
        double fNonEq;
        double alpha_NonEq;

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
