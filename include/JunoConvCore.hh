#ifndef JunoConvCore_h
#define JunoConvCore_h

#include "ReactorFlux.hh"
#include "JunoDetector.hh"

class JunoConvCore
{
    public:
        JunoConvCore();
        ~JunoConvCore();

    public:
        static void Initialize(int MO);
        static double fVisibleSpectrum(double* x, double* p);

    private:
        static int m_MO;

        static const int N_reactor = 10;
        static ReactorFlux* reactor[N_reactor];
        static JunoDetector* det;


};

#endif
