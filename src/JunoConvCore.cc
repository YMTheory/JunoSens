#include "JunoConvCore.hh"

using namespace std;

int JunoConvCore::m_MO;
ReactorFlux* JunoConvCore::reactor[N_reactor];
JunoDetector* JunoConvCore::det;

JunoConvCore::JunoConvCore()
{;}


JunoConvCore::~JunoConvCore()
{;}

void JunoConvCore::Initialize()
{
    //reactor = new ReactorFlux();
    //reactor->LoadCommonInputs();
    //reactor->SetMO(1);

    double baseline[N_reactor] = {52740, 52820, 52410, 52490, 52110, 52190, 52770, 52640, 215000, 265000};
    double power[N_reactor] = {2.9e9, 2.9e9, 2.9e9, 2.9e9, 2.9e9, 2.9e9, 4.6e9, 4.6e9, 17.4e9, 17.4e9}; 

    for(int i=0; i<N_reactor; i++) {
        reactor[i] = new ReactorFlux();
        reactor[i]->LoadCommonInputs();
        reactor[i]->SetBaseline(baseline[i]);
        reactor[i]->SetPower(power[i]);
        reactor[i]->SetMO(1);
    }

    det = new JunoDetector();
    det->LoadCommonInputs();
}

void JunoConvCore::SetMO(int MO) 
{
    for(int i=0; i<N_reactor; i++) {
        reactor[i] -> SetMO(m_MO);
    }
}



double JunoConvCore::fVisibleSpectrum(double* x, double* p)
{
    double Enu      = x[0];
    double costheta = x[1];
    double Ep       = x[2];

    //cout << ">>>>>>>>>>>>>> Current Integration <<<<<<<<<<<<<<<<<" << endl;
    //cout << Enu << " " << costheta << " " << Ep << endl;
    //cout << "ArrivedReactorFlux : " << reactor->ArrivedReactorFlux(Enu) << endl;
    //cout << "Nproton : " << det->GetNproton() << endl;
    //cout << "Efficiency : " << det->GetEfficiency() << endl;
    //cout << "IBDdiffXsec : " << det->IBDdiffXsec(Enu, costheta) << endl;
    //cout << "DetectorResponse : " << det->DetectorResponse(Enu, costheta, Ep) << endl;
    //cout << "----------------------------------------------------" << endl;
    
    double arrviedNuFlux = 0;
    for(int i=0; i<N_reactor; i++) {
        arrviedNuFlux += reactor[i]->ArrivedReactorFlux(Enu);
    }

    return arrviedNuFlux * det->IBDdiffXsec(Enu, costheta) * det->GetNproton() * det->GetEfficiency() * det->DetectorResponse(Enu, costheta, Ep);

}









