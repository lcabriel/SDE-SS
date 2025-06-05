#include "SDE-SS-NoisePlugs.h"

using namespace std;

//#################################################################################################################################################
//#################################### GENERIC NOISE PLUG CLASS ###################################################################################
//#################################################################################################################################################

//This function is used to compute the deterministic part of the SDE ( f(x) ) describing the noise accoding to Ito's formula.
//It has to be override by the different NoisePlugs.
float NoisePlug::deterministic_part(const vector<float> &x,float t){
    return 0.0;
}

//This function is used to compute the stochastic part of the SDE ( g(x) ) describing the noise accoding to Ito's formula.
//It has to be override by the different NoisePlugs.
float NoisePlug::stochastic_part(const vector<float> &x,float t){
    return 0.0;
}

//#################################################################################################################################################
//#################################### WIENER PROCESS NOISEPLUGS ##################################################################################
//#################################################################################################################################################

//Override of the NoisePlug deterministic_part. For the Wiener Process is simply f(x)=0.
float WienerProcessNP::deterministic_part(const vector<float> &x,float t){
    return 0.0;
}

//Override of the NoisePlug stochastic_part. For the Wiener Process is simply g(x)=1.
float WienerProcessNP::stochastic_part(const vector<float> &x,float t){
    return 1.0;
}