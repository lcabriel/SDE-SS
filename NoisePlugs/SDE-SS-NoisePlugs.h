#ifndef SDE_SS_NOISEPLUGS
#define SDE_SS_NOISEPLUGS

#include "SDE-SS.h"

using namespace std;

/*
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& INTRODUCTION &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
This sublibrary of SDE-SS contains the so-called NoisePlugs. The noise Plugs are used for systems that requires specific kind of noise which can
only be easily computed adding a variable to the system. This procedure will be still necessary using Plugs and they are not mandatory (most of them 
can be introduced by you directly) but the implementation is already made thus this should speed up your coding pipeline. 
How to use them is described below.
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

/*
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& HOW TO USE THE NOISEPLUGS &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Using the NoisePlugs is very simple. They are actually a little bit more than a normal function but with the possibility of storing the values of
the constants they use and, therefore, eventually change them during your work. To use them you must declare an extra variable for your system and
define the update of this variable in f_function_impl and g_function_impl with the deterministic_part and stochastic_part of your NoisePlug. It is
harder to explain than it seems therefore we have inserted a example code in this subdirectory that shows you how easy actually it is to use the 
NoisePlugs.
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
*/

//#################################################################################################################################################
//#################################### GENERIC NOISE PLUG CLASS ###################################################################################
//#################################################################################################################################################

//This class is used to represent a generic NoisePlug. This is used to represent compactly different kinds of complex noises. 
//It can be useful if your fields in the same codes can host different Plugs from time to time.
class NoisePlug{

public:

    //This function is used to compute the deterministic part of the SDE ( f(x) ) describing the noise accoding to Ito's formula.
    //It has to be override by the different NoisePlugs.
    virtual float deterministic_part(const vector<float> &x,float t);

    //This function is used to compute the stochastic part of the SDE ( g(x) ) describing the noise accoding to Ito's formula.
    //It has to be override by the different NoisePlugs.
    virtual float stochastic_part(const vector<float> &x,float t);

};

//#################################################################################################################################################
//#################################### WIENER PROCESS NOISEPLUGS ##################################################################################
//#################################################################################################################################################

/*
This NoisePlug allows us, adding an extra variable for the system, to simulate a Wiener Process in parallel with our simulations. This NoisePlug
can be very useful in some noises where there is a direct unremovable dependence from W(t) such as the Sine Wiener noise.
*/

//NoisePlug used to represent a simple Wiener process W(t).
class WienerProcessNP: public NoisePlug{

public:

    //Override of the NoisePlug deterministic_part. For the Wiener Process is simply f(x)=0.
    float deterministic_part(const vector<float> &x,float t) override;

    //Override of the NoisePlug stochastic_part. For the Wiener Process is simply g(x)=1.
    float stochastic_part(const vector<float> &x,float t) override;

};


























#endif