#include "SDE-SS-NoisePlugs.h"
#include "SDE-SS.h"

using namespace std;

/*
&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& NOISEPLUG EXAMPLE &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

In this small example we want to show you how you can introduce the NoisePlugs in your fields. It is fairly simple.abort
We would like in this example to simulate this trivial system:

dx/dt = x*(1-x) + x*(1-x)*gamma(t)

where gamma(t) is a Sine Wiener Bounded noise described as

gamma(t) = B*sin(sqrt(2/tau)*W(t))

The system above can be considered as an SDE where g(x) is equal to 0. Using the gamma(t) definition instead of gamma(t)
we have still the problem to simulate W(t). Therefore, this is the perfect case to add an extra variable to our system
and use the WienerProcessNP NoisePlug!
*/

//Thus as first thing to do we need to create our custom FieldClass as in "test_code.cpp":

class myField: public FieldClass{

    //THE PARAMETERS OF MY FIELD:
    float B;
    float tau;
    WienerProcessNP W; //We need simply to define it not having constructor

public:

    //CONSTRUCTOR: Here we need to set the parameters for the system and, obviously, the noise.
    myField(float b,float t,NoiseClass* N): B(b), tau(t)
    {
        setNoise(N);
    }

    //Then we have to override f(x) and g(x)

    //In the f(x) we have in [0] all the field of the main system while we need to add
    //an extra variable to manage the W(t). There we will use the NoisePlug!
    vector<float> f_function_impl(const vector<float>& x,float t) override{
        vector<float> y(x.size(),0.0);

        y[0] = x[0]*(1-x[0])+x[0]*(1-x[0])*B*sin(sqrt(2.0/tau)*x[1]); //x[1] == W(t)
        y[1] = W.deterministic_part(x,t);

        return y;
    }

    //Similarly, for the stochastic part we will have at [0] that g(x)=0 thus we can leave it
    //as it is. For [1] we will use the stochastic_part of the NoisePlug.
    vector<float> g_function_impl(const vector<float>& x,float t) override{
        vector<float> y(x.size(),0.0);

        y[1] = W.stochastic_part(x,t);

        return y;
    }

};

//We have no bound thus we can move to the main function....

int main(){

    WienerEuler noise; //We use a simple WienerEuler
    myField field(1.0,100.0,&noise); //B=1, tau=100

    SDE_SS_System system(2,&field); //2 because we need an extra variable for W(t)

    //We can, as example, produce and print a trajectory
    vector<vector<float>> traj = system.simulateTrajectory({0.5,0.0},500,0.01); //We need also to give an initial value to the extra variable.

    for(size_t i=0;i<traj.size();i++){
        cout << traj[i][0] << " " << traj[i][1] << " " << traj[i][2] << endl;// t x(t) W(t)
    }

    return 0;
}