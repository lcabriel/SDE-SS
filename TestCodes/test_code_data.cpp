#include "SDE-SS.h" //The library
#include <fstream> //Used to write the file

using namespace std;

/*
In this example we will show you how to implement a field where the data linker is needed. To do this we are going to implement a system to 
produce a trajectory and another one then based on a data linker.
*/

//######################################################### FIELD TO PRODUCE THE DATA ##################################################################################
/*1. We need at first to create our production field. To do this we need to create a custom version of the Field class. Defining as variable all 
the parameters we are going to use.

2. In the construction we need also to pass the noise. The class is made to accept any NoiseClass object however NoiseClass is purely virtual.
You can create one of your own looking to the standard ones defined in "SDE-SS.cpp" however the Gaussian White noise is already implemented both
using the Euler-Maruyama method and the Milstein one. They are called "WienerEuler" and "WienerMilstein". 
To set the Noise properly we suggest to use the internal FieldClass function "setNoise".

3. A generic system of SDEs can be written in Ito's formula as

        dx = f(x,t)dt + g(x,t)dW

If the noise is colored can be thought most of the time in this way supposing to add a dimension to define the noise propeties.
For other kinds of noises which cannot be expressed in this way, you are free to create your own class or, please, send it to us as a suggestion.
Coming back, the FieldClass is created around the Ito's formula thus we need to characterised the f(x,t) and g(x,t), to do this you need to override the
"f_function_impl" and "g_function_impl" as in the example.
*/

class myField: public FieldClass{

    //THE PARAMETERS
    float x_e; 
    float tau1;
    float tau2;
    float D;

public:

    //THE CONSTRUCTOR
    myField(float a,float b,float c,float d,NoiseClass* N){
        //Setting the values of the parameters
        x_e = a;
        tau1 = b;
        tau2 = c;
        D = d;

        //For the noise we use specifically this specific function  
        setNoise(N);
    }

    //The f_function definition: we are going to say as each variable of the system
    //should be update. This design is similar to the system of equations' shape.
    void f_function_impl(const vector<float> &x,float t,vector<float> &y) const override{

        y[0] = x[0]*(1-x[0])*(x_e - x[1]);
        y[1] = (x[2]-x[1])/tau2;
        y[2] = (x[0]-x[2])/tau1;
    }

    //The same can be said for the g_function.
    void g_function_impl(const vector<float> &x,float t,vector<float> &y) const override{

        y[0] = -x[0]*(1-x[0])*D;
        y[1] = 0;
        y[2] = 0;
    }

};


//The system described by the equations above it is also a bounded one in [0,1]. To this we need to define a bound function to pass to
//our system. The bound function expression is standard: it returns a boolean and requires a vector<float> as argument.
bool boundaries(vector<float> x){
    return ((x[0]>=0)&&(x[0]<=1));
}

//########################################################### CUSTOM DATALINKER #########################################################################################

//We have to define our custom version of the DataLinker.
class myDataLinker: public DataLinker{

public:

    //In the constructor we pass the data of the trajectory.
    myDataLinker(const Traj& traj,const bool AS): DataLinker(traj,AS){}

    //To get the first variable.
    float getData(const float t,const vector<float>& x) override{
        return data.getVars()[findTimeIndex(t)][0];
    }

    //To get the second variable.
    float getData2(const float t,const vector<float>& x){
        return data.getVars()[TOC][1];
    }

};

//########################################################### FIELD WITH DATA LINKER ####################################################################################

class DataField: public FieldClass{

    float tau; //a parameter
    myDataLinker* linker; //the reference to the data linker; we need to define it in the main the instance

public:

    //CONSTRUCTOR
    DataField(float t,myDataLinker* DL,NoiseClass* N){
        tau = t;
        linker = DL;
        setNoise(N);
    }

    //The f_function definition: let's suppose that the field with data linker is only 2D because is a different system.
    void f_function_impl(const vector<float> &x,float t,vector<float> &y) const override{
        y[0] = (linker->getData(t,x)-x[0])/tau;
        y[1] = (linker->getData2(t,x)-x[1])/tau;
    }

    //The same can be said for the g_function.
    void g_function_impl(const vector<float> &x,float t,vector<float> &y) const override{
        y[0] = 0;
        y[1] = 0;
    }

};

//####################################################################### MAIN ##########################################################################################
/*
In the main, we are going to show use and example of usage producing the PDF along the first variable for our system.
*/

int main(){
    WienerEuler noise1; //Define a white gaussian noise with Euler-Maruyama method of simulation
    myField field1(0.5,14.0,8.0,0.2,&noise1); //Define the parameters of the field above and pass the kind of noise
    
    //Now we can finally define our system as 3D one (first arg), with the created field (second arg) and we tell also that it is bounded (third arg).
    SDE_SS_System system1(3,&field1,true);

    //Having defined the system as "bounded" we need also to pass the bound function
    system1.setBoundFunction(boundaries);

    //Produce the trajectory of the first system
    Traj traj1{system1.simulateTrajectory({0.7,0.0,0.0},2250,0.01)};

    //Prepare the element of the second system
    WienerEuler noise2;
    myDataLinker linker(traj1,false);
    DataField field2(1.0,&linker,&noise2);

    //We can now create the second system
    SDE_SS_System system2(2,&field2,false);

    //Let's make a trajectory for the second system
    Traj traj2{system2.simulateTrajectory({0.0,0.0},100,0.01)};

    //THIS LAST PART IS USED SIMPLY TO SAVE THE RESULTS.
    ofstream FILE("traj_data.dat");

    for(size_t i=0;i<traj2.getLength();i++){
            FILE << traj2.getVars()[i][0] << " " << traj2.getVars()[i][1] << endl; 
    }

    FILE.close();

    cout << "FINISH" << endl;

    return 0;
}