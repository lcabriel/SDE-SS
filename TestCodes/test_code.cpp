#include "SDE-SS.h" //The library
#include <fstream> //Used to write the file

using namespace std;

//#######################################################################################################################################################################
/*1. We need at first to create our field. To do this we need to create a custom version of the Field class. Defining as variable all the parameters
we are going to use.

2. In the construction we need also to pass the noise. The class is made to accept any NoiseClass object however NoiseClass is purely virtual.
You can create one of your own looking to the standard ones defined in "SDE-SS.cpp" however the Gaussian White noise is already implemented both
using the Euler-Maruyama method and the Milstein one. They are calle "WienerEuler" and "WienerMilstein". 
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
    void f_function_impl(const valarray<float> &x,float t,valarray<float> &y) const override{

        y[0] = x[0]*(1-x[0])*(x_e - x[1]);
        y[1] = (x[2]-x[1])/tau2;
        y[2] = (x[0]-x[2])/tau1;
    }

    //The same can be said for the g_function.
    void g_function_impl(const valarray<float> &x,float t,valarray<float> &y) const override{

        y[0] = -x[0]*(1-x[0])*D;
        y[1] = 0;
        y[2] = 0;
    }

};


//The system described by the equations above it is also a bounded one in [0,1]. To this we need to define a bound function to pass to
//our system. The bound function expression is standard: it returns a boolean and requires a vector<float> as argument.
bool boundaries(valarray<float> x){
    return ((x[0]>=0)&&(x[0]<=1));
}

//When generating PDFs (e.g. using the functions PDF_1D), the initial conditions could be fixed or chosen according to a function.
//If we want to use the custom function we need to create a void -> vector<float> function to tell the system how the initial condition should be picked.
vector<float> initial_cond_function(){
    vector<float> x0(3,0.0);

    //For our problem we are going to set the variable 1 and 2 equal every time to 0
    //while the 0 variable will be chosen randomly with a uniform distribution in the system domain [0,1]

    //Define the engine and the distribution
    uintptr_t seed = reinterpret_cast<std::uintptr_t>(&x0); //Use the pointer as seed
    mt19937 eng =  mt19937(static_cast<long unsigned int>(time(0)*seed));
    uniform_real_distribution<double> distribution = uniform_real_distribution<double>(0,1);

    x0[0] = distribution(eng);

    return x0;
}

//####################################################################### MAIN ##########################################################################################
/*
In the main, we are going to show use and example of usage producing the PDF along the first variable for our system.
*/

int main(){
    WienerEuler noise; //Define a white gaussian noise with Euler-Maruyama method of simulation
    myField field(0.5,14.0,8.0,0.2,&noise); //Define the parameters of the field above and pass the kind of noise
    
    //Now we can finally define our system as 3D one (first arg), with the created field (second arg) and we tell also that it is bounded (third arg).
    SDE_SS_System system(3,&field,true);

    //Having defined the system as "bounded" we need also to pass the bound function
    system.setBoundFunction(boundaries);

    //Computing the PDF is heavy so we increase the number of threads used by heavy operations from 8 to 10
    system.setNumThreads(10);

    Traj T = system.simulateTrajectory({0.5,0.0,0.0},2000,0.01);

    /*We want our PDF in the instant equal to the last step thus we need to create first a TimePicture of the system.
    The points in the TimePicture are stored in a 2D array with along the row the trajectories and along the columns
    the variables. The values are taken at the time instant where the TimePicture is taken. 
    Wanting to have the last step for our PDF, we need only to specify certain arguments of this function:
    - 1°: it is the length of our simulations
    - 2°: the step size of the simulations
    - 3°: the number of simulations
    - 4°: are the initial condition not fixed?
    - 5°: being not fixed, we need the initial condition functions. If fixed an array x0 is required
    - 6° (NOT VISIBLE): is the time instant, if not specified it is the last step of the trajectories.
    */
    TimePicture picture=system.produceTimePicture(2000,0.01,50,true,{{}},initial_cond_function);

    //Created the time picture we can transform it into a 1D bin. We need to pass it to the PDF_1D function with:
    //2°: the number of bins
    //3:  the axis along which we are binning (0: time, n: the nth variable of the system)
    //4: the boolean is to have an auto adaptive domain. If true, the min and the max of the values in the picture
    //are used as extremes and the Nbins constructed between them.
    //5: Having chose a fixed domain we need to pass the domain as a vector.
    vector<vector<float>> bins = system.PDF_1D(picture,50,1,false,{0.0,1.0});

    //THIS LAST PART IS USED SIMPLY TO SAVE THE RESULTS.
    ofstream FILE("test_bins.dat");

    for(unsigned long int i=0;i<bins.size();i++){
            FILE << bins[i][0] << " " << bins[i][1] << endl; 
    }

    FILE.close();

    cout << "FINISH" << endl;

    return 0;
}