#ifndef SDE_SS_H
#define SDE_SS_H

#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <ctime>
#include <cstdlib>
#include <functional>
#include <stdexcept>
#include <cstdint>
#include <omp.h>

using namespace std;

//Function to produce in an easier way the error when needed.
inline void error(const string& s){
    throw runtime_error(s);
}

//#################################################### GENERIC NOISE CLASS ###########################################################

//Generic and virtual noise class: can be used to create new noises. compute_noise is used by the library to insert the noise in the SDE.
class NoiseClass{

public:
    
    //A standard compute noise that simply gives back the passed point in the system space.
    virtual vector<float> compute_noise(const vector<float> &x_i,float* h);    

};

//############################################### PREDEFINED WIENER PROCESS WITH EULER-MARUYAMA NOISE ###########################################################

//A useful premade class to produce values of a Wiener process simulating dW with an Euler-Maruyama method.
class WienerEuler: public NoiseClass{

    mt19937 eng; //Random number generator engine.
    normal_distribution<double> distribution; //Normal standard distribution.
    
public:

    //COSTRUCTOR: the distribution and the generator are initialized.
    //The seed of the random engine is created starting from the time and the point of the instance.
    WienerEuler();

    //Compute a dW kind of step for the Wiener process using the Euler-Maruyama method. 
    //Actually a vector of the same size of the system is returned and, in each slot, there is a different dW. 
    vector<float> compute_noise(const vector<float> &x_i,float* h) override;

};

//############################################### PREDEFINED WIENER PROCESS WITH MILSTEIN NOISE ################################################

//A useful premade class to produce values of a Wiener process simulating dW with a Milstein method.
class WienerMilstein: public NoiseClass{

    mt19937 eng; //Random number generator engine.
    normal_distribution<double> distribution; //Normal standard distribution.
    function<vector<float>(vector<float>)> D_g; //The derivative of the field g_function.

public:

    //CONTRUCTOR: the distribution and the generator are initialize. 
    //Moreover, requires as argument the derivative of the g_function used in the field class.
    //The seed of the random engine is created starting from the time and the point of the instance.
    WienerMilstein(function<vector<float>(const vector<float>&)> derivative = nullptr);

    //Compute a dW kind-of step for the Wiener process using the Milstein method. 
    //Actually a vector of the same size of the system is returned and, at each step there is a different dW. 
    vector<float> compute_noise(const vector<float> &x_i,float* h) override;

};

//################################################## GENERIC FIELD CLASS ###################################################################

//A generic FieldClass to express the field of SDE. Requires a NoiseClass to indentify the type of noise
//for which g_function is multiplied in the different differential equations of the system. 
//To perform some controls, the user has to define the system as in the examples overriding the _impl versions
//of the functions.
class FieldClass{

    NoiseClass* noise; //NoiseClass to generate the noise.
    bool noise_initialized{false}; //Flag to keep track if the noise is initialized.

public:

    //This public function allow to set the noise of the FieldClass in a controlled way.
    void setNoise(NoiseClass* N);

    //Function for the deterministic part of the field. To define, please
    //override "f_function_impl" as in documentation.
    vector<float> f_function(const vector<float> &x);

    //Function for the stochastic part of the field. To define, please
    //override "g_function_impl" as in documentation.
    vector<float> g_function(const vector<float> &x);

    //Function for the deterministic part of the field. Override this as in
    //documentation to implement your system.
    virtual vector<float> f_function_impl(const vector<float> &x);

    //Function for the deterministic part of the field. Override this as in
    //documentation to implement your system.
    virtual vector<float> g_function_impl(const vector<float> &x);

    //This function will give the result of the compute_noise of the local NoiseClass.
    vector<float> getNoise(const vector<float> &x_i,float* h);

};

//################################################# SYSTEM #################################################################################

//The main class of the library. It is used to define a system of SDEs and to produce trajectories of that.
//The class contain also some useful tool functions based on the trajectories.
class SDE_SS_System{

    unsigned int size; //Size of the system/problem.
    FieldClass* field; //Field of the system.
    bool bounded{false}; //Is the system bounded?
    function<bool(vector<float>)> bounds; //the function used to express the bounds of the domain, if present.
    unsigned int NumThreads{8}; //The number of threads used for the parallelizable operations.

    //################### SUPPORT FUNCTIONS #######################
    //These functions are used to light the code of other more complex functions.

    //This function will perform the checks of the parameters of computeTimePicture.
    void checkFunctionComputeTimePicture(unsigned int Nsim,bool random_initial,const vector<float> &x0,function<vector<float>()> random_f);

    //This function will perform the checks of the parameters of simulateTrajectory. 
    void checkTrajInput(const vector<float> &x0,float period, float h_0);

    //Given a vector of times and a time index this function will find the slot just before the given time instant.
    //This function is also STATIC.
    static size_t findTimeIndex(const vector<float> &times,float TI);

    //Given a vector of the shape of the one of the simulateTrajectory function,
    //this function will extract the column of times (column 0). This function is also STATIC.
    static vector<float> extractTimes(const vector<vector<float>> &traj);

    //This function will perform the checks of the parameters of PDF_1D.
    void checkFunctionPDF_1D(const vector<vector<float>> &picture,unsigned int Nbins,unsigned int axis,bool adaptive,const vector<float> &domain);

    //This function will perform the checks of the parameters of PDF_2D.
    void checkFunctionPDF_2D(const vector<vector<float>> &picture,vector<unsigned int> Nbins,vector<unsigned int> axis,bool adaptive,const vector<float> &domain);

    //This function will perform the checks of the parameters of computeAutocorrelation.
    void checkAutocorrelationInput(const vector<vector<float>> &traj,unsigned int axis,float tau);

    //################### EVOLUTION FUNCTIONS ####################
    //These functions are used to perform the evolutive aspects of the trajectory's computation.

    //Given the previous point and the step length, this internal function is the core function to evolve 
    //the last step in the new one of the trajectory. It will use the RK4 method.
    //The idea is to use a setup in the additive splitting way synergizing with RK4 and the noise method.
    vector<float> evolveTraj(const vector<float> &x_n,float h);

    //Given the starting point and the step, this function will return the RK4 update
    //of the deterministic part of the field.
    //REMEMBER: you have to multiply externally by h or your step eventually.
    vector<float> RK4_method(const vector<float> &x0,float h);

public:

    //Constructor: requires the size of the system, the Field Class of the system
    //and a bool to say if the system is bounded or not (rem: you need to set the "BoundFunction"
    //with the specificy method of this class, in this case).
    SDE_SS_System(unsigned int N,FieldClass* F,bool isBounded=false); 

    //####################### CORE FUNCTIONS ###################################

    //This core function simulate a single trajectory given the initial conditions, the time period and the standard step.
    //This function will return a 2D vector with all the points of the trajectory. 
    //In the returning vector the times and the values are conjointed. The first column is used for times.
    vector<vector<float>> simulateTrajectory(const vector<float> &x0,float period, float h_0);

    //###################### PUBLIC UTILITY FUNCTIONS ##########################

    //This function will set the bound function. The bound function should return a boolean
    //and has as input a vector<float> (the point). The function should return "true" when
    //the point is in the domain.
    //The system has to be bounded (construction) to this function to work.
    void setBoundFunction(function<bool(const vector<float>&)> f);

    //This function is used to set the number of threads used by the heavy functions of the class.
    //The standard value is 8.
    void setNumThreads(unsigned int N);

    //####################### TOOL FUNCTIONS #####################################
    //These function are made to automatize some typical tests and procedure used in the analysis of the SDE.

    //This function will automatically produce a group of the trajectories to obtain the value in a time instant. 
    //It requires, in order,:
    //- The period, the step size and the number of trajectories (MANDATORY).
    //- A bool to say if the initial condition are random (DEFAULT = false = non random initial conditions).
    //- If the initial conditions are not random, a valid initial condition vector is necessary.
    //- If they are random, a function to compute the initial condition is required. This function has to be
    //  a zero argument - returning vector<float> function.
    //- Eventually, a time index where the PDF has to be computed. If it is not given, the PDF will be computed
    //  at the last step.
    //The output of the function will be a 2D matrix with the values of each trajectory at the time instant in each row.   
    vector<vector<float>> produceTimePicture(float period,float h_0,unsigned int Nsim,
                                            bool random_initial=false,const vector<float> &x0 = {},
                                            function<vector<float>()> random_f = nullptr,float time_instant = -1.0f);

    //This function will produce starting from a TimePicture such the one produced by "produceTimePicture" a 1D bin
    //system useful to obtain PDFs. This function require:
    //- The TimePicture style vector<vector<float>> (MANDATORY).
    //- The number of bins (MANDATORY).
    //- The axis along which doing the bins (0: times,n: the n variable of the system) (MANDATORY).
    //- A bool to express if the binning domain is given or adaptive (true = adaptive).
    //- If false a vector with upper and lower domain is required.
    //The function will return a 2D vector with in each row the central value (first column) and the bin value (second column).
    vector<vector<float>> PDF_1D(const vector<vector<float>> &picture,unsigned int Nbins,unsigned int axis,
                                bool adaptive = false,vector<float> domain = {0.0,0.0}); 
    
    //This function will produce starting from a Time Picutre such the one produced by "produceTimePicture" a 2D bin
    //system useful to obtain PDFs, This function require:
    //- The TimePicture style vector<vector<float>> (MANDATORY).
    //- A 2D vector for the number of bins along the axis (MANDATORY).
    //- A 2D vector of the axis along which doing the bins (0: times,n: the n variable of the system) (MANDATORY).
    //- A bool to express if the binning domain is given or adaptive (true = adaptive).
    //- If false a vector with upper and lower domain is required. It should be a 4 slot vector (lb,ub,lb,ub).
    //The function will return a 2D vector with in each row the central value coordinates (firsts 2 column) and the bin value (last column).
    vector<vector<float>> PDF_2D(const vector<vector<float>> &picture,vector<unsigned int> Nbins,vector<unsigned int> axis,
                                bool adaptive = false, vector<float> domain = {0.0,0.0,0.0,0.0});

    //This function will produce, starting from a trajectory such the one produced by "produceTimePicture" the
    //autocorrelation of the trajectory for a certain time delay. This function require:
    //- A trajectory in style vector<vector<float>> (MANDATORY).
    //- The axis (variable) along which computing the autocorrelation (MANDATORY).
    //- The time delay (\tau) of the autocorrelation (MANDATORY).
    //The output will be the autocorrelation value computed as covariance/variance.
    //Also the time delay is converted in the nearest below number of steps.
    float computeAutocorrelation(const vector<vector<float>> &traj,unsigned int axis,float tau);
    
};

#endif