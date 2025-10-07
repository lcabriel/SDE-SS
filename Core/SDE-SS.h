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
#include <algorithm>
#include <valarray>

using namespace std;

//Function to produce in an easier way the error when needed.
[[noreturn]] inline void error(const string& msg){
    cerr << "Fatal error: " << msg << endl;
    exit(EXIT_FAILURE);
}

//#################################################### GENERIC NOISE CLASS ###########################################################

//Generic and virtual noise class: can be used to create new noises. compute_noise is used by the library to insert the noise in the SDE.
class NoiseClass{

public:
    
    //A standard compute noise that simply gives back the passed point in the system space.
    virtual valarray<float> compute_noise(const valarray<float> &x_i,const float* h);    

};

//############################################### PREDEFINED WIENER PROCESS WITH EULER-MARUYAMA NOISE ###########################################################

//A useful premade class to produce values of a Wiener process simulating dW with an Euler-Maruyama method.
class WienerEuler: public NoiseClass{

    mt19937 eng; //Random number generator engine.
    normal_distribution<double> distribution; //Normal standard distribution.
    
public:

    //CONSTRUCTOR: initializes the distribution and random generator.
    //The seed of the random engine is created starting from the time and the point of the instance.
    WienerEuler();

    //Compute a dW kind of step for the Wiener process using the Euler-Maruyama method. 
    //Actually a vector of the same size of the system is returned and, in each slot, there is a different dW. 
    valarray<float> compute_noise(const valarray<float> &x_i, const float* h) override;

};

//############################################### PREDEFINED WIENER PROCESS WITH MILSTEIN NOISE ################################################

//A useful premade class to produce values of a Wiener process simulating dW with a Milstein method.
class WienerMilstein: public NoiseClass{

    mt19937 eng; //Random number generator engine.
    normal_distribution<double> distribution; //Normal standard distribution.
    function<valarray<float>(valarray<float>)> D_g; //The derivative of the field g_function.

public:

    //CONSTRUCTOR: the distribution and the generator are initialize. 
    //Moreover, requires as argument the derivative of the g_function used in the field class.
    //The seed of the random engine is created starting from the time and the point of the instance.
    WienerMilstein(const function<valarray<float>(const valarray<float>&)>& derivative = nullptr);

    //Compute a dW kind-of step for the Wiener process using the Milstein method. 
    //Actually a vector of the same size of the system is returned and, at each step there is a different dW. 
    valarray<float> compute_noise(const valarray<float> &x_i, const float* h) override;

};

//################################################## GENERIC FIELD CLASS ###################################################################

//A generic FieldClass to express the field of SDE. Requires a NoiseClass to indentify the type of noise
//for which g_function is multiplied in the different differential equations of the system. 
//To perform some controls, the user has to define the system as in the examples overriding the _impl versions
//of the functions.
class FieldClass{

    NoiseClass* noise; //NoiseClass to generate the noise.
    bool noise_initialized{false}; //Flag to keep track if the noise is initialized.

protected:

    //This public function allow to set the noise of the FieldClass in a controlled way.
    void setNoise(NoiseClass* const N);

    //Function for the deterministic part of the field. To define, please
    //override "f_function_impl" as in documentation.
    void f_function(const valarray<float> &x,float t,valarray<float> &y);

    //Function for the stochastic part of the field. To define, please
    //override "g_function_impl" as in documentation.
    void g_function(const valarray<float> &x,float t,valarray<float> &y);

    friend class SDE_SS_System;

public:

    //Implementation function for the deterministic part of the field. Override this as in
    //documentation to implement your system.
    virtual void f_function_impl(const valarray<float> &x,float t,valarray<float> &y) const;

    //Implementation function for the stochastic part of the field. Override this as in
    //documentation to implement your system.
    virtual void g_function_impl(const valarray<float> &x,float t,valarray<float> &y) const;

    //This function will give the result of the compute_noise of the local NoiseClass.
    valarray<float> getNoise(const valarray<float> &x_i,const float* h);

};

//################################################# TRAJECTORY #############################################################################

//The Traj class is used to represent a simulated trajectory of a system of SDEs. It is purely an object class to manage easily the values of
//the variables and the time instants keeping them all in the same place. It is the product of the simulateTrajectory function.
class Traj{

    vector<float> times; //The vector of the time instant of the trajectory
    vector<vector<float>> vars; //The 2D array of the system variables. Rows: time; Columns: variables.
    size_t step_num; //The time length/num of steps of the trajectory.

public:

    //CONSTRUCTOR 1: the time vector and the 2D array of variable values using std::move.
    Traj(vector<float>&& t,vector<vector<float>>&& v);

    //CONSTRUCTOR 2: passing another Traj instance.
    Traj(const Traj& t);

    //####################### GET FUNCTIONS ##############################

    //Return the reference to the vector of the time instants considered in the trajectory.
    const vector<float>& getTimes() const{
        return times;
    }

    //Return the reference to the 2D array of the variables in every time instant of the trajectory.
    const vector<vector<float>>& getVars() const{
        return vars;
    }

    //Return the number of steps/length of the trajectory.
    const size_t& getLength() const{
        return step_num;
    }

    //Given a certain time index, this function will return the situation as a vector of times+vars. 
    vector<float> getInstant(const size_t index) const;

};

//################################################# TIME PICTURE ###########################################################################

//This TimePicture class is used to represent a so-called time picture: a set of N points taken from N simulated trajectory of a SDEs system
//at a given time instant T. As Traj, this is purely an object class to manage easily the values of the trajectories at that time instant.
//It is the product of the computeTimePicture function.
class TimePicture{

    vector<vector<float>> points; //The values of the N simulations at the time instant T. (row,columns)=(trajs,variables).
    size_t N; //The number of simulations.
    float T; //The time instant.

public:

    //CONSTRUCTOR 1: build using the 2D array of values and the time instant. The set of values must have the different trajectories
    //along the rows and the variables along the columns.
    TimePicture(vector<vector<float>> values,float t);

    //CONSTURCTOR 2: build copying another TimePicture instance.
    TimePicture(const TimePicture& TP);

    //####################### GET FUNCTIONS ##############################

    //Return the reference to the 2D array of the values of the different trajectories at the time instant of the picture.   
    const vector<vector<float>>& getAllPoints() const{
        return points;
    }

    //Return the reference to the array of the values of a specific trajectory "p" at the time instant of the picture. 
    const vector<float>& getPoint(size_t p) const{
        return points[p];
    }

    //Return the reference to the total number of trajectories composing the TimePicture.
    const size_t& getNumSim() const{
        return N;
    }

    //Return the reference to the value identifying the time instant of the TimePicture.  
    const float& getTimeInstant() const{
        return T;
    }

};

//################################################# DATALINKER ############################################################################

//The generic DataLinker is an element that can be used by the FieldClass when the field presents a forcing factor which is given by a data set (vector<vector<float>>) previously produced.
//The function to obtain get data values is "getData" which should be override in a child class as explained in the documentation.
//Moreover, the class is kept generic to add in the child possible variations of "getData" or if some parameters are required.
class DataLinker{

protected:

    Traj data; //The array of data used in the simulation.
    const bool adapting_set; //Used to express if the set can change according to the variables of the system where the DataLinker is used.
    size_t TOC{0}; //TimeOptimizationCounter: used for optimization, indicate the position on data of the previously picked time instant.

    //Given a certain time instant, the variable will find the a time-index near the value (the next one). The TimeOptimizationCounter
    //is used to improve performances allowing the bypass of this function when multiple getData access to the sime time index (see examples).
    size_t findTimeIndex(const float t);

    //This function is made to be used inside of getData. It allows, instead of using the time index found with findTimeIndex, to obtain a sort
    //in-between value in the dataset using a linear interpolator between the previous and the next point along the given variable.
    //WARNING: this function could increase significatly your computation times and it is made for dataset characterised by a wide grid.
    float interpolData(const float& t,const unsigned int& variable);

public:

    //CONSTRUCTOR 1: passing a traj and specifying if is adaptive
    DataLinker(const Traj& t,const bool AS);

    //CONSTRUCTOR 2: copying a DataLinker
    DataLinker(const DataLinker& DL);

    //Given a certain time instant (and sys vars for adapting sets), the getData returns a child-defined float value near to the given time.
    virtual float getData(const float t,const valarray<float>& x);

    //This function allow, if the DataLinker is adaptive, to pass a new set of data to the DataLinker.
    void setNewData(const Traj& t);

    //This function is used to reset the TimeOptimizationCounter which is used to optimize the finding process of the time instant during the
    //simulations. 
    void resetTimeOptimizationCounter();

};

//################################################# SET OF POINTS ##########################################################################

//The SetOfPoints is conceptually very similar to the trajectory. It shares a lot of the characteristics and public functions but instead
//of having all the points of the trajectory it contains only a set of it that are meaningful for the user. It is the product of
//simulateTrajectorySOP function.
class SetOfPoints{

    vector<float> times; //The vector of the time instants of the set
    vector<vector<float>> vars; //The 2D array of the system variables. Rows: time instants; Columns: variables.

public:

    //CONSTRUCTOR 1: the time instants vector and the 2D array of variable values using std::move.
    SetOfPoints(vector<float>&& t,vector<vector<float>>&& v);

    //CONSTRUCTOR 2: passing another SetOfPoints instance.
    SetOfPoints(const SetOfPoints& sop);

    //####################### GET FUNCTIONS ##############################

    //Return the reference to the vector of the time instants considered in the set.
    const vector<float>& getTimes() const{
        return times;
    }

    //Return the reference to the 2D array of the variables in every time instant of the set.
    const vector<vector<float>>& getVars() const{
        return vars;
    }

    //Given a certain index, this function will return the situation as a vector of times+vars 
    vector<float> getInstant(const size_t index) const;

};

//################################################# SYSTEM #################################################################################

//The main class of the library. It is used to define a system of SDEs and to produce trajectories of that.
//The class contain also some useful tool functions based on the trajectories.
class SDE_SS_System{

    unsigned int size; //Size of the system/problem.
    FieldClass* field; //Field of the system.
    bool bounded{false}; //Is the system bounded?
    function<bool(valarray<float>)> bounds; //the function used to express the bounds of the domain, if present.
    unsigned int NumThreads{8}; //The number of threads used for the parallelizable operations.

    //################### SUPPORT FUNCTIONS #######################
    //Helper functions to simplify the logic of higher-level routines.

    //This function will perform the checks of the parameters of computeTimePicture.
    void checkFunctionComputeTimePicture(const unsigned int Nsim, const bool random_initial,const vector<vector<float>> &x0,
                                        const function<vector<float>()>& random_f);

    //This function will perform the checks of the parameters of simulateTrajectory. 
    void checkTrajInput(const vector<float> &x0,const float period,const float h_0);

    //This function will perform the checks of the parameters of PDF_1D.
    void checkFunctionPDF_1D(const unsigned int Nbins,const unsigned int axis,const bool adaptive,const vector<float> &domain);

    //This function will perform the checks of the parameters of PDF_2D.
    void checkFunctionPDF_2D(const vector<unsigned int> Nbins,const vector<unsigned int> axis,const bool adaptive,const vector<float> &domain);

    //This function will perform the checks of the parameters of computeAutocorrelation.
    void checkAutocorrelationInput(const Traj& traj,const unsigned int axis,const float tau);

    //################### EVOLUTION FUNCTIONS ####################
    //These functions are used to perform the evolutive aspects of the trajectory's computation.

    //Given the previous point and the step length, this internal function is the core function to evolve 
    //the last step in the new one of the trajectory. It will use the RK4 method.
    //The idea is to use a setup in the additive splitting way synergizing with RK4 and the noise method.
    valarray<float> evolveTraj(const valarray<float> &x_n,float h,float t,vector<valarray<float>> &k,
        valarray<float> &x1,valarray<float> &g,valarray<float> &n,valarray<float> &x_temp);

    //Given the starting point and the step, this function will return the RK4 update
    //of the deterministic part of the field.
    //REMEMBER: you have to multiply externally by h or your step eventually.
    valarray<float> RK4_method(const valarray<float> &x0,float h,float t,vector<valarray<float>> &k,
        valarray<float> &x_temp);

public:

    //Constructor: requires the size of the system, the Field Class of the system
    //and a bool to say if the system is bounded or not (rem: you need to set the "BoundFunction"
    //with the specificy method of this class, in this case).
    SDE_SS_System(unsigned int N,FieldClass* F,bool isBounded=false); 

    //####################### CORE FUNCTIONS ###################################

    //This core function simulate a single trajectory given the initial conditions, the time period and the standard step.
    //This function will return an element of Traj class.
    Traj simulateTrajectory(const vector<float> &x0,const float period,const float h_0);

    //This core function acts as simulateTrajectory, however the internal points of the trajectory are not saved and only
    //the last value of the trajectory is returned as a vector of shape (time,[coords]).
    vector<float> simulateTrajectoryLastPoint(const vector<float> &x0,const float period,const float h_0);

    //This core function acts as simulateTrajectory but it will not return the entire trajectory. In fact, this function
    //asks for a set of time instants and the output will be the values of the trajectory in those instants (actually the
    //immediately before point).
    SetOfPoints simulateTrajectorySOP(const vector<float> &x0,const float period,const float h_0,const vector<float> &instants);

    //This function will automatically produce a group of the trajectories to obtain the value in a time instant. 
    //It requires, in order,:
    //- The period, the step size and the number of trajectories (MANDATORY).
    //- A bool to say if the initial condition are random (DEFAULT = false = non random initial conditions).
    //- If the initial conditions are not random, a vector of  valid initial condition vectors is necessary.
    //  It should have as many initial condition points as Nsim.
    //- If they are random, a function to compute the initial condition is required. This function has to be
    //  a zero argument - returning vector<float> function.
    //- Eventually, a time index where the PDF has to be computed. If it is not given, the PDF will be computed
    //  at the last step.
    //The output of the function will be a 2D matrix with the values of each trajectory at the time instant in each row.     
    TimePicture produceTimePicture(const float period,const float h_0,unsigned int Nsim,
                                            const bool random_initial=false,const vector<vector<float>> &x0 = {{}},
                                            const function<vector<float>()>& random_f = nullptr,const float time_instant = -1.0f);

    //###################### PUBLIC UTILITY FUNCTIONS ##########################

    //This function will set the bound function. The bound function should return a boolean
    //and has as input a vector<float> (the point). The function should return "true" when
    //the point is in the domain.
    //The system has to be bounded (construction) to this function to work.
    void setBoundFunction(const function<bool(const valarray<float>&)>& f);

    //This function is used to set the number of threads used by the heavy functions of the class.
    //The standard value is 8.
    void setNumThreads(unsigned int N);

    //Given a vector of times and a time index this function will find the slot just before the given time instant.
    //This function is also STATIC.
    static size_t findTimeIndex(const vector<float> &times,const float TI);

    //####################### TOOL FUNCTIONS #####################################
    //These function are made to automatize some typical tests and procedure used in the analysis of the SDE.

    //This function will produce starting from a TimePicture such the one produced by "produceTimePicture" a 1D bin
    //system useful to obtain PDFs. This function require:
    //- A TimePicture (MANDATORY).
    //- The number of bins (MANDATORY).
    //- The axis along which doing the bins (0: times,n: the n variable of the system) (MANDATORY).
    //- A bool to express if the binning domain is given or adaptive (true = adaptive).
    //- If false a vector with upper and lower domain is required.
    //The function will return a 2D vector with in each row the central value (first column) and the bin value (second column).
    vector<vector<float>> PDF_1D(const TimePicture& picture,unsigned int Nbins,unsigned int axis,
                                bool adaptive = false,vector<float> domain = {0.0,0.0}); 
    
    //This function will produce starting from a Time Picutre such the one produced by "produceTimePicture" a 2D bin
    //system useful to obtain PDFs, This function require:
    //- A TimePicture (MANDATORY).
    //- A 2D vector for the number of bins along the axis (MANDATORY).
    //- A 2D vector of the axis along which doing the bins (0: times,n: the n variable of the system) (MANDATORY).
    //- A bool to express if the binning domain is given or adaptive (true = adaptive).
    //- If false a vector with upper and lower domain is required. It should be a 4 slot vector (lb,ub,lb,ub).
    //The function will return a 2D vector with in each row the central value coordinates (firsts 2 column) and the bin value (last column).
    vector<vector<float>> PDF_2D(const TimePicture& picture,vector<unsigned int> Nbins,vector<unsigned int> axis,
                                bool adaptive = false, vector<float> domain = {0.0,0.0,0.0,0.0});

    //This function will produce, starting from a trajectory such the one produced by "produceTimePicture" the
    //autocorrelation of the trajectory for a certain time delay. This function require:
    //- A trajectory in style vector<vector<float>> (MANDATORY).
    //- The axis (variable) along which computing the autocorrelation (MANDATORY).
    //- The time delay (tau) of the autocorrelation (MANDATORY).
    //The output will be the autocorrelation value computed as covariance/variance.
    //Also the time delay is converted in the nearest below number of steps.
    float computeAutocorrelation(const Traj& traj,unsigned int axis,float tau);
    
};

#endif