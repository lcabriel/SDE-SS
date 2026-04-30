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
#include <cstring>

using namespace std;

//Function to produce in an easier way the error when needed.
[[noreturn]] inline void error(const string& msg){
    cerr << "Fatal error: " << msg << endl;
    exit(EXIT_FAILURE);
}

//#################################################### RANDOM NUMBER GENERATOR #######################################################
//In this section we implement some class and methods that are used to generate random numbers in our library. In particular we opted
//a PCG32 random engine and for the gaussian distribution a Ziggurat Algorithm.

/*
 High-performance Random Number Generator (Permuted Congruential Generator).
 Based on the algorithm by Melissa O'Neill, PCG32 is specifically optimized for SDE simulations:
 - MINIMAL STATE: It requires only 16 bytes of memory (compared to 2.5KB for mt19937),
   ensuring excellent L1 cache locality, which is critical in multi-threaded (OpenMP) contexts.
 - XSH-RR TRANSFORMATION: It applies a fixed bit-shift and a dynamic "Random Rotate"
   to the 64-bit internal LCG state. This eliminates the linear correlations typical
   of standard LCGs and produces statistically robust 32-bit outputs.
*/

class PCG32 {
private:
    uint64_t state;
    uint64_t inc;

    // Ziggurat tables for the Standard Normal Distribution (128 layers)
    // Parameters: R = 3.442619855899, Area = 0.00991256303526217
    static const uint32_t zig_kn[128];
    static const float zig_wn[128];
    static const float zig_fn[128];

public:
    // Constructor: seed sets the initial state, seq allows different sequences
    PCG32(uint64_t seed = 0x4d595df4d0f33173ULL, uint64_t seq = 1u) {
        state = 0U;
        inc = (seq << 1u) | 1u;
        next();
        state += seed;
        next();
    }

    // Standard PCG32 output transformation (XSH-RR)
    uint32_t next() {
        uint64_t oldstate = state;
        // Linear Congruential Generator step
        state = oldstate * 6364136223846793005ULL + inc;
        // Output transformation: Xorshift followed by Random Rotate
        uint32_t xorshifted = static_cast<uint32_t>(((oldstate >> 18u) ^ oldstate) >> 27u);
        uint32_t rot = static_cast<uint32_t>(oldstate >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
    }

    // Generates a random float in the range [0, 1)
    float nextFloat() {
        return (next() >> 8) * 0x1.0p-24f;
    }

    /**
     * Generates a random number from a standard normal distribution N(0,1)
     * using the Ziggurat Algorithm.
     * * This implementation uses a 128-layer Ziggurat with signed integer logic
     * to provide a fast, symmetric distribution.
     */
    float nextGaussian() {
        for (;;) {
            uint32_t k = next();
            
            // GSL logic: i is the first 8 bits
            uint32_t i_full = (k & 0xFF);
            // j is the next 24 bits (from bit 8 to 31)
            uint32_t j = (k >> 8) & 0xFFFFFF; 

            // Extract sign from the 8th bit of i_full
            float sign = (i_full & 0x80) ? 1.0f : -1.0f;
            // The actual layer index is the first 7 bits (0-127)
            uint32_t i = i_full & 0x7F;

            // STEP 1: Fast Path
            // x is calculated using 24-bit resolution j
            float x = j * zig_wn[i];

            if (j < zig_kn[i]) {
                return sign * x;
            }

            // STEP 2: Boundary and Tail
            float y;
            if (i < 127) {
                // Wedge sampling
                float y0 = zig_fn[i];
                float y1 = zig_fn[i + 1];
                y = y1 + (y0 - y1) * nextFloat();
            } else {
                // Tail sampling (Layer 127 in GSL)
                float U1 = 1.0f - nextFloat();
                float U2 = nextFloat();
                // PARAM_R for 128 steps is 3.442619855899
                x = 3.442619855899f - std::log(U1) / 3.442619855899f;
                y = std::exp(-3.442619855899f * (x - 0.5f * 3.442619855899f)) * U2;
            }

            // Final Rejection Test
            if (y < std::exp(-0.5f * x * x)) {
                return sign * x;
            }
        }
    }
};

// --- Verified Table Initializations ---

inline const uint32_t PCG32::zig_kn[128] = {
    0, 12590644, 14272653, 14988939, 15384584, 15635009, 15807561, 15933577, 16029594, 16105155, 16166147, 16216399, 16258508, 16294295, 16325078, 16351831,
    16375291, 16396026, 16414479, 16431002, 16445880, 16459343, 16471578, 16482744, 16492970, 16502368, 16511031, 16519039, 16526459, 16533352, 16539769, 
    16545755, 16551348, 16556584, 16561493, 16566101, 16570433, 16574511, 16578353, 16581977, 16585398, 16588629, 16591685, 16594575, 16597311, 16599901, 
    16602354, 16604679, 16606881, 16608968, 16610945, 16612818, 16614592, 16616272, 16617861, 16619363, 16620782, 16622121, 16623383, 16624570, 16625685, 
    16626730, 16627708, 16628619, 16629465, 16630248, 16630969, 16631628, 16632228, 16632768, 16633248, 16633671, 16634034, 16634340, 16634586, 16634774,
    16634903, 16634972, 16634980, 16634926, 16634810, 16634628, 16634381, 16634066, 16633680, 16633222, 16632688, 16632075, 16631380, 16630598, 16629726, 
    16628757, 16627686, 16626507, 16625212, 16623794, 16622243, 16620548, 16618698, 16616679, 16614476, 16612071, 16609444, 16606571, 16603425, 16599973, 
    16596178, 16591995, 16587369, 16582237, 16576520, 16570120, 16562917, 16554758, 16545450, 16534739, 16522287, 16507638, 16490152, 16468907, 16442518, 
    16408804, 16364095, 16301683, 16207738, 16047994, 15704248, 15472926
};

inline const float PCG32::zig_wn[128] = {
    1.62318314817e-08, 2.16291505214e-08, 2.54246305087e-08, 2.84579525938e-08, 3.10340022482e-08, 3.33011726243e-08, 3.53439060345e-08, 3.72152672658e-08,
    3.8950989572e-08, 4.05763964764e-08, 4.21101548915e-08, 4.35664624904e-08, 4.49563968336e-08, 4.62887864029e-08, 4.75707945735e-08, 4.88083237257e-08,
    5.00063025384e-08, 5.11688950428e-08, 5.22996558616e-08, 5.34016475624e-08, 5.44775307871e-08, 5.55296344581e-08, 5.65600111659e-08, 5.75704813695e-08,
    5.85626690412e-08, 5.95380306862e-08, 6.04978791776e-08, 6.14434034901e-08, 6.23756851626e-08, 6.32957121259e-08, 6.42043903937e-08, 6.51025540077e-08,
    6.59909735447e-08, 6.68703634341e-08, 6.77413882848e-08, 6.8604668381e-08, 6.94607844804e-08, 7.03102820203e-08, 7.11536748229e-08, 7.1991448372e-08,
    7.2824062723e-08, 7.36519550992e-08, 7.44755422158e-08, 7.52952223703e-08, 7.61113773308e-08, 7.69243740467e-08, 7.77345662086e-08, 7.85422956743e-08,
    7.93478937793e-08, 8.01516825471e-08, 8.09539758128e-08, 8.17550802699e-08, 8.25552964535e-08, 8.33549196661e-08, 8.41542408569e-08, 8.49535474601e-08,
    8.57531242006e-08, 8.65532538723e-08, 8.73542180955e-08, 8.8156298059e-08, 8.89597752521e-08, 8.97649321908e-08, 9.05720531451e-08, 9.138142487e-08,
    9.21933373471e-08, 9.30080845407e-08, 9.38259651738e-08, 9.46472835298e-08, 9.54723502847e-08, 9.63014833769e-08, 9.71350089201e-08, 9.79732621669e-08,
    9.88165885297e-08, 9.96653446693e-08, 1.00519899658e-07, 1.0138063623e-07, 1.02247952126e-07, 1.03122261554e-07, 1.04003996769e-07, 1.04893609795e-07,
    1.05791574313e-07, 1.06698387725e-07, 1.07614573423e-07, 1.08540683296e-07, 1.09477300508e-07, 1.1042504257e-07, 1.11384564771e-07, 1.12356564007e-07,
    1.13341783071e-07, 1.14341015475e-07, 1.15355110887e-07, 1.16384981291e-07, 1.17431607977e-07, 1.18496049514e-07, 1.19579450872e-07, 1.20683053909e-07,
    1.21808209468e-07, 1.2295639141e-07, 1.24129212952e-07, 1.25328445797e-07, 1.26556042658e-07, 1.27814163916e-07, 1.29105209375e-07, 1.30431856341e-07,
    1.31797105598e-07, 1.3320433736e-07, 1.34657379914e-07, 1.36160594606e-07, 1.37718982103e-07, 1.39338316679e-07, 1.41025317971e-07, 1.42787873535e-07,
    1.44635331499e-07, 1.4657889173e-07, 1.48632138436e-07, 1.50811780719e-07, 1.53138707402e-07, 1.55639532047e-07, 1.58348931426e-07, 1.61313325908e-07,
    1.64596952856e-07, 1.68292495203e-07, 1.72541128694e-07, 1.77574279496e-07, 1.83813550477e-07, 1.92166040885e-07, 2.05295471952e-07, 2.22600839893e-07
};

inline const float PCG32::zig_fn[128] = {
    1.0f, 0.963598623011f, 0.936280813353f, 0.913041104253f, 0.892278506696f, 0.873239356919f, 0.855496407634f, 0.838778928349f,
    0.822902083699f, 0.807732738234f, 0.793171045519f, 0.779139726505f, 0.765577436082f, 0.752434456248f, 0.739669787677f, 0.727249120285f,
    0.715143377413f, 0.703327646455f, 0.691780377035f, 0.68048276891f, 0.669418297233f, 0.65857233912f, 0.647931876189f, 0.637485254896f,
    0.62722199145f, 0.617132611532f, 0.607208517467f, 0.597441877296f, 0.587825531465f, 0.578352913803f, 0.569017984198f, 0.559815170911f,
    0.550739320877f, 0.541785656682f, 0.532949739145f, 0.524227434628f, 0.515614886373f, 0.507108489253f, 0.498704867478f, 0.490400854812f,
    0.482193476986f, 0.47407993601f, 0.466057596125f, 0.458123971214f, 0.450276713467f, 0.442513603171f, 0.434832539473f, 0.427231532022f,
    0.419708693379f, 0.41226223212f, 0.404890446548f, 0.397591718955f, 0.390364510382f, 0.383207355816f, 0.376118859788f, 0.369097692334f,
    0.362142585282f, 0.355252328834f, 0.348425768415f, 0.341661801776f, 0.334959376311f, 0.328317486588f, 0.321735172063f, 0.31521151497f,
    0.308745638367f, 0.302336704338f, 0.29598391232f, 0.289686497571f, 0.283443729739f, 0.27725491156f, 0.271119377649f, 0.265036493387f,
    0.259005653912f, 0.253026283183f, 0.247097833139f, 0.241219782932f, 0.235391638239f, 0.229612930649f, 0.223883217122f, 0.218202079518f,
    0.212569124201f, 0.206983981709f, 0.201446306496f, 0.195955776745f, 0.190512094256f, 0.185114984406f, 0.179764196185f, 0.174459502324f,
    0.169200699492f, 0.1639876086f, 0.158820075195f, 0.153697969964f, 0.148621189348f, 0.143589656295f, 0.138603321143f, 0.133662162669f,
    0.128766189309f, 0.123915440582f, 0.119109988745f, 0.114349940703f, 0.10963544023f, 0.104966670533f, 0.100343857232f, 0.0957672718266f,
    0.0912372357329f, 0.0867541250127f, 0.082318375932f, 0.0779304915295f, 0.0735910494266f, 0.0693007111742f, 0.065060233529f, 0.0608704821745f,
    0.056732448584f, 0.05264727098f, 0.0486162607163f, 0.0446409359769f, 0.0407230655415f, 0.0368647267386f, 0.0330683839378f, 0.0293369977411f,
    0.0256741818288f, 0.0220844372634f, 0.0185735200577f, 0.0151490552854f, 0.0118216532614f, 0.00860719483079f, 0.00553245272614f, 0.00265435214565f
};

//#################################################### GENERIC NOISE CLASS ###########################################################

//Generic and virtual noise class: can be used to create new noises. compute_noise is used by the library to insert the noise in the SDE.
class NoiseClass{

public:
    
    //This function has to be override by the child. It is used in the CompFuncManager.
    //To override in YourNoise: -> return new YourNoise(*this);
    //To change seed see WienerEuler or WienerMilstein.
    virtual NoiseClass* clone() const = 0;

    //A standard compute noise that simply gives back the passed point in the system space.
    virtual void compute_noise(const valarray<float> &x_i,const float* h,valarray<float> &x_out);    

};

//############################################### PREDEFINED WIENER PROCESS WITH EULER-MARUYAMA NOISE ###########################################################

//A useful premade class to produce values of a Wiener process simulating dW with an Euler-Maruyama method.
class WienerEuler: public NoiseClass{

    PCG32 eng; //Random number generator engine.
    
public:

    //Cloning function. Used in the CompFuncManager.
    NoiseClass* clone() const override;

    //CONSTRUCTOR: initializes the distribution and random generator.
    //The seed of the random engine is created starting from the time and the point of the instance.
    WienerEuler();

    //Compute a dW kind of step for the Wiener process using the Euler-Maruyama method. 
    //Actually a vector of the same size of the system is returned and, in each slot, there is a different dW. 
    void compute_noise(const valarray<float> &x_i, const float* h,valarray<float> &x_out) override;

};

//############################################### PREDEFINED WIENER PROCESS WITH MILSTEIN NOISE ################################################

//A useful premade class to produce values of a Wiener process simulating dW with a Milstein method.
class WienerMilstein: public NoiseClass{

    PCG32 eng; //Random number generator engine.
    function<valarray<float>(valarray<float>)> D_g; //The derivative of the field g_function.

    valarray<float> dg_eval; //Preallocation of the valarray used by compute_noise

public:

    //Cloning function. Used in the CompFuncManager.
    NoiseClass* clone() const override;

    //CONSTRUCTOR: the distribution and the generator are initialize. 
    //Moreover, requires as argument the derivative of the g_function used in the field class.
    //The seed of the random engine is created starting from the time and the point of the instance.
    WienerMilstein(const function<valarray<float>(const valarray<float>&)>& derivative = nullptr);

    //Compute a dW kind-of step for the Wiener process using the Milstein method. 
    //Actually a vector of the same size of the system is returned and, at each step there is a different dW. 
    void compute_noise(const valarray<float> &x_i, const float* h, valarray<float> &x_out) override;

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
    friend class CompFuncManager;

public:

    //This function has to be override by the child. It is used in the CompFuncManager.
    //To override in YourField: -> return new YourField(*this);
    virtual FieldClass* clone() const = 0;

    //Implementation function for the deterministic part of the field. Override this as in
    //documentation to implement your system.
    virtual void f_function_impl(const valarray<float> &x,float t,valarray<float> &y) const;

    //Implementation function for the stochastic part of the field. Override this as in
    //documentation to implement your system.
    virtual void g_function_impl(const valarray<float> &x,float t,valarray<float> &y) const;

    //This function is optional and has not to be implemented. Can be override to add some
    //procedure that must not be done inside the RK4 steps but only at the beggining
    //(e.g. sorting a common random variable for the RK4 steps).
    virtual void start_step_procedure(const valarray<float> &x,float t){}

    //This function will give the result of the compute_noise of the local NoiseClass.
    void getNoise(const valarray<float> &x_i,const float* h, valarray<float> &x_out);

    //This function simply return the pointer to its associated NoiseClass object.
    NoiseClass* getNoiseClass(){
        return noise;
    }

};

//################################################# TRAJECTORY #############################################################################

//The Traj class is used to represent a simulated trajectory of a system of SDEs. It is purely an object class to manage easily the values of
//the variables and the time instants keeping them all in the same place. It is the product of the simulateTrajectory function.
class Traj{

    vector<float> times; //The vector of the time instant of the trajectory
    vector<float> vars; //The flat-2D array of the system variables. Rows: time; Columns: variables.
    unsigned short n_vars; //Number of variables for each time instant
    size_t step_num; //The time length/num of steps of the trajectory.

public:

    //CONSTRUCTOR 1: the time vector, the flat 2D array of variable values and the number of vars using std::move.
    Traj(vector<float>&& t,vector<float>&& v, unsigned short nv);

    //CONSTRUCTOR 2: passing another Traj instance.
    Traj(const Traj& t);

    //Given the time index and the variable return the position in the flat vector.
    const size_t flatNotation(const size_t& TI,const unsigned short& var) const{
        return (TI*n_vars+var);
    }

    //####################### GET FUNCTIONS ##############################

    //Return the reference to the vector of the time instants considered in the trajectory.
    const vector<float>& getTimes() const{
        return times;
    }

    //Return the reference to the 2D array of the variables in every time instant of the trajectory.
    const vector<float>& getVars() const{
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

    vector<float> points; //The values of the N simulations at the time instant T in a flat 2D matrix. (row,columns)=(trajs,variables).
    size_t N; //The number of simulations.
    unsigned int n_vars; //The number of variables.
    float T; //The time instant.

public:

    //CONSTRUCTOR 1: build using the flat 2D array of values, the time instant and the number of vars. The set of values must have the 
    //different trajectories along the rows and the variables along the columns.
    TimePicture(vector<float> values,float t,unsigned int nv);

    //CONSTURCTOR 2: build copying another TimePicture instance.
    TimePicture(const TimePicture& TP);

    //####################### GET FUNCTIONS ##############################

    //Return the reference to the 2D array of the values of the different trajectories at the time instant of the picture.   
    const vector<float>& getAllPoints() const{
        return points;
    }

    //Return the reference to the array of the values of a specific trajectory "p" at the time instant of the picture. 
    const vector<float> getPoint(size_t p) const{
        vector<float> output(n_vars,0.0f);

        std::memcpy(output.data(),points.data()+p*n_vars,n_vars);

        return output;
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
    vector<float> vars; //The flat 2D array of the system variables. Rows: time instants; Columns: variables.

    unsigned short n_vars; //Number of variables for each time instant

public:

    //CONSTRUCTOR 1: the time instants vector, the flat 2D array of variable values and the number of vars using std::move.
    SetOfPoints(vector<float>&& t,vector<float>&& v,unsigned short nv);

    //CONSTRUCTOR 2: passing another SetOfPoints instance.
    SetOfPoints(const SetOfPoints& sop);

    //Given the time index and the variable return the position in the flat vector.
    const size_t flatNotation(const size_t& TI,const unsigned short& var) const{
        return (TI*n_vars+var);
    }

    //####################### GET FUNCTIONS ##############################

    //Return the reference to the vector of the time instants considered in the set.
    const vector<float>& getTimes() const{
        return times;
    }

    //Return the reference to the 2D array of the variables in every time instant of the set.
    const vector<float>& getVars() const{
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

    //################### SUPPORT FUNCTIONS #######################
    //Helper functions to simplify the logic of higher-level routines.

    //This function will perform the checks of the parameters of simulateTrajectory. 
    void checkTrajInput(const vector<float> &x0,const float period,const float h_0);

    //################### EVOLUTION FUNCTIONS ####################
    //These functions are used to perform the evolutive aspects of the trajectory's computation.

    //Given the previous point and the step length, this internal function is the core function to evolve 
    //the last step in the new one of the trajectory. It will use the RK4 method.
    //The idea is to use a setup in the additive splitting way synergizing with RK4 and the noise method.
    void evolveTraj(const valarray<float> &x_n,float h,float t,vector<valarray<float>> &k,
        valarray<float> &x_next,valarray<float> &g,valarray<float> &n,valarray<float> &x_temp);

    //Given the starting point and the step, this function will return the RK4 update
    //of the deterministic part of the field.
    //REMEMBER: you have to multiply externally by h or your step eventually.
    void RK4_method(const valarray<float> &x0,float h,float t,vector<valarray<float>> &k,
        valarray<float> &update);

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

    //###################### PUBLIC UTILITY FUNCTIONS ##########################

    //This function will set the bound function. The bound function should return a boolean
    //and has as input a vector<float> (the point). The function should return "true" when
    //the point is in the domain.
    //The system has to be bounded (construction) to this function to work.
    void setBoundFunction(const function<bool(const valarray<float>&)>& f);

    //Given a vector of times and a time index this function will find the slot just before the given time instant.
    //This function is also STATIC.
    static size_t findTimeIndex(const vector<float> &times,const float TI);
    
};

//############################################# COMPLEX FUNCTIONS MANAGER #######################################################################

//The Complex Function Manager (CompFuncManager) is used to perform in an optimized way complex or heavy function. Most of them requires a lot
//of simulation. This is done via parallel computing. A reference system features has to be passed to create than copies for each thread.
//Aside from this, CompFuncManager is a merely box class with some useful yet heavy operations inside.
class CompFuncManager{

    const float PAR_THR_INDEX{0.1f}; //Ratio Nsim/NumThread threshold to opt parallelization of some functions.

    unsigned int ref_size; //Size of the reference system/problem.
    FieldClass* ref_field; //Field of the reference system.
    bool ref_bounded{false}; //Is the reference system bounded?
    function<bool(valarray<float>)> ref_bounds; //the function used to express the bounds of the domain of the reference system, if present.

    unsigned int NumThreads{8}; //The number of threads used for the parallelizable operations.

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ INTERNAL CHECK FUNCTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    //This function will perform the checks of the parameters of computeTimePicture.
    void checkFunctionComputeTimePicture(const unsigned int Nsim, const bool random_initial,const vector<vector<float>> &x0,
                                        const function<vector<float>()>& random_f);

    //This function will perform the checks of the parameters of PDF_1D.
    void checkFunctionPDF_1D(const unsigned int Nbins,const unsigned int axis,const bool adaptive,const vector<float> &domain);

    //This function will perform the checks of the parameters of PDF_2D.
    void checkFunctionPDF_2D(const vector<unsigned int> Nbins,const vector<unsigned int> axis,const bool adaptive,const vector<float> &domain);

    //This function will perform the checks of the parameters of computeAutocorrelation.
    void checkAutocorrelationInput(const Traj& traj,const unsigned int axis,const float tau);

public:

    //Constructor: requires the size of the reference system, the Field Class of the reference system
    //and a bool to say if the reference system is bounded or not. If bounder you need also to pass a
    //bound function valarray<float>->bool.
    CompFuncManager(unsigned int N,FieldClass* F,bool isBounded=false,const function<bool(const valarray<float>&)>& f = nullptr); 

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ TOOL FUNCTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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
    //The output of the function will be a TimePicture class object.     
    TimePicture produceTimePicture(const float period,const float h_0,unsigned int Nsim,
                                   const bool random_initial=false,const vector<vector<float>> &x0 = {{}},
                                   const function<vector<float>()>& random_f = nullptr,const float time_instant = -1.0f);

    //This function will produce starting from a TimePicture such the one produced by "produceTimePicture" a 1D bin
    //system useful to obtain PDFs. This function require:
    //- A TimePicture (MANDATORY).
    //- The number of bins (MANDATORY).
    //- The axis along which doing the bins (0: times,n: the n variable of the system) (MANDATORY).
    //- A bool to express if the binning domain is given or adaptive (true = adaptive).
    //- If false a vector with upper and lower domain is required.
    //The function will return a 2D vector with in each row the central value (first column) and the bin value (second column).
    //If there are too much threads than simulations (less than 10 simulation per thread) the code will be executed serially!!!
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
    //If there are too much threads than simulations (less than 10 simulation per thread) the code will be executed serially!!!
    vector<vector<float>> PDF_2D(const TimePicture& picture,vector<unsigned int> Nbins,vector<unsigned int> axis,
                                bool adaptive = false, vector<float> domain = {0.0,0.0,0.0,0.0});

    //This function will produce, starting from a trajectory such the one produced by "produceTimePicture" the
    //autocorrelation of the trajectory for a certain time delay. This function require:
    //- A trajectory in style vector<vector<float>> (MANDATORY).
    //- The axis (variable) along which computing the autocorrelation (MANDATORY).
    //- The time delay (tau) of the autocorrelation (MANDATORY).
    //The output will be the autocorrelation value computed as covariance/variance.
    //Also the time delay is converted in the nearest below number of steps.
    //If there are too much threads than simulations (less than 10 simulation per thread) the code will be executed serially!!!
    float computeAutocorrelation(const Traj& traj,unsigned int axis,float tau);

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ UTILITY FUNCTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    //This function is used to set the number of threads used by the heavy functions of the class.
    //The standard value is 8.
    void setNumThreads(unsigned int N);

    //Given a vector of times and a time index this function will find the slot just before the given time instant.
    //This function is also STATIC.
    static size_t findTimeIndex(const vector<float> &times,const float TI);

};

#endif