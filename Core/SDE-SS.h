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

public:

    //CONSTRUCTOR: The first parameter (seed) set the state while seq (the second parameter)
    //allow to have different sequences with the same seed.
    PCG32(uint64_t seed = 0x4d595df4d0f33173ULL, uint64_t seq = 1u) {
        state = 0U;
        inc = (seq << 1u) | 1u;
        next();
        state += seed;
        next();
    }

    //Obtain the next random 64-bit state.
    uint32_t next() {
        uint64_t oldstate = state;
        //Linear LCG: multiplication and increment
        state = oldstate * 6364136223846793005ULL + inc;
        //Output Transformation  (XSH RR)
        uint32_t xorshifted = static_cast<uint32_t>(((oldstate >> 18u) ^ oldstate) >> 27u);
        uint32_t rot = static_cast<uint32_t>(oldstate >> 59u);
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31u));
    }

    //Generate a random float in [0,1) with 24 bit precision (float's mantiss)
    float nextFloat() {
        return (next() >> 8) * 0x1.0p-24f;
    }

    //Generate a gaussian distributed number using the Ziggurat Algorithm.
    float nextGaussian();
};

//ZIGGURAT TABLES (128 BLOCKS):
//Internal bounds
static const uint32_t zig_kn[128] = {
    0, 3749721382, 3866291672, 3925762886, 3965905151, 3996227494, 4020584400, 4040941893,
    4058434698, 4073749504, 4087353986, 4099591461, 4110727915, 4120970347, 4130477141, 4139368817,
    4147735395, 4155649987, 4163171329, 4170348701, 4177224217, 4183834382, 4190209930, 4196377852,
    4202361958, 4208183577, 4213861895, 4219414212, 4224855909, 4230199321, 4235454924, 4240631627,
    4245736932, 4250777033, 4255757049, 4260681146, 4265552599, 4270373998, 4275147817, 4279876359,
    4284561803, 4289206214, 4293811559, 4298379634, 4302912117, 4307410537, 4311876307, 4316310723,
    4320715014, 4325090333, 4329437783, 4333758398, 4338053158, 4342322998, 4346568813, 4350791452,
    4354991733, 4359170442, 4363328334, 4367466133, 4371584534, 4375684200, 4379765759, 4383829813,
    4387876939, 4391907693, 4395922610, 4399922204, 4403906967, 4407877372, 4411833871, 4415776900,
    4419706877, 4423624207, 4427529278, 4431422461, 4435304113, 4439174577, 4443034188, 4446883269,
    4450722131, 4454551075, 4458370392, 4462180362, 4465981251, 4469773315, 4473556799, 4477331940,
    4481098967, 4484858097, 4488609539, 4492353495, 4496089163, 4499818731, 4503541400, 4507257358,
    4510966787, 4514669865, 4518366768, 4522057667, 4525742728, 4529422114, 4533095984, 4536764500,
    4540427820, 4544086100, 4547739493, 4551388152, 4555032228, 4558671869, 4562307222, 4565938431,
    4569565640, 4573188989, 4576808617, 4580424660, 4584037254, 4587646532, 4591252627, 4594855669,
    4598455787, 4602053106, 4605647754, 4609239855, 4612829532, 4616416905, 4619992019, 4623564998
};
//Scaling factors
static const float zig_wn[128] = {
    0.0008003102f, 0.0007559146f, 0.0007217031f, 0.0006941198f, 0.0006711467f, 0.0006515250f, 0.0006344212f, 0.0006192774f,
    0.0006057053f, 0.0005934177f, 0.0005821935f, 0.0005718642f, 0.0005623010f, 0.0005534024f, 0.0005450871f, 0.0005372895f,
    0.0005299534f, 0.0005230303f, 0.0005164778f, 0.0005102591f, 0.0005043425f, 0.0004987011f, 0.0004933111f, 0.0004881515f,
    0.0004832039f, 0.0004784518f, 0.0004738804f, 0.0004694770f, 0.0004652294f, 0.0004611270f, 0.0004571603f, 0.0004533206f,
    0.0004496001f, 0.0004459914f, 0.0004424881f, 0.0004400000f, 0.0004358137f, 0.0004326173f, 0.0004295055f, 0.0004264738f,
    0.0004235183f, 0.0004206354f, 0.0004178220f, 0.0004150752f, 0.0004123924f, 0.0004097712f, 0.0004072092f, 0.0004047043f,
    0.0004022543f, 0.0003998573f, 0.0003975113f, 0.0003952146f, 0.0003929654f, 0.0003907621f, 0.0003886034f, 0.0003864877f,
    0.0003844137f, 0.0003823800f, 0.0003803856f, 0.0003784291f, 0.0003765095f, 0.0003746257f, 0.0003727766f, 0.0003709611f,
    0.0003691783f, 0.0003674272f, 0.0003657071f, 0.0003640169f, 0.0003623561f, 0.0003607237f, 0.0003591189f, 0.0003575412f,
    0.0003559897f, 0.0003544637f, 0.0003529626f, 0.0003514857f, 0.0003500325f, 0.0003486022f, 0.0003471944f, 0.0003458085f,
    0.0003444439f, 0.0003431001f, 0.0003417765f, 0.0003404727f, 0.0003391881f, 0.0003379224f, 0.0003366751f, 0.0003354458f,
    0.0003342340f, 0.0003330393f, 0.0003318613f, 0.0003306997f, 0.0003295541f, 0.0003284242f, 0.0003273096f, 0.0003262100f,
    0.0003251250f, 0.0003240545f, 0.0003229980f, 0.0003219553f, 0.0003209262f, 0.0003199103f, 0.0003189075f, 0.0003179174f,
    0.0003169397f, 0.0003159744f, 0.0003150210f, 0.0003140795f, 0.0003131496f, 0.0003122310f, 0.0003113236f, 0.0003104270f,
    0.0003095412f, 0.0003086659f, 0.0003078008f, 0.0003069459f, 0.0003061009f, 0.0003052657f, 0.0003044400f, 0.0003036237f,
    0.0003028167f, 0.0003020188f, 0.0003012297f, 0.0003004495f, 0.0002996778f, 0.0002989146f, 0.0002981597f, 0.0002974129f
};
//Function values
static const float zig_fn[128] = {
    1.0000000f, 0.9918931f, 0.9839305f, 0.9761096f, 0.9684281f, 0.9608833f, 0.9534731f, 0.9461951f,
    0.9390471f, 0.9320271f, 0.9251329f, 0.9183625f, 0.9117140f, 0.9051855f, 0.8987748f, 0.8924803f,
    0.8863001f, 0.8802324f, 0.8742755f, 0.8684277f, 0.8626873f, 0.8570526f, 0.8515222f, 0.8460944f,
    0.8407677f, 0.8355407f, 0.8304120f, 0.8253802f, 0.8204440f, 0.8156020f, 0.8108530f, 0.8061957f,
    0.8016289f, 0.7971515f, 0.7927622f, 0.7884599f, 0.7842436f, 0.7801121f, 0.7760645f, 0.7721000f,
    0.7682173f, 0.7644158f, 0.7606944f, 0.7570524f, 0.7534888f, 0.7500028f, 0.7465936f, 0.7432599f,
    0.7400013f, 0.7368171f, 0.7337063f, 0.7306682f, 0.7277021f, 0.7248073f, 0.7219832f, 0.7192288f,
    0.7165437f, 0.7139270f, 0.7113781f, 0.7088962f, 0.7064808f, 0.7041312f, 0.7018468f, 0.6996270f,
    0.6974712f, 0.6953787f, 0.6933491f, 0.6913817f, 0.6894759f, 0.6876313f, 0.6858473f, 0.6841235f,
    0.6824593f, 0.6808543f, 0.6793080f, 0.6778199f, 0.6763896f, 0.6750168f, 0.6737009f, 0.6724416f,
    0.6712386f, 0.6700914f, 0.6689997f, 0.6679632f, 0.6669814f, 0.6660541f, 0.6651811f, 0.6643621f,
    0.6635967f, 0.6628848f, 0.6622261f, 0.6616205f, 0.6610677f, 0.6605676f, 0.6601199f, 0.6597246f,
    0.6593816f, 0.6590906f, 0.6588517f, 0.6586649f, 0.6585300f, 0.6584472f, 0.6584164f, 0.6584377f,
    0.6585112f, 0.6586369f, 0.6588149f, 0.6590453f, 0.6593283f, 0.6596640f, 0.6600525f, 0.6604941f,
    0.6609889f, 0.6615372f, 0.6621392f, 0.6627953f, 0.6635058f, 0.6642709f, 0.6650912f, 0.6659668f,
    0.6668984f, 0.6678864f, 0.6689312f, 0.6700334f, 0.6711936f, 0.6724121f, 0.6736897f, 0.6750269f
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
    CompFuncManager(unsigned int N,FieldClass* F,bool isBounded=false,const function<bool(const valarray<float>&)>& f); 

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