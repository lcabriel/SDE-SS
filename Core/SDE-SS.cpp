#include "SDE-SS.h"

using namespace std;

//###################################################################################################################################################

//OVERRIDING SUM OPERATOR BETWEEN VECTOR<FLOAT>
vector<float> operator+(const vector<float>& a, const vector<float>& b) {
    if (a.size() != b.size()) {
        throw invalid_argument("Vectors should have the same dimension");
    }

    vector<float> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }

    return result;
}

//OVERRIDING MULT WITH FLOAT OPERATION FOR VECTOR<FLOAT>
vector<float> operator*(const vector<float>& vec, float scalar) {
    vector<float> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * scalar;
    }
    return result;
}

// scalar * vector
vector<float> operator*(float scalar, const vector<float>& vec) {
    return vec * scalar;  //We use the overriding above
}

//OVERRIDING MULTI BETWEEN VECTOR<FLOAT>: IT IS THE MULT BETWEEN HOMOLOGOUS COMPONENTS (NOT SCALAR PRODUCT)
vector<float> operator*(const vector<float>& a,const vector<float>& b) {
    if (a.size() != b.size()) {
        throw invalid_argument("Vectors should have the same dimension");
    }

    vector<float> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i]*b[i];
    }

    return result;
}

//####################################################################################################################################################

//Generate a gaussian distributed number using the Ziggurat Algorithm.
float PCG32::nextGaussian(){
    for (;;) {
        uint32_t u = next();
        int32_t i = u & 127; // Choose one of the blocks
        
        // If the number falls in the internal part of the block, return it.
        if (u < zig_kn[i]) return static_cast<float>(u) * zig_wn[i];

        // Otherwise, manage the blocks border or the tail situation.
        if (i == 0) {
            float x, y;
            do {
                x = -logf(nextFloat()) * 0.2904764f;
                y = -logf(nextFloat());
            } while (y + y < x * x);
            return (u & 128) ? (3.442620f + x) : -(3.442620f + x);
        }

        //Refinement for the block borders
        float x = static_cast<float>(u) * zig_wn[i];
        if (zig_fn[i] + nextFloat() * (zig_fn[i - 1] - zig_fn[i]) < expf(-0.5f * x * x)) {
            return x;
        }
    }    
}

//##############################################################################################################################################################
//############################################################# GENERIC NOISE CLASS ########################################################################################
//##############################################################################################################################################################

//A standard compute noise that simply gives back the passed point in the system space.
void NoiseClass::compute_noise(const valarray<float>& x_i,const float* h,valarray<float>& x_out){}

//##############################################################################################################################################################
//############################################################# PREDEFINED WIENER PROCESS ################################################################################
//##############################################################################################################################################################

//CONSTRUCTOR: the distribution and the generator are initialized.
//The seed of the random engine is created starting from the time and the point of the instance.
WienerEuler::WienerEuler():
    eng(static_cast<unsigned long>(std::time(0) * reinterpret_cast<uintptr_t>(this)))
{}

//Compute a dW kind of step for the Wiener process using the Euler-Maruyama method. 
//Actually a vector of the same size of the system is returned and, in each slot, there is a different dW. 
void WienerEuler::compute_noise(const valarray<float>& x_i,const float* h,valarray<float> &x_out){
    for(size_t i=0;i<x_i.size();i++){ //dW -> N(0,1)*sqrt(h)
        x_out[i] = eng.nextGaussian()*sqrt(*h);
    }
}

NoiseClass* WienerEuler::clone() const{
    WienerEuler* n_noise = new WienerEuler(*this);

    auto new_seed{static_cast<unsigned int>(std::time(nullptr)) ^ static_cast<unsigned int>(reinterpret_cast<uintptr_t>(n_noise))};
    n_noise->eng = PCG32(new_seed);

    return n_noise;
}

//##############################################################################################################################################################
//############################################################# PREDEFINED WIENER PROCESS ################################################################################
//##############################################################################################################################################################

//CONSTRUCTOR: initializes the distribution and random generator.
//Moreover, requires as argument the derivative of the g_function used in the field class.
//The seed of the random engine is created starting from the time and the point of the instance.
WienerMilstein::WienerMilstein(const function<valarray<float>(const valarray<float>&)>& derivative):
    eng(static_cast<unsigned long>(time(0) * reinterpret_cast<uintptr_t>(this)))
{
    //Check if the argument is valid.
    if(!derivative) error("WienerMilstein: Constructor: Runtime error: Derivative of g_function argument is not valid.");

    D_g = derivative;
}

//Compute a dW kind-of step for the Wiener process using the Milstein method. 
//Actually a vector of the same size of the system is returned and, at each step there is a different dW. 
void WienerMilstein::compute_noise(const valarray<float> &x_i,const float* h, valarray<float> &x_out){
    dg_eval = D_g(x_i);
    for(size_t i=0;i<x_i.size();i++){
        const float dW = eng.nextGaussian()*sqrt(*h);

        x_out[i] = dW + (0.5f)*dg_eval[i]*(dW*dW-*h);
    }
}

NoiseClass* WienerMilstein::clone() const{
    WienerMilstein* n_noise = new WienerMilstein(*this);

    auto new_seed{static_cast<unsigned int>(std::time(nullptr)) ^ static_cast<unsigned int>(reinterpret_cast<uintptr_t>(n_noise))};
    n_noise->eng = PCG32(new_seed);

    return n_noise;    
}

//##############################################################################################################################################################
//########################################################### GENERIC FIELD CLASS ############################################################################
//##############################################################################################################################################################

//This public function allow to set the noise of the FieldClass in a controlled way.
void FieldClass::setNoise(NoiseClass* const N){

    if(N==nullptr) error("FieldClass: setNoise: Runtime error: Invalid NoiseClass passed to the function.");

    noise = N;
    noise_initialized = true;
}

//Function for the deterministic part of the field. To define, please
//override "f_function_impl" as in documentation.
void FieldClass::f_function(const valarray<float> &x,float t,valarray<float> &y){

    if(!noise_initialized) error("FieldClass: f_function: Runtime error: Noise in the FieldClass not initialized.");

    f_function_impl(x,t,y); 
}

//Function for the stochastic part of the field. To define, please
//override "g_function_impl" as in documentation.
void FieldClass::g_function(const valarray<float> &x,float t,valarray<float> &y){

    if(!noise_initialized) error("FieldClass: g_function: Runtime error: Noise in a FieldClass not initialized.");

    g_function_impl(x,t,y);
}

//Implementation function for the deterministic part of the field. Override this as in
//documentation to implement your system.
void FieldClass::f_function_impl(const valarray<float> &x,float t,valarray<float> &y) const{}

//Implementation function for the stochastic part of the field. Override this as in
//documentation to implement your system.
void FieldClass::g_function_impl(const valarray<float> &x,float t,valarray<float> &y) const{}

//This function will give the result of the compute_noise of the local NoiseClass.
void FieldClass::getNoise(const valarray<float> &x_i,const float* h,valarray<float> &x_out){
    
    if(!noise_initialized) error("FieldClass: getNoise: Runtime error: Noise in a FieldClass not initialized");

    return noise->compute_noise(x_i,h,x_out);
}

//##############################################################################################################################################################
//########################################################### DATALINKER #############################################################################
//##############################################################################################################################################################

//CONSTRUCTOR 1: passing a traj and specifying if is adaptive
DataLinker::DataLinker(const Traj& t,const bool AS): data(t), adapting_set(AS){}

//CONSTRUCTOR 2: copying a DataLinker
DataLinker::DataLinker(const DataLinker& DL): data(DL.data), adapting_set(DL.adapting_set){}

//Given a certain time instant, the variable will find the a time-index near the value (the next one). The TimeOptimizationCounter
//is used to improve performances allowing the bypass of this function when multiple getData access to the sime time index (see examples).
size_t DataLinker::findTimeIndex(const float t){
    //1. Try to look if is the same time
    if((t>=data.getTimes()[TOC])&&(t<data.getTimes()[TOC+1])){
        return TOC;
    }
    //2. Try to find if it is after the counter
    for(size_t i=TOC+1;i<data.getLength();i++){
        if((t>=data.getTimes()[i])&&(t<data.getTimes()[i+1])){
            TOC = i;
            return i;
        }
    }
    //3. Otherwise is before
    for(size_t i=0;i<TOC;i++){
        if((t>=data.getTimes()[i])&&(t<data.getTimes()[i+1])){
            TOC = i;
            return i;
        }        
    }

    cout << "DataLinker: findTimeIndex: Runtime error: time instant outside of the domain." << endl;
    exit(401);
}

//This function is made to be used inside of getData. It allows, instead of using the time index found with findTimeIndex, to obtain a sort
//in-between value in the dataset using a linear interpolator between the previous and the next point along the given variable.
//WARNING: this function could increase significatly your computation times and it is made for dataset characterised by a wide grid.
float DataLinker::interpolData(const float& t,const unsigned int& var){
    size_t t_i = findTimeIndex(t);

    const vector<float> p1 = data.getInstant(t_i);
    const vector<float> p2 = data.getInstant(t_i+1);

    return (p1[var]+(t-p1[0])*(p2[var]-p1[var])/(p2[0]-p1[0]));
}


//Given a certain time instant the getData returns a child-defined float value near to the given time.
float DataLinker::getData(const float t,const valarray<float>& x){    
    return 0.0;
}

//This function allow, if the DataLinker is adaptive, to pass a new set of data to the DataLinker.
void DataLinker::setNewData(const Traj& t){
    if(adapting_set) data=t;
    else{
        cout << "DataLinker: setNewData: Runtime error: trying to change the data set of a non-adaptive DataLinker instance." << endl;
        exit(402);
    }
}

//This function is used to reset the TimeOptimizationCounter which is used to optimize the finding process of the time instant during the
//simulations. 
void DataLinker::resetTimeOptimizationCounter(){
    TOC=0;
}

//##############################################################################################################################################################
//########################################################### TRAJECTORY #############################################################################
//##############################################################################################################################################################

//CONSTRUCTOR 1: the time instants vector, the flat 2D array of variable values and the number of vars using std::move.
Traj::Traj(vector<float>&& t, vector<float>&& v, unsigned short nv): times(move(t)), vars(move(v)), n_vars(nv) {
    if(times.size() != vars.size()/nv)
        error("Traj: Constructor: Runtime error: Times vector has different number of steps than the variables 2D vector");
    step_num = times.size();
}

//CONSTRUCTOR 2: passing another Traj instance.
Traj::Traj(const Traj& t): times(t.getTimes()),vars(t.getVars()),step_num(t.getLength()){}

//Given a certain time index, this function will return the situation as a vector of times+vars 
vector<float> Traj::getInstant(const size_t index) const{
    vector<float> output(n_vars+1,0.0);

    output[0] = times[index];

    for(unsigned int i=0;i<n_vars;i++){
        output[i+1] = vars[flatNotation(index,i)];
    }

    return output;
}

//##############################################################################################################################################################
//########################################################### TIME PICTURE #############################################################################
//##############################################################################################################################################################

//CONSTRUCTOR 1: build using the set of values and the time instant. The set of values must have the different trajectories
//along the rows and the variables along the columns.
TimePicture::TimePicture(vector<float> values,float t,unsigned int nv): points(values),T(t),n_vars(nv){
    N = values.size();
}

//CONSTURCTOR 2: build copying another TimePicture instance.
TimePicture::TimePicture(const TimePicture& TP): points(TP.getAllPoints()),T(TP.getTimeInstant()),N(TP.getNumSim()){}

//##############################################################################################################################################################
//########################################################### SET OF POINTS #############################################################################
//##############################################################################################################################################################

//CONSTRUCTOR 1: the time instants vector, the flat 2D array of variable values and the number of vars using std::move.
SetOfPoints::SetOfPoints(vector<float>&& t,vector<float>&& v,unsigned short nv): times(move(t)), vars(move(v)), n_vars(nv) {
    if(times.size() != vars.size())
        error("SetOfPoints: Constructor: Runtime error: Times vector has different number of instants than the variables 2D vector");
}

//CONSTRUCTOR 2: passing another SetOfPoints instance.
SetOfPoints::SetOfPoints(const SetOfPoints& sop): times(sop.getTimes()),vars(sop.getVars()){}

//Given a certain index, this function will return the situation as a vector of times+vars 
vector<float> SetOfPoints::getInstant(const size_t index) const{
    vector<float> output(n_vars+1,0.0);

    output[0] = times[index];

    for(unsigned int i=0;i<n_vars;i++){
        output[i+1] = vars[flatNotation(index,i)];
    }

    return output;
}

//##############################################################################################################################################################
//######################################################### SYSTEM #############################################################################################
//##############################################################################################################################################################

//Constructor: requires the size of the system, the Field Class of the system
//and a bool to say if the system is bounded or not (rem: you need to set the "BoundFunction"
//with the specificy method of this class, in this case).
SDE_SS_System::SDE_SS_System(unsigned int N,FieldClass* F,bool isBounded):
    field(F), bounded(isBounded), size(N)
{
    //Check if the size of the problem is meaningful
    if(F == nullptr) error("SDE_SS_System: Constructor: Runtime error: Passed FieldClass is invalid (nullptr).");
    if(N<1) error("SDE_SS_System: Constructor: Runtime error: Given size is <1. It has no sense.");

}

//################### EVOLUTION FUNCTIONS ####################
//These functions are used to perform the evolutive aspects of the trajectory's computation.

//Given the previous point and the step length, this internal function is the core function to evolve 
//the last step in the new one of the trajectory. It will use the RK4 method.
//The idea is to use a setup in the Strang splitting way synergizing with RK4 and the noise method.
void SDE_SS_System::evolveTraj(const valarray<float> &x_n,float h,float t,vector<valarray<float>> &k,
    valarray<float> &x_next,valarray<float> &g,valarray<float> &n,valarray<float> &x_temp){

    //1. We simulate the deterministic part of the field with RK4 for half step
    //Now we compute X^1 (the first intermediate step)
    RK4_method(x_n,h/2,t,k,x_next);
    x_next *= (h/2.0f);
    x_next += x_n;

    //2. Add the stochastic part
    field->g_function(x_next,t,g);
    field->getNoise(x_next,&h,n);

    g *= n;
    x_next += g; //Compute the noisy part; is multiplied than for g_function
    
    //3. Evolve the last half step with RK4 for half step
    RK4_method(x_next,h/2,t,k,x_temp);
    x_temp *= (h/2.0f);
    x_next += x_temp;
}

//Given the starting point and the step, this function will return the RK4 update
//of the deterministic part of the field.
//REMEMBER: you have to multiply externally by h or your step eventually.
void SDE_SS_System::RK4_method(const valarray<float> &x0,float h,float t,vector<valarray<float>> &k,
    valarray<float> &update){

    update=x0;

    field->f_function(update,t,k[0]); //compute k1
    update = k[0];
    update *= (h/2.0f);
    update += x0;

    field->f_function(update,t+h/2.0f,k[1]); //compute k2
    update = k[1];
    update *= (h/2.0f);
    update += x0;

    field->f_function(update,t+h/2.0f,k[2]); //compute k3
    update = k[2];
    update *= h;
    update += x0;

    field->f_function(update,t+h,k[3]); //compute k4

    //Note: the final for is more optimize for this operation instead of using the operators:
    const float sixth{1.0f/6.0f};
    for(unsigned int i=0;i<size;i++){
        update[i] = (k[0][i]+2.0f*k[1][i]+2.0f*k[2][i]+k[3][i])*sixth;
    }
}

//##################### CORE FUNCTIONS ##########################################

//This core function simulate a single trajectory given the initial conditions, the time period and the standard step.
//This function will return a 2D vector with all the points of the trajectory. 
//In the returning vector the times and the values are conjointed. The first column is used for times.
Traj SDE_SS_System::simulateTrajectory(const vector<float> &x0,const float period,const float h_0){
    checkTrajInput(x0,period,h_0); //Check if the inputs are meaningful

    //Setup the initial point of the trajectory
    vector<float> values;
    vector<float> times;

    //Estimation of the size (for memory opt)
    values.reserve((period*1.1f/h_0)*size);
    times.reserve(period*1.1f/h_0);

    values.insert(values.end(),x0.begin(),x0.end());
    times.push_back(0.0f);

    float h{h_0}; //The h used step-by-step. Can change to adapt to the boundary

    //k vector used by the RK4 method to avoid reallocation
    vector<valarray<float>> ks(4,valarray<float>(0.0f,size));

    valarray<float> x1,g,n,x_temp; //Used for pooling in evolveTraj and RK4_method

    g = valarray<float>(0.0f,size);
    n = valarray<float>(0.0f,size);
    x_temp = valarray<float>(0.0f,size);

    valarray<float> current_point(0.0f,size);
    valarray<float> new_point(0.0f,size);

    do{
        std::copy(values.end()-size,values.end(),begin(current_point));

        evolveTraj(current_point,h,times.back(),ks,new_point,g,n,x_temp); //define the new point

        //if out of bound
        if(bounded){
            while(!bounds(new_point)){
                h = h/10.0f;
                evolveTraj(current_point,h,times.back(),ks,new_point,g,n,x_temp);
            }
        }

        //Add the new point
        values.insert(values.end(),begin(new_point),end(new_point));
        times.emplace_back(times.back()+h);

        //Reset the step
        h = h_0;

    }while((times.back()+h)<=period);

    //Return the Traj object
    return Traj(move(times),move(values),size);
}

//This core function acts as simulateTrajectory, however the internal points of the trajectory are not saved and only
//the last value of the trajectory is returned as a vector of shape (time,[coords]).
vector<float> SDE_SS_System::simulateTrajectoryLastPoint(const vector<float> &x0,const float period,const float h_0){
    checkTrajInput(x0,period,h_0);

    //Setup the initial point of the trajectory
    valarray<float> old_point = valarray<float>(x0.data(),size);
    float time{0.0f};

    float h{h_0}; //The h used step-by-step. Can change to adapt to the boundary

    vector<valarray<float>> ks(4,valarray<float>(0.0f,size)); //k vector used by the RK4 method to avoid reallocation

    valarray<float> x1,g,n,x_temp; //Used for pooling in evolveTraj and RK4_method

    g = valarray<float>(0.0f,size);
    n = valarray<float>(0.0f,size);
    x_temp = valarray<float>(0.0f,size);

    valarray<float> new_point(0.0f,size);

    do{

        evolveTraj(old_point,h,time,ks,new_point,g,n,x_temp); //define the new point

        //if out of bound
        if(bounded){
            while(!bounds(new_point)){
                h = h/10.0f;
                evolveTraj(old_point,h,time,ks,new_point,g,n,x_temp);
            }
        }

        //Substitute the old point
        old_point.swap(new_point);
        time += h;

        //Reset the step
        h = h_0;

    }while((time+h)<=period);

    vector<float> point(size + 1);
    point[0] = time;
    for(unsigned int i = 0; i < size; ++i) {
        point[i+1] = old_point[i];
    }
    return point;
}

//This core function acts as simulateTrajectory but it will not return the entire trajectory. In fact, this function
//asks for a set of time instants and the output will be the values of the trajectory in those instants (actually the
//immediately before point).
SetOfPoints SDE_SS_System::simulateTrajectorySOP(const vector<float> &x0,const float period,const float h_0,const vector<float> &instants){
    checkTrajInput(x0,period,h_0);

    //Setup the initial point of the trajectory
    valarray<float> old_point = valarray<float>(x0.data(),size);
    float time{0.0f};

    vector<float> values;
    vector<float> times;

    values.reserve(instants.size()*size);
    times.reserve(instants.size());

    //Setup all the other useful variables

    vector<float> sort_inst = instants;
    sort(sort_inst.begin(),sort_inst.end()); //Sort the instants to improve performances

    float h{h_0}; //The h used step-by-step. Can change to adapt to the boundary
    unsigned int counter{0}; //Used to keep track of the last found point (to avoid striving at every iteration)

    vector<valarray<float>> ks(4,valarray<float>(0.0f,size)); //k vector used by the RK4 method to avoid reallocation
    valarray<float> g,n,x_temp; //Used for pooling in evolveTraj and RK4_method

    g = valarray<float>(0.0f,size);
    n = valarray<float>(0.0f,size);
    x_temp = valarray<float>(0.0f,size);

    valarray<float> new_point(0.0f,size);

    do{

        evolveTraj(old_point,h,time,ks,new_point,g,n,x_temp); //define the new point

        //if out of bound
        if(bounded){
            while(!bounds(new_point)){
                h = h/10.0f;
                evolveTraj(old_point,h,time,ks,new_point,g,n,x_temp);
            }
        }

        //Check if it is in the set the previous point
        if(counter < instants.size()){
            if((time+h) >= sort_inst[counter]){
                counter++;
                values.insert(values.end(),begin(old_point),end(old_point));
                times.push_back(time);
            }
        }

        //Update time
        time += h;

        //Substitute the old point
        old_point.swap(new_point);

        //Reset the step
        h = h_0;

    }while((time+h)<=period);

    return SetOfPoints(move(times),move(values),size);
}

//################### PUBLIC UTILITY FUNCTIONS ####################################

//This function will set the bound function. The bound function should return a boolean
//and has as input a vector<float> (the point). The function should return "true" when
//the point is in the domain.
//The system has to be bounded (construction) to this function to work.
void SDE_SS_System::setBoundFunction(const function<bool(const valarray<float>&)>& f = nullptr){
    //Check that the system is defined as bounded and the function defined
    if(!bounded) error("SDE_SS_System: setBoundFunction: Runtime error: SDE_SS_System is not bounded but setBoundFunction is called.");
    if(!f) error("SDE_SS_System: setBoundFunction: Runtime error: setBoundFunction argument is not valid.");

    bounds = f;
}

//##################################################### SUPPORT FUNCTIONS #########################################################################
//Helper functions to simplify the logic of higher-level routines.

//This function will perform the checks of the parameters of simulateTrajectory. 
void SDE_SS_System::checkTrajInput(const vector<float> &x0,const float period,const float h_0){
    
    //Check correct size of the initial condition vector, the positivity of the period and the meaningfulness of the step
    if(x0.size()!=size) error("SDE_SS_System: simulateTrajectory: Runtime error: Initial condition vector have more or less dimensions than the system.");
    if(period<=0) error("SDE_SS_System: simulateTrajectory: Runtime error: Negative or null time period for the trajectory.");
    if((h_0<=0)||(h_0>period)) error("SDE_SS_System: simulateTrajectory: Runtime error: Unmeaningful time step for the trajectory.");
}

//Given a vector of times and a time index this function will find the slot just before the given time instant.
//This function is also STATIC.
size_t SDE_SS_System::findTimeIndex(const vector<float> &times,const float TI){
    for(size_t i=1;i<times.size();i++){
        if(TI<times[i]){
            return i-1;
        }
    }

    cout << "SDE_SS_System: findTimeIndex: Runtime error: given time instant is not present!" << endl;
    exit(205);
}

//##############################################################################################################################################################
//######################################################### COMPLEX FUNCTION MANAGER ###########################################################################
//##############################################################################################################################################################

//Constructor: requires the size of the reference system, the Field Class of the reference system
//and a bool to say if the reference system is bounded or not. If bounder you need also to pass a
//bound function valarray<float>->bool.
CompFuncManager::CompFuncManager(unsigned int N,FieldClass* F,bool isBounded=false,const function<bool(const valarray<float>&)>& f = nullptr):
    ref_field(F),ref_bounded(isBounded),ref_size(N){
        //Check if the size of the problem is meaningful
        if(F == nullptr) error("CompFuncManager: Constructor: Runtime error: Passed FieldClass is invalid (nullptr).");
        if(N<1) error("CompFuncManager: Constructor: Runtime error: Given size is <1. It has no sense.");

        //If bounded check that a bound function is passed:
        if(isBounded){
            if(f==nullptr) error("CompFuncManager: Constructor: Runtime error: Passed bound function missing or invalid.");
        }
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ PUBLIC UTILITY FUNCTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//This function is used to set the number of threads used by the complex functions.
//The standard value is 8.
void CompFuncManager::setNumThreads(unsigned int N){
    NumThreads = N;
}

//Given a vector of times and a time index this function will find the slot just before the given time instant.
//This function is also STATIC.
size_t CompFuncManager::findTimeIndex(const vector<float> &times,const float TI){
    for(size_t i=1;i<times.size();i++){
        if(TI<times[i]){
            return i-1;
        }
    }

    cout << "CompFuncManager: findTimeIndex: Runtime error: given time instant is not present!" << endl;
    exit(205);
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ INTERNAL CHECK FUNCTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//This function will perform the checks of the parameters of computeTimePicture.
void CompFuncManager::checkFunctionComputeTimePicture(const unsigned int Nsim,const bool random_initial,const vector<vector<float>> &x0,
                                                    const function<vector<float>()>& random_f){
    //1.0: we do not check period and h_0 because it is done in simulateTrajectory

    //1.1: we briefly check that Nsim is not 0, the function is defined and the size is meaningful.
    if(Nsim==0) error("CompFuncManager: computeTimePicture: Runtime error: produceTimePicture cannot accept Nsim=0.");
    if(random_initial && !random_f) error("CompFuncManager: computeTimePicture: Runtime error: random initial condition function not defined.");
    if(!random_initial && (x0.size()!=Nsim || x0[0].size()!=ref_size)) error("CompFuncManager: computeTimePicture: Runtime error: initial condition has an invalid size.");
}

//This function will perform the checks of the parameters of PDF_1D.
void CompFuncManager::checkFunctionPDF_1D(const unsigned int Nbins,const unsigned int axis,const bool adaptive,const vector<float> &domain){
    //Check if the number of bins meaningful, the axis is not bigger than the dimensionality and the domain is valid.
    if(Nbins<1) error("CompFuncManager: PDF_1D: Runtime error: Number of bins passed is less than 1.");
    if(axis>=(ref_size+1)) error("CompFuncManager: PDF_1D: Runtime error: No variable associated with that axis of the ref system.");
    if(domain[0]==domain[1]) error("CompFuncManager: PDF_1D: Runtime error: Undefined or domain with equal extremes.");
}

//This function will perform the checks of the parameters of PDF_2D.
void CompFuncManager::checkFunctionPDF_2D(const vector<unsigned int> Nbins,const vector<unsigned int> axis,
                                        const bool adaptive,const vector<float> &domain){
    //1.1 Check if the number of bins meaningful
    //the axis is not bigger than the dimensionality and the domain is valid. This is done for both the axis.
    if(Nbins[0]<1 || Nbins[1]<1) error("CompFuncManager: PDF_2D: Runtime error: Number of bins passed is less than 1.");
    if(axis[0]>=(ref_size+1) || axis[1]>=(ref_size+1)) error("CompFuncManager: PDF_2D: Runtime error: No variable associated with that axis of the reference system.");
    if(domain[0]==domain[1]||domain[2]==domain[3]) error("CompFuncManager: PDF_2D: Runtime error: Undefined or domain with equal extremes.");
}

//This function will perform the checks of the parameters of computeAutocorrelation.
void CompFuncManager::checkAutocorrelationInput(const Traj& traj,const unsigned int axis,const float tau){

    //Check that the axis actually exist in comparison with the system size. Then check if the time window is lesser than the trajectory length.
    //At last, verify also that the time window is meaningful thus strictly positive.
    if(axis>=ref_size) error("CompFuncManager: computeAutocorrelation: Runtime error: No variable associated with the given autocorrelation axis.");
    if(tau>(traj.getTimes()).back()) error("CompFuncManager: computeAutocorrelation: Runtime error: Time delay bigger than the trajectory length.");
    if(tau<=0) error("CompFuncManager: computeAutocorrelation: Runtime error: Time delay should be bigger than 0.");
}

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ TOOL FUNCTIONS $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

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
TimePicture CompFuncManager::produceTimePicture(const float period, const float h_0, unsigned int Nsim, const bool random_initial,
                                                const vector<vector<float>> &x0,const function<vector<float>()> &random_f,
                                                const float time_instant)
{
    // 1. Check if the input are meaningful.
    checkFunctionComputeTimePicture(Nsim, random_initial, x0, random_f);

    // 2. Prepare the flat 2D output vector.
    vector<float> points(Nsim*(ref_size + 1),0.0f);

    float mean_time_instant{0.0f};

    // 3. Perform the parallel computations
    #pragma omp parallel num_threads(NumThreads) reduction(+:mean_time_instant)
    {
        const unsigned short ID{omp_get_thread_num()};
        const unsigned short LNT{NumThreads};

        //Compute the size of the chunk for the thread.
        const unsigned int trajectories_per_thread{Nsim/LNT};
        const unsigned int remainder{Nsim%LNT};
        const unsigned int step{ref_size+1};

        const unsigned int my_start_traj{ID*trajectories_per_thread+min((unsigned int)ID,remainder)};
        const unsigned int my_end_traj{my_start_traj+trajectories_per_thread+(ID<remainder ? 1 : 0)};

        const unsigned int chunk_begin{my_start_traj*step};
        const unsigned int chunk_end{my_end_traj*step};
        
        //Create the thread's own system.
        FieldClass* th_field{ref_field->clone()};
        NoiseClass* th_noise{ref_field->getNoiseClass()->clone()};
        th_field->setNoise(th_noise);
        SDE_SS_System system(ref_size,th_field,ref_bounded);

        if(ref_bounded){
            system.setBoundFunction(ref_bounds);
        }

        //Some preallocation
        vector<float> init(ref_size,0.0f);
        vector<float> local_p(ref_size+1,0.0f);

        for(unsigned int i=chunk_begin;i<chunk_end;i+=step){
            if(random_initial){
                init = random_f();
            }
            else{
                init = x0[i/step];
            }
            
            if(time_instant<=0.0f){
                local_p = system.simulateTrajectoryLastPoint(init,period,h_0);
                mean_time_instant += local_p[0];
            }
            else{
                local_p = (system.simulateTrajectorySOP(init,period,h_0,{time_instant})).getInstant(0);
            }

            std::memcpy(points.data()+i,local_p.data(),step*sizeof(float));
        }

        delete th_noise;
        delete th_field;
    }    

    // 4. Final time instant
    if (time_instant <= 0.0f)
        mean_time_instant /= static_cast<float>(Nsim);
    else
        mean_time_instant = time_instant;

    return TimePicture(std::move(points), mean_time_instant, ref_size);
}

//This function will produce starting from a TimePicture such the one produced by "produceTimePicture" a 1D bin
//system useful to obtain PDFs. This function require:
//- A TimePicture (MANDATORY).
//- The number of bins (MANDATORY).
//- The axis along which doing the bins (0: times,n: the n variable of the system) (MANDATORY).
//- A bool to express if the binning domain is given or adaptive (true = adaptive).
//- If false a vector with upper and lower domain is required.
//The function will return a 2D vector with in each row the central value (first column) and the bin value (second column).
//If there are too much threads than simulations (less than 10 simulation per thread) the code will be executed serially!!!
vector<vector<float>> CompFuncManager::PDF_1D(const TimePicture& picture,unsigned int Nbins,unsigned int axis,
                                            bool adaptive,vector<float> domain)
{
    //0. Perform some checks on the arguments:
    checkFunctionPDF_1D(Nbins,axis,adaptive,domain);

    //1. Import data and features
    const vector<float>& data{picture.getAllPoints()};
    const unsigned int Nsim{picture.getNumSim()};
    const unsigned int step{ref_size+1};

    //2. Compute the limits of the domain along the given axis
    float v_max,v_min;

    if(adaptive){
        v_max = -1e30f;
        v_min = 1e30f;

        #pragma omp parallel for num_threads(NumThreads) reduction(min:v_min) reduction(max:v_max) if(Nsim/NumThreads>=PAR_THR_INDEX)
        for(size_t i=0;i<data.size();i+=step){
            float v{data[i+axis]};
            if(v<v_min) v_min = v;
            if(v>v_max) v_max = v;
        }
    }
    else{
        v_min = domain[0];
        v_max = domain[1];
    }

    //3. Setup the bins
    const float delta{(v_max-v_min)/Nbins};
    vector<unsigned int> bins(Nbins,0);

    //4. Compute the bins in a parallel way
    #pragma omp parallel num_threads(NumThreads) if(Nsim/NumThreads>=PAR_THR_INDEX)
    {
        vector<unsigned int> local_bins(Nbins,0);

        #pragma omp for nowait
        for(size_t i=0;i<data.size();i+=step){
            float v{data[i+axis]};
            if(v>=v_min && v<=v_max){
                unsigned int b_index{static_cast<unsigned int>((v-v_min)/delta)};
                if(b_index==Nbins && v==v_max) b_index = Nbins -1;
                if(b_index<Nbins) local_bins[b_index]++;
            }
        }

        #pragma omp critical
        {
            for(size_t j=0;j<Nbins;j++){
                bins[j] += local_bins[j];
            }
        }
         
    }

    //5. Forming the output
    vector<vector<float>> res(Nbins,vector<float>(2,0.0f));
    float norm_factor{1.0f/(Nsim*delta)};

    for(unsigned int i=0;i<Nbins;i++){
        res[i][0] = v_min + (i + 0.5f)*delta;
        res[i][1] = bins[i]*norm_factor;
    }

    return res;
} 

//This function will produce starting from a Time Picutre such the one produced by "produceTimePicture" a 2D bin
//system useful to obtain PDFs, This function require:
//- A TimePicture (MANDATORY).
//- A 2D vector for the number of bins along the axis (MANDATORY).
//- A 2D vector of the axis along which doing the bins (0: times,n: the n variable of the system) (MANDATORY).
//- A bool to express if the binning domain is given or adaptive (true = adaptive).
//- If false a vector with upper and lower domain is required. It should be a 4 slot vector (lb,ub,lb,ub).
//The function will return a 2D vector with in each row the central value coordinates (firsts 2 column) and the bin value (last column).
//If there are too much threads than simulations (less than 10 simulation per thread) the code will be executed serially!!!
vector<vector<float>> CompFuncManager::PDF_2D(const TimePicture& picture,vector<unsigned int> Nbins,vector<unsigned int> axis,
                            bool adaptive, vector<float> domain)
{
    //0. Perform some checks on the arguments:
    checkFunctionPDF_2D(Nbins,axis,adaptive,domain);

    //1. Import data and features
    const vector<float>& data{picture.getAllPoints()};
    const unsigned int Nsim{picture.getNumSim()};
    const unsigned int step{ref_size+1};

    //2. Compute the limits of the domain along the given axis
    float x_max,x_min,y_max,y_min;

    if(adaptive){
        x_max = y_max = -1e30f;
        x_min = y_min = 1e30f;

        #pragma omp parallel for num_threads(NumThreads) reduction(min:x_min,y_min) reduction(max:x_max,y_max) if(Nsim/NumThreads>=PAR_THR_INDEX)
        for(size_t i=0;i<data.size();i+=step){
            float vx{data[i+axis[0]]};
            float vy{data[i+axis[1]]};
            if(vx<x_min) x_min = vx;
            if(vx>x_max) x_max = vx;
            if(vy<y_min) y_min = vy;
            if(vy>y_max) y_max = vy;
        }
    }
    else{
        x_min = domain[0];
        x_max = domain[1];
        y_min = domain[2];
        y_max = domain[3];
    }

    //3. Setup the bins
    float delta_x{(x_max-x_min)/Nbins[0]};
    float delta_y{(y_max-y_min)/Nbins[1]};
    vector<unsigned int> bins(Nbins[0]*Nbins[1],0);

    //4. Compute the bins in a parallel way
    #pragma omp parallel num_threads(NumThreads) if(Nsim/NumThreads>=PAR_THR_INDEX)
    {
        vector<unsigned int> local_bins(Nbins[0]*Nbins[1],0);

        #pragma omp for nowait
        for(size_t i=0;i<data.size();i+=step){
            float vx{data[i+axis[0]]};
            float vy{data[i+axis[1]]};
            if(vx>=x_min && vx<=x_max && vy>=y_min && vy<=y_max){
                unsigned int x_index{static_cast<unsigned int>((vx-x_min)/delta_x)};
                unsigned int y_index{static_cast<unsigned int>((vy-y_min)/delta_y)};

                if(x_index==Nbins[0] && vx==x_max) x_index = Nbins[0]-1;
                if(y_index==Nbins[1] && vy==y_max) y_index = Nbins[1]-1;

                if(x_index<Nbins[0] && y_index<Nbins[1]){
                    local_bins[x_index*Nbins[1]+y_index]++;
                }
            }
        }

        #pragma omp critical
        {
            for(size_t j=0;j<bins.size();j++){
                bins[j] += local_bins[j];
            }
        }
         
    }

    //5. Forming the output
    vector<vector<float>> res(Nbins[0]*Nbins[1],vector<float>(3,0.0f));
    float norm_factor{1.0f/(Nsim*delta_x*delta_y)};

    for(unsigned int i=0;i<Nbins[0];i++){
        for(unsigned int j=0;j<Nbins[1];j++){
            size_t index = i*Nbins[1] + j; 
            res[index][0] = x_min + (i + 0.5f)*delta_x;
            res[index][1] = y_min + (j + 0.5f)*delta_y;
            res[index][2] = bins[index]*norm_factor;
        }
    }

    return res;
}

//This function will produce, starting from a trajectory such the one produced by "produceTimePicture" the
//autocorrelation of the trajectory for a certain time delay. This function require:
//- A trajectory in style vector<vector<float>> (MANDATORY).
//- The axis (variable) along which computing the autocorrelation (MANDATORY).
//- The time delay (tau) of the autocorrelation (MANDATORY).
//The output will be the autocorrelation value computed as covariance/variance.
//Also the time delay is converted in the nearest below number of steps.
//If there are too much threads than simulations (less than 10 simulation per thread) the code will be executed serially!!!
float CompFuncManager::computeAutocorrelation(const Traj& traj,unsigned int axis,float tau){
    checkAutocorrelationInput(traj,axis,tau);

    //1. We need to convert the time delay in the step delay
    size_t window{findTimeIndex(traj.getTimes(),tau)};

    //2. We need then to compute the time mean of the trajectory for the axis
    float mean{0.0f};
    #pragma omp parallel for reduction(+:mean) if(traj.getLenght()/NumThreads>=PAR_THR_INDEX)
    for(size_t i=0;i<traj.getLength();i++){
        mean += traj.getVars()[traj.flatNotation(i,axis)];
    }

    mean = mean/traj.getLength();

    //3. We can compute the pseudo-covariance and the pseudo-variance at the same time for the common part
    float cov{0.0};
    float var{0.0};
    #pragma omp parallel for reduction(+:cov,var) if((traj.getLength()-1-window)/NumThreads>=PAR_THR_INDEX)
    for(size_t i=0;i<traj.getLength()-1-window;i++){
        cov += (traj.getVars()[traj.flatNotation(i,axis)]-mean)*(traj.getVars()[traj.flatNotation(i+window,axis)]-mean);
        var += (traj.getVars()[traj.flatNotation(i,axis)]-mean)*(traj.getVars()[traj.flatNotation(i,axis)]-mean);
    }
    //We need also to complete the variance
    float temp_var{0.0f};
    #pragma omp parallel for reduction(+:temp_var) if((1+window)/NumThreads>=PAR_THR_INDEX)
    for(size_t i=traj.getLength()-1-window;i<traj.getLength()-1;i++){
        temp_var += (traj.getVars()[traj.flatNotation(i,axis)]-mean)*(traj.getVars()[traj.flatNotation(i,axis)]-mean);
    }

    return cov/(var+temp_var);
}
