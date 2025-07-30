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

//##############################################################################################################################################################
//############################################################# GENERIC NOISE CLASS ########################################################################################
//##############################################################################################################################################################

//A standard compute noise that simply gives back the passed point in the system space.
vector<float> NoiseClass::compute_noise(const vector<float>& x_i,float* h){
    return x_i;
}


//##############################################################################################################################################################
//############################################################# PREDEFINED WIENER PROCESS ################################################################################
//##############################################################################################################################################################

//COSTRUCTOR: the distribution and the generator are initialized.
//The seed of the random engine is created starting from the time and the point of the instance.
WienerEuler::WienerEuler():
    eng(static_cast<unsigned long>(time(0) * reinterpret_cast<uintptr_t>(this))),
    distribution(0.0, 1.0)
{}

//Compute a dW kind of step for the Wiener process using the Euler-Maruyama method. 
//Actually a vector of the same size of the system is returned and, in each slot, there is a different dW. 
vector<float> WienerEuler::compute_noise(const vector<float>& x_i,float* h){
    vector<float> noise_out(x_i.size(),0.0);

    for(size_t i=0;i<x_i.size();i++){ //dW -> N(0,1)*sqrt(h)
        noise_out[i] = distribution(eng)*sqrt(*h);
    }

    return noise_out;
}


//##############################################################################################################################################################
//############################################################# PREDEFINED WIENER PROCESS ################################################################################
//##############################################################################################################################################################

//CONTRUCTOR: the distribution and the generator are initialize. 
//Moreover, requires as argument the derivative of the g_function used in the field class.
//The seed of the random engine is created starting from the time and the point of the instance.
WienerMilstein::WienerMilstein(function<vector<float>(const vector<float>&)> derivative):
    eng(static_cast<unsigned long>(time(0) * reinterpret_cast<uintptr_t>(this))),
    distribution(0.0, 1.0)
{
    //Check if the argument is valid.
    try{
        if(!derivative) error("Derivative of g_function argument is not valid.");
    }
    catch(const runtime_error& e){
        cout << "WienerMilstein: Constructor: Runtime error:" << e.what() << endl;
        exit(301);
    }

    D_g = derivative;
}

//Compute a dW kind-of step for the Wiener process using the Milstein method. 
//Actually a vector of the same size of the system is returned and, at each step there is a different dW. 
vector<float> WienerMilstein::compute_noise(const vector<float> &x_i,float* h){
    vector<float> noise_out(x_i.size(),0.0);

    for(size_t i=0;i<x_i.size();i++){
        float dW = distribution(eng)*sqrt(*h);

        noise_out[i] = dW + (0.5f)*(D_g(x_i)[i])*(dW*dW-*h);
    }

    return noise_out;
}


//##############################################################################################################################################################
//########################################################### GENERIC FIELD CLASS ############################################################################
//##############################################################################################################################################################

//This public function allow to set the noise of the FieldClass in a controlled way.
void FieldClass::setNoise(NoiseClass* N){

    try{
        if(N==nullptr) error("Invalid NoiseClass passed to the function.");
    }
    catch(const runtime_error& e){
        cout << "FieldClass: setNoise: Runtime error: " << e.what() << endl;
        exit(102);
    }

    noise = N;
    noise_initialized = true;
}

//Function for the deterministic part of the field. To define, please
//override "f_function_impl" as in documentation.
void FieldClass::f_function(const vector<float> &x,float t,vector<float> &y){
    try{
        if(!noise_initialized) error("Noise in the FieldClass not initialized.");
    }
    catch(const runtime_error& e){
        cout << "FieldClass: f_function: Runtime error: " << e.what() << endl;
        exit(101);
    }
    f_function_impl(x,t,y); 
}

//Function for the stochastic part of the field. To define, please
//override "g_function_impl" as in documentation.
void FieldClass::g_function(const vector<float> &x,float t,vector<float> &y){
    try{
        if(!noise_initialized) error("Noise in a FieldClass not initialized.");
    }
    catch(const runtime_error& e){
        cout << "FieldClass: g_function: Runtime error: " << e.what() << endl;
        exit(101);
    }
    g_function_impl(x,t,y);
}

//Function for the deterministic part of the field. Override this as in
//documentation to implement your system.
void FieldClass::f_function_impl(const vector<float> &x,float t,vector<float> &y){}

//Function for the stochastic part of the field. Override this as in
//documentation to implement your system.
void FieldClass::g_function_impl(const vector<float> &x,float t,vector<float> &y){}

//This function will give the result of the compute_noise of the local NoiseClass.
vector<float> FieldClass::getNoise(const vector<float> &x_i,float* h){
    try{
        if(!noise_initialized) error("Noise in a FieldClass not initialized");
    }
    catch(const runtime_error& e){
        cout << "FieldClass: getNoise: Runtime error: " << e.what() << endl;
        exit(101);
    }

    return noise->compute_noise(x_i,h);
}

//##############################################################################################################################################################
//########################################################### DATALINKER #############################################################################
//##############################################################################################################################################################

//Given a certain time instant, the variable will find the a time-index near the value (the previous one). The TimeOptimizationCounter
//is used to improve performances.
size_t DataLinker::findTimeIndex(float t){
    //1. Try to look if is the same time
    if((t>=data[counter][0])&&(t<data[counter+1][0])){
        return counter;
    }
    //2. Try to find if it is after the counter
    for(size_t i=counter+1;i<data.size();i++){
        if((t>=data[i][0])&&(t<data[i+1][0])){
            counter = i;
            return i;
        }
    }
    //3. Otherwise is before
    for(size_t i=0;i<counter;i++){
        if((t>=data[i][0])&&(t<data[i+1][0])){
            counter = i;
            return i;
        }        
    }

    cout << "DataLinker: findTimeInstant: Runtime error: time instant not found" << endl;
    exit(401);
}

//Given a certain time instant the getData returns a child-defined float value near to the given time.
float DataLinker::getData(float t){
    return 0.0;
}

//This function is used to reset the TimeOptimizationCounter which is used to optimize the finding process of the time instant during the
//simulations. 
void DataLinker::resetTimeOptimizationCounter(){
    counter=0;
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
    try{
        if(F == nullptr) error("Passed FieldClass is invalid (nullptr).");
        if(N<1) error("Given size is <1. It has no sense.");
    }
    catch(const runtime_error& e){
        cout << "SDE_SS_System: Constructor: Runtime error: " << e.what() << endl;
        exit(201);
    }
}

//################### EVOLUTION FUNCTIONS ####################
//These functions are used to perform the evolutive aspects of the trajectory's computation.

//Given the previous point and the step length, this internal function is the core function to evolve 
//the last step in the new one of the trajectory. It will use the RK4 method.
//The idea is to use a setup in the Strang splitting way synergizing with RK4 and the noise method.
vector<float> SDE_SS_System::evolveTraj(const vector<float> &x_n,float h,float t,vector<vector<float>> &k){

    //1. We simulate the deterministic part of the field with RK4 for half step
    //Now we compute X^1 (the first intermediate step)
    vector<float> x1{x_n+RK4_method(x_n,h/2,t,k)*(h/2.0f)};

    //2. Add the stochastic part
    vector<float> g(size,0.0);
    field->g_function(x1,t,g);

    x1 = x1+g*(field->getNoise(x1,&h)); //Compute the noisy part; is multiplied than for g_function

    //3. Evolve the last half step with RK4 for half step
    x1 = x1+RK4_method(x1,h/2,t,k)*(h/2.0f);

    return x1;
}

//Given the starting point and the step, this function will return the RK4 update
//of the deterministic part of the field.
//REMEMBER: you have to multiply externally by h or your step eventually.
vector<float> SDE_SS_System::RK4_method(const vector<float> &x0,float h,float t,vector<vector<float>> &k){

    vector<float> new_x{x0};

    field->f_function(new_x,t,k[0]); //compute k1
    new_x = x0+(h/2.0f)*k[0];

    field->f_function(new_x,t,k[1]); //compute k2
    new_x = x0+(h/2.0f)*k[1];

    field->f_function(new_x,t,k[2]); //compute k3
    new_x = x0+h*k[2];

    field->f_function(new_x,t,k[3]); //compute k4

    return (k[0]+2.0f*k[1]+2.0f*k[2]+k[3])*(1.0f/6.0f);
}

//##################### CORE FUNCTIONS ##########################################

//This core function simulate a single trajectory given the initial conditions, the time period and the standard step.
//This function will return a 2D vector with all the points of the trajectory. 
//In the returning vector the times and the values are conjointed. The first column is used for times.
vector<vector<float>> SDE_SS_System::simulateTrajectory(const vector<float> &x0,float period,float h_0){
    checkTrajInput(x0,period,h_0); //Check if the inputs are meaningful

    //Setup the initial point of the trajectory
    vector<vector<float>> values;
    vector<float> times;

    values.push_back(x0);
    times.push_back(0.0);

    float h{h_0}; //The h used step-by-step. Can change to adapt to the boundary

    vector<vector<float>> ks(4,vector<float>(size,0.0)); //k vector used by the RK4 method to avoid reallocation

    do{

        vector<float> new_point{evolveTraj(values.back(),h,times.back(),ks)}; //define the new point

        //if out of bound
        if(bounded){
            while(!bounds(new_point)){
                h = h/10.0;
                new_point = evolveTraj(values.back(),h,times.back(),ks);
            }
        }

        //Add the new point
        values.push_back(new_point);
        times.push_back(times.back()+h);

        //Reset the step
        h = h_0;

    }while((times.back()+h)<=period);

    //Prepare the final output
    vector<vector<float>> output(times.size(),vector<float>(size+1,0.0));

    for(unsigned long int i=0;i<times.size();i++){
        output[i][0] = times[i];
        for(size_t j=0;j<size;j++){
            output[i][j+1] = values[i][j];
        }
    }

    return output;
}

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
vector<vector<float>> SDE_SS_System::produceTimePicture(float period,float h_0,unsigned int Nsim,
                                                        bool random_initial,const vector<vector<float>> &x0 ,
                                                        function<vector<float>()> random_f,float time_instant)
{
    //1. Perform some checks on the arguments:
    checkFunctionComputeTimePicture(Nsim,random_initial,x0,random_f);

    //2. Now we can start to produce the simulations and extract the points.
    //We will work using OpenMP to optimize the performances.
    vector<vector<float>> picture;

    //2.1 Assign at each threads its quota. The master will have also the remainder
    vector<unsigned int> quota(NumThreads,0);
    vector<unsigned int> cumulative_quota(NumThreads,0); //Used for initial conditions

    quota[0] = (Nsim / NumThreads) + (Nsim % NumThreads);
    cumulative_quota[0] = 0;
    for(size_t i=1;i<NumThreads;i++){
        quota[i] = Nsim / NumThreads;
        cumulative_quota[i] += quota[i-1] + cumulative_quota[i-1];
    }

    #pragma omp parallel num_threads(NumThreads)
    {
        const unsigned int ID = omp_get_thread_num(); //Get thread ID

        //2.2 We create the local storaging thus each thread can 
        vector<vector<float>> local_picture(quota[ID],vector<float>(size+1,0.0));

        //2.3 Start to produce the trajectories
        for(size_t i=0;i<quota[ID];i++){

            //2.3.1: We need to produce the initial condition
            vector<float> init;
            if(random_initial){ //random initial condition
                init = random_f();
            }
            else{ //non random initial condition
                init = x0[cumulative_quota[ID]+i];
            }

            //2.3.2: Compute the trajectory
            vector<vector<float>> traj{simulateTrajectory(init,period,h_0)};

            //2.3.3: Find the time index associated to the time instant and extract
            size_t time_index;
            if(time_instant<=0){ //last point
                time_index = traj.size()-1;
            }
            else{ //If is not the last we do it externally with two functions
                time_index = findTimeIndex(extractTimes(traj),time_instant);
            }

            local_picture[i] = traj[time_index];
        }

        //Reunite all the values

        #pragma omp critical
        {
            for(size_t i=0;i<quota[ID];i++){
                picture.push_back(local_picture[i]);
            }
        }

    }

    return picture;
}

//################### PUBLIC UTILITY FUNCTIONS ####################################

//This function will set the bound function. The bound function should return a boolean
//and has as input a vector<float> (the point). The function should return "true" when
//the point is in the domain.
//The system has to be bounded (construction) to this function to work.
void SDE_SS_System::setBoundFunction(function<bool(const vector<float>&)> f = nullptr){
    
    //Check that the system is defined as bounded and the function defined
    try{
        if(!bounded) error("SDE_SS_System is not bounded but setBoundFunction is called.");
        if(!f) error("setBoundFunction argument is not valid.");
    }
    catch(const runtime_error& e){
        cout << "SDE_SS_System: setBoundFunction: Runtime error: " << e.what() << endl;
        exit(203);
    }

    bounds = f;
}

//This function is used to set the number of threads used by the heavy functions of the class.
//The standard value is 8.
void SDE_SS_System::setNumThreads(unsigned int N){
    NumThreads = N;
}

//##################### TOOL FUNCTIONS ##########################################
//These function are made to automatize some typical tests and procedure used in the analysis of the SDE.

//This function will produce starting from a TimePicture such the one produced by "produceTimePicture" a 1D bin
//system useful to obtain PDFs. This function require:
//- The TimePicture style vector<vector<float>> (MANDATORY).
//- The number of bins (MANDATORY).
//- The axis along which doing the bins (0: times,n: the n variable of the system) (MANDATORY).
//- A bool to express if the binning domain is given or adaptive (true = adaptive).
//- If false a vector with upper and lower domain is required.
//The function will return a 2D vector with in each row the central value (first column) and the bin value (second column).
vector<vector<float>> SDE_SS_System::PDF_1D(const vector<vector<float>> &picture,unsigned int Nbins,unsigned int axis,
                                            bool adaptive,vector<float> domain)
{
    //0. Create the bin vector
    vector<vector<float>> bins(Nbins,vector<float>(2,0.0));

    //1. Perform some checks on the arguments:
    checkFunctionPDF_1D(picture,Nbins,axis,adaptive,domain);

    //2. If the domain is adaptive we need to find the extremes
    vector<float> extremes(2,0.0);
    
    if(adaptive){
        extremes[0] = picture[0][axis]; //The min
        extremes[1] = extremes[0]; //The max

        for(unsigned long int i=0;i<picture.size();i++){
            if(picture[i][axis]<extremes[0]) extremes[0]=picture[i][axis];
            if(picture[i][axis]>extremes[1]) extremes[1]=picture[i][axis];
        }
    }
    else{ //Domain is passed
        extremes = domain;
    }

    //3. Produce the bins limits
    vector<float> limits(Nbins,0.0);
    const float delta{(extremes[1]-extremes[0])/Nbins};

    for(size_t i=0;i<Nbins-1;i++){
        limits[i] = extremes[0] + (i+1)*delta;
    }
    limits[Nbins-1] = extremes[1];

    //4. Produce the bins mid-values
    for(size_t i=0;i<Nbins;i++){
        bins[i][0] = limits[i] - delta/2;
    }

    //5. We can start to work on the bins. It will be parallelized with OpenMP.
    //5.1 This time we need to compute the quotas before:
    vector<int> quotas(NumThreads+1,0);
    quotas[0] = 0;
    quotas[1] = picture.size()/NumThreads + picture.size()%NumThreads;
    for(size_t i=1;i<NumThreads;i++){
        quotas[i+1] = quotas[i] + picture.size()/NumThreads;
    }

    //5.2 We can work on the bins in a parallel way
    #pragma omp parallel num_threads(NumThreads)
    {
        const unsigned int ID = omp_get_thread_num();

        vector<float> local_bins(Nbins,0.0);

        for(size_t i=quotas[ID];i<quotas[ID+1];i++){
            for(size_t j=0;j<Nbins;j++){
                if(picture[i][axis]<=limits[j]){
                    local_bins[j] += 1.0;
                    break;
                }
            }
        }

        #pragma omp critical
        {
            for(size_t j=0;j<Nbins;j++){
                bins[j][1] += local_bins[j];
            }
        }
         
    }

    return bins;

}  

//This function will produce starting from a Time Picutre such the one produced by "produceTimePicture" a 2D bin
//system useful to obtain PDFs, This function require:
//- The TimePicture style vector<vector<float>> (MANDATORY).
//- A 2D vector for the number of bins along the axis (MANDATORY).
//- A 2D vector of the axis along which doing the bins (0: times,n: the n variable of the system) (MANDATORY).
//- A bool to express if the binning domain is given or adaptive (true = adaptive).
//- If false a vector with upper and lower domain is required. It should be a 4 slot vector (lb,ub,lb,ub).
//The function will return a 2D vector with in each row the central value coordinates (firsts 2 column) and the bin value (last column).
vector<vector<float>> SDE_SS_System::PDF_2D(const vector<vector<float>> &picture,vector<unsigned int> Nbins,vector<unsigned int> axis,
                            bool adaptive, vector<float> domain)
{
    //0. Create the bin vector
    vector<vector<float>> bins(Nbins[0]*Nbins[1],vector<float>(3,0.0));

    //1. Perform some checks on the arguments:
    checkFunctionPDF_2D(picture,Nbins,axis,adaptive,domain);

    //2. If the domain is adaptive we need to find the extremes
    vector<float> extremes(4,0.0);
    
    if(adaptive){
        extremes[0] = picture[0][axis[0]]; //The min
        extremes[1] = extremes[0]; //The max
        extremes[2] = picture[0][axis[1]]; //The min
        extremes[3] = extremes[2]; //The max

        for(unsigned long int i=0;i<picture.size();i++){
            if(picture[i][axis[0]]<extremes[0]) extremes[0]=picture[i][axis[0]];
            if(picture[i][axis[0]]>extremes[1]) extremes[1]=picture[i][axis[0]];
            if(picture[i][axis[1]]<extremes[2]) extremes[2]=picture[i][axis[1]];
            if(picture[i][axis[1]]>extremes[3]) extremes[3]=picture[i][axis[1]];
        }
    }
    else{ //Domain is passed
        extremes = domain;
    }

    //3. Produce the bins limits
    vector<float> limits_0(Nbins[0],0.0);
    vector<float> limits_1(Nbins[1],0.0);
    const float delta_0{(extremes[1]-extremes[0])/Nbins[0]};
    const float delta_1{(extremes[3]-extremes[2])/Nbins[1]};

    for(size_t i=0;i<Nbins[0]-1;i++){
        limits_0[i] = extremes[0] + (i+1)*delta_0;
    }
    limits_0[Nbins[0]-1] = extremes[1];

    for(size_t i=0;i<Nbins[1]-1;i++){
        limits_1[i] = extremes[2] + (i+1)*delta_1;
    }
    limits_1[Nbins[1]-1] = extremes[3];

    //4. Produce the bins mid-values
    for(size_t i=0;i<Nbins[0];i++){
        for(size_t j=0;j<Nbins[1];j++){
            bins[i*Nbins[1]+j][0] = limits_0[i] - delta_0/2;
        }
    }

    for(size_t i=0;i<Nbins[1];i++){
        for(size_t j=0;j<Nbins[0];j++){
            bins[j*Nbins[1]+i][1] = limits_1[i] - delta_1/2;
        }
    }

    //5. We can start to work on the bins. It will be parallelized with OpenMP.
    //5.1 This time we need to compute the quotas before:
    vector<int> quotas(NumThreads+1,0);
    quotas[0] = 0;
    quotas[1] = picture.size()/NumThreads + picture.size()%NumThreads;
    for(size_t i=1;i<NumThreads;i++){
        quotas[i+1] = quotas[i] + picture.size()/NumThreads;
    }

    //5.2 We can work on the bins in a parallel way
    #pragma omp parallel num_threads(NumThreads)
    {
        const unsigned int ID = omp_get_thread_num();

        vector<float> local_bins(Nbins[0]*Nbins[1],0.0);

        for(size_t i=quotas[ID];i<quotas[ID+1];i++){
            for(size_t j=0;j<Nbins[0];j++){
                for(size_t k=0;k<Nbins[1];k++){
                    if((picture[i][axis[0]]<=limits_0[j])&&(picture[i][axis[1]]<=limits_1[k])){
                        local_bins[j*Nbins[1]+k] += 1.0;
                        break;
                    }
                }
                break;
            }
        }

        #pragma omp critical
        {
            for(size_t j=0;j<Nbins[0];j++){
                for(size_t k=0;k<Nbins[1];k++){
                    bins[j*Nbins[1]+k][2] += local_bins[j*Nbins[1]+k];
                }
            }
        }
         
    }

    return bins;

}

//This function will produce, starting from a trajectory such the one produced by "produceTimePicture" the
//autocorrelation of the trajectory for a certain time delay. This function require:
//- A trajectory in style vector<vector<float>> (MANDATORY).
//- The axis (variable) along which computing the autocorrelation (MANDATORY).
//- The time delay (\tau) of the autocorrelation (MANDATORY).
//The output will be the autocorrelation value computed as covariance/variance.
//Also the time delay is converted in the nearest below number of steps.
float SDE_SS_System::computeAutocorrelation(const vector<vector<float>> &traj,unsigned int axis,float tau){
    checkAutocorrelationInput(traj,axis,tau);

    //1. We need to convert the time delay in the step delay
    size_t window{findTimeIndex(extractTimes(traj),tau)};

    //2. We need then to compute the time mean of the trajectory for the axis
    float mean{0.0};
    for(size_t i=0;i<traj.size();i++){
        mean += traj[i][axis];
    }
    mean = mean/traj.size();

    //3. We can compute the pseudo-covariance and the pseudo-variance at the same time for the common part
    float cov{0.0};
    float var{0.0};
    for(size_t i=0;i<traj.size()-1-window;i++){
        cov += (traj[i][axis]-mean)*(traj[i+window][axis]-mean);
        var += (traj[i][axis]-mean)*(traj[i][axis]-mean);
    }
    //We need also to complete the variance
    for(size_t i=traj.size()-1-window;i<traj.size()-1;i++){
        var += (traj[i][axis]-mean)*(traj[i][axis]-mean);
    }

    return cov/var;
}

//##################################################### SUPPORT FUNCTIONS #########################################################################
//These functions are used to light the code of other more complex functions.

//This function will perform the checks of the parameters of computeTimePicture.
void SDE_SS_System::checkFunctionComputeTimePicture(unsigned int Nsim,bool random_initial,const vector<vector<float>> &x0,function<vector<float>()> random_f){
    //1.0: we do not check period and h_0 because it is done in simulateTrajectory

    //1.1: we briefly check that Nsim is not 0, the function is defined and the size is meaningful.
    try{
        if(Nsim==0) error("produceTimePicture cannot accept Nsim=0.");
        if(random_initial && !random_f) error("random initial condition function not defined.");
        if(!random_initial && x0.size()!=Nsim && x0[0].size()!=size) error("initial condition has an invalid size.");
    }
    catch(const runtime_error& e){
        cout << "SDE_SS_System: computeTimePicture: Runtime error: " << e.what() << endl;
        exit(204);
    }
}

//This function will perform the checks of the parameters of simulateTrajectory. 
void SDE_SS_System::checkTrajInput(const vector<float> &x0,float period,float h_0){
    
    //Check correct size of the initial condition vector, the positivity of the period and the meaningfulness of the step
    try{
        if(x0.size()!=size) error("Initial condition vector have more or less dimensions than the system.");
        if(period<=0) error("Negative or null time period for the trajectory.");
        if((h_0<=0)||(h_0>period)) error("Unmeaningful time step for the trajectory.");
    }
    catch(const runtime_error& e){
        cout << "SDE_SS_System: simulateTrajectory: Runtime error: " << e.what() << endl;
        exit(202);
    }
}

//Given a vector of times and a time index this function will find the slot just before the given time instant.
//This function is also STATIC.
size_t SDE_SS_System::findTimeIndex(const vector<float> &times,float TI){
    for(size_t i=0;i<times.size();i++){
        if(TI<times[i]){
            return i-1;
        }
    }

    cout << "SDE_SS_System: findTimeIndex: Runtime error: given time instant is not present!" << endl;
    exit(205);
}

//Given a vector of the shape of the one of the simulateTrajectory function,
//this function will extract the column of times (column 0). This function is also STATIC.
vector<float> SDE_SS_System::extractTimes(const vector<vector<float>> &traj){
    vector<float> times(traj.size(),0.0);

    for(size_t i=0;i<traj.size();i++){
        times[i] = traj[i][0];
    }

    return times;
}

//This function will perform the checks of the parameters of PDF_1D.
void SDE_SS_System::checkFunctionPDF_1D(const vector<vector<float>> &picture,unsigned int Nbins,unsigned int axis,bool adaptive,const vector<float> &domain){
    //Check the picture size is equal to the system size, the number of bins meaningful,
    //the axis is not bigger than the dimensionality and the domain is valid.
    try{
        if(picture[0].size()!=(size+1)) error("Time Picture passed has a size incompatible with the problem.");
        if(Nbins<1) error("Number of bins passed is less than 1.");
        if(axis>=(size+1)) error("No variable associated with that axis of the system.");
        if(domain[0]==domain[1]) error("Undefined or domain with equal extremes.");
    }
    catch(const runtime_error& e){
        cout << "SDE_SS_System: PDF_1D: Runtime error: " << e.what() << endl;
        exit(205);
    }
}

//This function will perform the checks of the parameters of PDF_2D.
void SDE_SS_System::checkFunctionPDF_2D(const vector<vector<float>> &picture,vector<unsigned int> Nbins,vector<unsigned int> axis,bool adaptive,const vector<float> &domain){
    //1.1 Check the picture size is equal to the system size, the number of bins meaningful
    //the axis is not bigger than the dimensionality and the domain is valid. This is done for both the axis.
    try{
        if(picture[0].size()!=(size+1)) error("Time Picture passed has a size incompatible with the problem.");
        if(Nbins[0]<1 || Nbins[1]<1) error("Number of bins passed is less than 1.");
        if(axis[0]>=(size+1) || axis[1]>=(size+1)) error("No variable associated with that axis of the system.");
        if(domain[0]==domain[1]||domain[2]==domain[3]) error("Undefined or domain with equal extremes.");
    }
    catch(const runtime_error& e){
        cout << "SDE_SS_System: PDF_2D: Runtime error: " << e.what() << endl;
        exit(205);
    }
}

//This function will perform the checks of the parameters of computeAutocorrelation.
void SDE_SS_System::checkAutocorrelationInput(const vector<vector<float>> &traj,unsigned int axis,float tau){
    try{
        if(axis>=size+1) error("No variable associated with the given autocorrelation axis.");
        if(tau>traj.back()[0]) error("Time delay bigger than the trajectory length.");
        if(tau>0) error("Time delay should be bigger than 0.");
    }
    catch(const runtime_error& e){
        cout << "SDE_SS_System: computeAutocorrelation: Runtime error: " << e.what() << endl;
        exit(206);
    }
}
