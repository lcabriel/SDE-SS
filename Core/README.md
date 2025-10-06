# A small documentation for you:

This library has still not a proper documentation file that explains all the features and the function (it will be done in the future).
Nonetheless, in this README I will provide a little description of the class methods (from now on called "functions") and how to use it.
Only the public methods will be described but, if you like to implement something on your own, I had profusely left comment in the codes.

## Generic NoiseClass ("NoiseClass"):

This class is used to represent the noise aspect of the fields. This class is meant to be virtual and the noise classes such as *WienerEuler* and
*WienerMilstein* (see below) are child of this one. Feel free to implement your own *NoiseClass*-es: to do this the only function you need to
override is

- *vector&lt;float&gt; compute_noise(const vector&lt;float&gt; &x_i,float\* h))*: this class is used by the field to compute the noise. Here you need to
	put the actual "form" of your noise. As example you can look to *WienerMilstein*. The class can accept the point *x_i* if your noise is
	based on the coordinates or if you are reserving a variable to keep track of the noise (e.g. Ornstein-Uhlenbeck). It will require also
	the step length of your simulation's step.

## Wiener noise with Euler-Maruyama ("WienerEuler"):

This child of *NoiseClass* implements the gaussian white noise or Wiener process increment using the Euler-Maruyama approach. The scalability
of the noise with a constant is leaved for the *g_function* thus this class do not require any parameter upon construction. The random engine
is a mt19937 with a seed based on the instantiation time and the pointer of the entity.

- *vector&lt;float&gt; compute_noise(const vector&lt;float&gt; &x_i,float\* h))*: the virtual function is override in this version and the function will
	return a vector with the size equal to the dimensionality of the system with, in every slot, a *dW* computed separately in case you have
	multiple SDEs or you required multiple increments for some reasons.

## Wiener noise with Milstein ("WienerMilsten"):

Conceptually similar than before this class do the same thing as the *WienerEuler* but using the Milstein method. Due to this, it requires as
an argument a function which should be the derivative of your *g_function* that you use in the field. The other features are the same as
Wiener-Euler.

- *vector&lt;float&gt; compute_noise(const vector&lt;float&gt; &x_i,float\* h))*: as before, the function will return a vector of increments but computed
	using the Milstein method instead of the Euler-Maruyama.

## Generic FieldClass ("FieldClass"):

This class is used to represent your fields. This class it is meant as a base that you have to personalize via child class to express your system's
field. Upon initialization it requires a pointer to a *NoiseClass* (or its children) instance and your child constructor should have this too. To
set then the noise you have to insert in your constructor the *setNoise* protected function passing the pointer. 
Upon implementation you only need to modify the virtual functions *f_function_impl* and *g_function_impl*, while the other functions are only used 
by the other classes of the library. To see how to override these two, please look to *test_code.cpp* in the main directory.

- *void f_function_impl(const vector&lt;float&gt; &x,float t,vector&lt;float&gt; &y) const*: to personalize your system you need to override this virtual function
	inside your field. *f_function* is used to express the deterministic function of the Ito's formula ($dx = f(x)dt+g(x)dW$). To understand how to override correctly, please look at the test code. However, this function will require a system's point and will return the f(x) evaluation vector of the same size 
	of the system.
- *void g_function_impl(const vector&lt;float&gt; &x,float t,vector&lt;float&gt; &y) const*: Same as the function above. You need to override it to describe the 
	*g(x)* part of Ito's formula.
- *vector&lt;float&gt; getNoise(const vector&lt;float&gt; &x_i,float\* h)*: Upon call, this function will call the *compute_noise* function of the *NoiseClass*
	of your field that you should have set at construction.

## Trajectory ("Traj"):

This class is used to represent a trajectory of a system of SDEs. This is the type of output produced by the *simulateTrajectory* function of a *SDE_SS_System* instance. This a purely object and it is made more for memory optimization, tidiness of the code and unambiguousness. The class contains three elements: the vector of the time instants considered in the simulation, a 2D matrix of the values of the variables with the rows symbolizing the different time instants and the number of steps. It can be built both passing another *Traj* element or passing separately the times vector and the 2D matrix (in this case some checks are performed). Once a *Traj* is built there are no way to change the internal values.

- *getTimes()*, *getVars()*, *getLenght()*: these 3 get functions return respectively the three internal variables of the trajectory.
- *vector&lt;float&gt; getInstant(size_t TI)*: given a certain time index, this function will return the features of the system at that index. The output is a vector 		with the time at the 0th slot and then the values of the variables in the same order in which are saved along the columns of the 2D matrix.

## Time Picture ("TimePicture"):

This class is used to represent the concept of "time picture": this is a set of N points taken from N simulated trajectory of the SDEs system at a given
time instant T. It is like taking a picture of the system situation after a certain time T collecting the values that the set of N trajectory had at that
moment. As the Trajectory class, this is purely an object class to manage easily th values of the trajectories at that time instant. It is the product
of the *computeTimePicture* function. In *TimePicture*, there are stored the values in a 2D array with along the rows the different simulations and along
the columns the variables of the SDEs system. There are also saved the time instant of reference and the number of simulations.

- *getAllPoints()*, *getNumSim()*, *getTimeInstant()*: these 3 get functions return respectively the three internal variables of the trajectory.
- *vector&lt;float&gt;& getPoint(size_t p)*: given a specific simulation by index with p, this function return the values of that simulation at the time instant of
	the time picture.

## Set of Points ("SetOfPoints"):

This class is used to optimize memory usage when you are not interested in all the intermediate points of a trajectory. This is very very similar in its
structure to a Traj elements but it is made to contain only a set of points of a trajectory. The structure of the two classes are very similar with the
only exception of not having the *getLenght* amd *getInstant* functions. (MAYBE IN THE FUTURE WILL BE UNIFIED AS CHILD OF A COMMON CLASS...)

## Data Linker ("DataLinker"):

This small generic class is implemented when your fields presents some parameters or elements which are not fixed in time and are based on the values of a trajectory
*Traj* as the one produced by the *SDE_SS_System::simulateTrajectory*. This class (as *FieldClass*) has to be personalize via child
class giving a definition to the virtual function *getData* which will be use to retrive the data from the array when computing the field. This class has already implemented a constructor that accept a *Traj* instance or a time vector and a 2D values with in each row a time instant and the along the columns the values of the variable (in this case a same dimension check on the two arguments is performed) but you can expand these on your children classes. A DataLinker entity should be theoretically used inside your own custom FieldClass, however for more information about the implementation please look to the *test_code_data.cpp*. Both the array (called *data*) and the TimeOptimizationCounter (called *TOC*; see below) can be accessed by child classes.

- *float getData(const float t,const vector&lt;float&gt;& x)*: this virtual function, when defined in the child class, should gather one of the variables of the
	trajectory to the nearest previous time instant defined by *t*. Therefore, this function should be override by the child classes as shown in *test_code_data.cpp* having the possibility to access directly to the data array *data*. If necessary, you can create multiple copies of this function in your child if you need to access different data in different points. Obviously the array of data is expressed as a *Traj* trajectory, thus *t* has to be converted in an index to access on *Traj* arrays, calling in the child-defined *getData* the protected function *findTimeIndex*. The vector *x* is needed in case of an Adaptive DataLinker in which the data set expressed by the *Traj* can change accordingly to the system state *x*.
- *void setNewData(const Traj& t)*: this is the function to change properly the data set in case of an adaptive_set. The declaration if the DataLinker is adaptive or
	not is made upon construction specifying as second argument (as a boolean) if the DataLinker is adaptive or not. The *setNewData* function will check this upon call.
- *About the TimeOptimizationCounter*: this is not exactly a function but can be useful for the final user to know how the DataLinker works under the hood to exploit
	as better as possible its potential. The TimeOptimizationCounter (called *counter* in the code) is defined upon instancing with value equal to 0. This is used by the data linker to point where the previous time instant pointed in the trajectory. Usually, the data linked is used in the system simulations which evolve instant by instant the trajectory in time and the TimeOptimizationCounter is used by *findTimeIndex* to reduce the necessary time avoiding to repeat the searching process in its entirety. Moreover, I opted to let *counter* accessible in the class to allow the user in the new child-class function to avoid multiple searches if your process require to access to multiple variables of a given time instant in order. In any moment you can reset the counter to 0 with the public function *void resetTimeOptimizationCounter()*.   
- *float interpolData(const float& t,const unsigned int& variable)*: this protected function has to be used inside your custom *getData* function. It is not mandatory
	and it makes your code heavier under the computational aspect. Using this function, instead of picking the nearest previous data value according to your time instant *t*, the DataLinker will perform a linear interpolation between the previous and the next data value weighting with the time instant *t*. Calling the function you need to specify along which variable you want to interpolate (*variable*) with the usual notation of the *Traj* objects. *interpolData* can be useful when your data have a wide grid and there are a lot of point between the two same values during your simulation.   

## Stochastic Differential Equation - Super Solver - System ("SDE_SS_System"):

This is the main class of the library and the one you will use the most. It has three objective: representing your system, allowing you to produce 
trajectories and giving some useful tool to elaborate trajectories in an efficient way. We can split the public functions in three macro areas:

- **Core functions**: the fundamental ones, what you are going to use the most.
- **Public Utility functions**: small functions used to manage some aspect of your system and the simulations.
- **Tool functions**: also called *heavy*, this functions perform more complex tasks based on the simulations produced by the others.

We are going to explore then more deeply the three areas. The *SDE_SS_System* is also completely parallelized using OpenMp in the mechanics of the
heavy functions to grant good performances.
Upon construction, a *SDE_SS_System* will require as arguments the size of your system (a.k.a. the number of equations/variables), the pointer to
a FieldClass (or its children) instance and, optionally, a boolean value to express if the system is bounded (if it is, please take a look to *setBoundFunction*).

### Core functions:

The Core functions are the bulk of the *SDE_SS_System* and the functions you will use the most, probably. They are mostly linked to the "produce trajectories"
aspect of your system.

- *Traj simulateTrajectory(const vector&lt;float&gt; &x0,const float period,const float h_0)*: this function will simply produce a trajectory
	using an evolution method based on a Strang splitting treating the deterministic part with a RK4 method and the deterministic part depending on the *NoiseClass*
	*compute_noise* method of implementation (e.g. Euler-Maruyama for *WienerEuler*). The trajectory will be produced as a Traj object (see above).
	It requires as input a vector representing the initial condition, the length of the simulation and the base time step of the trajectory.
- *vector&lt;float&gt; simulateTrajectoryLastPoint(const vector&lt;float&gt; &x0,const float period,const float h_0)*: this function will act as simulate trajectory
	but it will not keep all the intermediate point. Under the runtime POV, this function uses lesser memory not saving the intermediate results. The output will
	not be a *Traj* element but only the last point of the trajectory expressed as (time,[variables]).
- *SetOfPoints simulateTrajectorySOP(const vector&lt;float&gt; &x0,const float period,const float h_0,const vector&lt;float&gt; &instants)*: a sort of in-between
	of the previous twos, this function will not keep all the intermediate values of the trajectories but only the ones indicated by the *instants* input. Actually
	it is difficult and consuming to obtain the exact values so the points immediately before is taken.

A little more articulate is the function to produce TimePictures

- *TimePicture produceTimePicture(\[TOO MANY\])*: this complex function produced a so-called Time Picture. A Time Picture is, given a set of trajectories, 
	the set of the values of the trajectories in a given time instant. This function requires a lot of arguments but produces its own trajectories internally. Therefore it needs:

	- The length of the simulations.
	- The base step of the simulations.
	- The number of simulations to produce.
	- A boolean value (optional) to say that the initial condition are random.
	- If the initial condition are not random, you need to pass the initial condition for each one in a (Nsim x init_cond) vector&lt;vector&lt;float&gt;&gt;.
	- Otherwise, if the initial condition are random, you need to pass a function to compute them. It should be a void -&gt; vector&lt;float&gt; function.
	- At last, if you want a particular time instant inside of the trajectory you need to pass it to the function, otherwise the last step values
		are used.

	Doing this, this function will return a TimePicture instance.

### Public Utility functions:

These small functions are actually quite important and are implemented to set some important features of the system, such as the bounds of the process (if bounded)
and the number of threads used in the parallel operations.

- *void setBoundFunction(function&lt;bool(const vector&lt;float&gt;&)&gt; f)*: if, during the construction procedure of a *SDE_SS_System* instance, the system is
	characterized as "bounded", you have also to set the bound function manually using this function. This function is used to check if the system is
	within the limit at every new time step (if it is not, the evolution is tried again with a smaller step). Considering this, it should be a *vector&lt;float&gt; -&gt; bool*
	function.
- *void setNumThreads(unsigned int N)*: all the Tool functions (see below) are parallelized on a certain degree. The number of threads used for the parallelisation
	is set equal to **8** by default. However, it is possible to change the number of threads used by these heavy functions with this public utility function.
- *static size_t findTimeIndex(const vector&lt;float&gt; &times,float TI)*: this static function is actually made to be used by certain Tool functions however it does
	not depend on the system characteristics and can be useful. Given a vector of times this function will found the last slot of the vector which time value is
	immedialtely lesser than the passed time index *TI*.

### Tool functions:

The Tool functions are functions (most of them heavy) that elaborates trajectories and time pictures to obtain common type of results used in papers and in research. All of them exploit parallelization using OpenMP, therefore they should be able to work perfectly even linearly, obviously with worse performances.

- *vector&lt;vector&lt;float&gt;&gt; PDF_1D(\[TOO MANY\])*: the idea is fairly simple: given a time picture and chosen a variable, this function will produce a bins
	system of same width along that variable ("axis"). Certainly, there are some features to define with the arguments:

	- As first argument, you have to pass a time picture as the one produced by *produceTimePicture* (*TimePicture*).
	- The number of bins to create.
	- The axis along which the bins are created (0th: time, nth: the nth-variable).
	- A boolean variable to tell if the bins system should be adaptive. An adaptive bin domains means that the function will find the minimum and the maximum value
		along the variable and then will build the bins between the two.
	- If is not adaptable, a domain expressed as a vector of two slot should be passed. The first value is the lower extreme while the second is the upper one.

	This function will return a *vector&lt;vector&lt;float&gt;&gt;* with number of rows equal to the number of bins and on each row the central value and the bins value.

- *vector&lt;vector&lt;float&gt;&gt; PDF_2D(\[TOO MANY\])*: very similar to the previous one, this function make a bin not along an axis but on the plane described
	by two of the variables. The arguments are almost the same but all is duplicated.

	- The first argument is the same.
	- Now the number of bins is expressed by a *vector&lt;unsigned int&gt;* with in first position the number of bins along the first axis and in the second position
		the number of bins along the second axis.
	- As for the bins, now the axis are expressed with a *vector&lt;unsigned int&gt;* of two slot, one for each variable of the 2D plane.
	- The boolean variable works the same.
	- The domain vector is now made by four slots instead of two and it has shape \[Axis 1 LB,Axis 1 UB,Axis 2 LB,Axis 2 UB\].

	The output is similar but with an extra column between the two because the central value is now a point in the plane and its expressed by two coordinates.

- *float computeAutocorrelation(const Traj& traj,unsigned int axis,float tau)*: the computation of the autocorrelation is often
	necessary especially in situations of quasi-periodicity. This functions requires as input a trajectory, the axis (*axis=n-1* -&gt; nth variable) along which you want to compute the autocorrelation and for which time delay. The returning value will be the autocorrelation value for the trajectory for that time delay computed as covariance over variance.
