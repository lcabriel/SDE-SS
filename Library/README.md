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

- *vector&lt;float&gt; f_function_impl(const vector&lt;float&gt; &x)*: to personalize your system you need to override this virtual function inside your field.
	*f_function* is used to express the deterministic function of the Ito's formula ($dx = f(x)dt+g(x)dW$). To understand how to override correctly,
	please look at the test code. However, this function will require a system's point and will return the f(x) evaluation vector of the same size 
	of the system.
- *vector&lt;float&gt; g_function_impl(const vector&lt;float&gt; &x)*: Same as the function above. You need to override it to describe the *g(x)* part of Ito's
	formula.
- *vector&lt;float&gt; getNoise(const vector&lt;float&gt; &x_i,float\* h)*: Upon call, this function will call the *compute_noise* function of the *NoiseClass*
	of your field that you should have set at construction.

## Stochastic Differential Equation - Super Solver - System ("SDE_SS_System"):

This is the main class of the library and the one you will use the most. It has three objective: representing your system, allowing you to produce 
trajectories and giving some useful tool to elaborate trajectories in an efficient way. We can split the public functions in three macro areas:

- **Core functions**: the fundamental, what you are going to use the most.
- **Public Utility functions**: small functions used to manage some aspect of your system and the simulations.
- **Tool functions**: also called *heavy*, this functions perform more complex tasks based on the simulations produced by the others.

We are going to explore then more deeply the three areas. The *SDE_SS_System* is also completely parallelized using OpenMp in the mechanics of the
heavy functions to grant good performances.
Upon construction, a *SDE_SS_System* will require as arguments the size of your system (a.k.a. the number of equations/variables), the pointer to
a FieldClass (or its children) instance and, optionally, a boolean value to express if the system is bounded (if it is please look to *setBoundFunction*).

### Core functions:

The Core functions are the bulk of the *SDE_SS_System* and the function you will use the most, probably. They are mostly linked to the "produce trajectories"
aspect of your system.

- *vector&lt;vector&lt;float&gt;&gt; simulateTrajectory(const vector&lt;float&gt; &x0,float period, float h_0)*: this function will simply produce a trajectory
	using an evolution method based on a Strang splitting treating the deterministic part with a RK4 method and the deterministic part depending on the *NoiseClass*
	*compute_noise* method of implementation (e.g. Euler-Maruyama for *WienerEuler*). The trajectory will be represented as a *vector&lt;vector&lt;float&gt;&gt;*
	with a the time steps along the row and the point variables along the column. Actually, the first column is reserved for the time.
	It requires as input a vector representing the initial condition, the length of the simulation and the base time step of the trajectory.
- *vector&lt;vector&lt;float&gt;&gt; produceTimePicture(\[TOO MANY\])*: this complex function produced a so-called Time Picture (maybe in the future it will have
	its own class). A Time Picture is, given a set of trajectories, the set of the values of the trajectories in a given time instant. This function requires a
	lot of arguments but produces its own trajectories internally. Therefore it needs:

	- The length of the simulations.
	- The base step of the simulations.
	- The number of simulations to produce.
	- A boolean value (optional) to say that the initial condition are random.
	- If the initial condition are not random, you need to pass the initial condition for each one in a (Nsim x init_cond) vector&lt;vector&lt;float&gt;&gt;.
	- Otherwise, if the initial condition are random, you need to pass a function to compute them. It should be a void -&gt; vector&lt;float&gt; function.
	- At last, if you want a particular time instant inside of the trajectory you need to pass it to the function, otherwise the last step values
		are used.

	Doing this, this function will return a vector&lt;vector&lt;float&gt;&gt; with in each row the values of the time and the variables at the time instant of one of the
	trajectories.

### Public Utility functions:

This small functions are actually quite important and are implemented to set some important features of the system such as the bounds of the process (if bounded)
and the number of threads used in the parallel operations.

- *void setBoundFunction(function&lt;bool(const vector&lt;float&gt;&)&gt; f)*: if, during the construction procedure of a *SDE_SS_System* instance, the system is
	characterized as "bounded", you have also to set the bound function manually using this function. This function is used to check if the system is
	within the limit at every new time step (if it is not the evolution is tried again with a smaller step). Considering this, it should be a *vector&lt;float&gt; -&gt; bool*
	function.
- *void setNumThreads(unsigned int N)*: all the Tool functions (see below) are parallelized on a certain degree. The number of threads used for the parallelisation
	is set equal to 8 by default. However, it is possible to change the number of threads used by these heavy functions with this public utility function.
- *static size_t findTimeIndex(const vector&lt;float&gt; &times,float TI)*: this static function is actually made to be used by certain Tool functions however it does
	not depend on the system characteristics and can be useful. Given a vector of times this function will found the last slot of the vector which time value is
	lesser than the passed time index *TI*. 
- *static vector&lt;float&gt; extractTimes(const vector&lt;vector&lt;float&gt;&gt; &traj)*: similarly to the previous one in terms of usage inside the library, this
	function can be useful because, given a trajectory-like *vector&lt;vector&lt;float&gt;&gt;* this function will extract the times as a *vector&lt;float&gt;*
	(a.k.a. the 0th column).

### Tool functions:

The Tool functions are functions (most of them heavy) that elaborates trajectories and time pictures to obtain common type of results used in papers and in research. All of them exploit parallelization using OpenMP therefore they should be able to work perfectly even linearly, obviously with worse performances.

- *vector&lt;vector&lt;float&gt;&gt; PDF_1D(\[TOO MANY\])*: the idea is fairly simple: given a time picture and chosen a variable this function will produce a bins
	system of same width along that variable ("axis"). Certainly, there are some features to define with the arguments:

	- As first argument, you have to pass a time picture as the one produced by *produceTimePicture* (*vector&lt;vector&lt;float&gt;&gt;*).
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

- *float computeAutocorrelation(const vector&lt;vector&lt;float&gt;&gt; &traj,unsigned int axis,float tau)*: the computation of the autocorrelation is often
	necessary especially in situations of quasi-periodicity. This functions requires as input a trajectory, the axis along which you want to compute the
	autocorrelation and for which time delay. The returning value will be the autocorrelation value for the trajectory for that time delay computed as covariance
	over variance.
