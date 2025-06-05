# SDE-SS (Stochastic Differential Equation - Super Solver) v. 1.1.0
This repository contains the C++ library SDE-SS made mainly to compute trajectories of systems of SDEs. Moreover, useful related function are also integrated.
The main part of the library is contained inside the "Core" subdirectory. Here there are all the most important and fundamental part that allow you to simulate
SDEs.
Nonetheless, the library provides also some plug-ins that are created to simplify your coding pipeline with some common structures that you can encounter 
working with systems of SDEs:

- **NoisePlugs**: this structures, a little more complex than a standard functions, are used to introduce rapidly common types of errors that requires to be
	treated as extra variables (e.g. Ornstein-Uhlenbeck, Wiener processes, ...).

## How to use:

To use the library in its core, you have to simply add at beginning of your .cpp the header writing

```
#include "SDE-SS.h"
```

The library is made completely on standard libraries therefore no other particulars are required.
You need then to copy the header and the .cpp file into your directory and compile your test file $test.cpp$ as

```
g++ test.cpp SDE-SS.cpp -o test.exe -lm -fopenmp
```

The descriptor *-lm* is needed to link the **cmath** standard library used to compute certain values. About *-fopenmp*, this is optional
but is strongly suggested to improve your performance especially when using heavy functions such as *PDF_1D*.

## Example of usage:

For a better example, please look to the file *test_code.cpp*.
The central class of the library is the *SDE_SS_System*, this class is the representation of your dynamic system.

To be characterised, a system requires a certain dimensionality and a field. The field has to be represented creating a child class of the general *FieldClass*.
The *FieldClass* requires also a Noise. We have already implemented a couple of noise in the library however custom one can be created modifying the general *NoiseClass*.

Thus defined a certain *field* as in the *test_code.cpp* we can build a 3D system as

```
SDE_SS_System system(3,&field,true);
```

The last boolean argument of the constructor is used to tell the system that the problem is bounded. To understand how to use the bounds, please look to *test_code.cpp* and to the README of the "Library" directory.

Defined the system producing a trajectory is simple:

```
vector<vector<float>> traj{system.simulateTrajectory(x0,period,h)};
```

where *x0* are the initial conditions, *period* is the length of the simulation and $h$ the base time step.

## Info and improvements:

This is the first bulk of the library and new features will be implemented. 

If you have some suggestions or function you would like to see in the library, please write me at:

EMAIL: lorenzo.cabriel@phd.units.it 



