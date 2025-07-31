# NoisePlugs:

This sublibrary of SDE-SS contains the so-called **NoisePlugs**. The *NoisePlugs* are used for systems that requires specific kind of noise which can
only be easily computed adding a variable to the system. This procedure will be still necessary using Plugs and they are not mandatory (most of them 
can be introduced by you directly) but the implementation is already made thus this should speed up your coding pipeline. 
How to use them is described below.

## How to use NoisePlugs:

Using the *NoisePlugs* is very simple. They are actually a little bit more than a normal function but with the possibility of storing the values of
the constants they use and, therefore, eventually change them during your work. To use them you must declare an extra variable for your system and
define the update of this variable in *f_function_impl* and *g_function_impl* with the *deterministic_part* and *stochastic_part* of your *NoisePlug*.
It is harder to explain than it seems, therefore we have inserted an example code in this subdirectory that shows you how easy actually it is to use the 
NoisePlugs. Nonetheless, I have added a test code to show how to use the NoisePlugs in the *TestCodes* directory.

## What NoisePlugs have been defined?

Until now we have defined only the following NoisePlugs, however feel free to suggest possible additions writing to the email in the main README.

- Wiener Process (*WienerProcessNP*): This NoisePlug can allow you, adding an extra variable, to use this extra variable of your system to simulate
	in parallel with your field a Wiener Process. This *NoisePlug* is created to be used in systems with Sine Wiener noises or similar which cannot
	be rewritten in a Ito's formula without the presence of *W(t)*. 