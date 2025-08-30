### #This is a re-simulation of a published paper, "Neural heterogeneity controls computations in spiking neural networks."

While joining the Journal Club of the Computational Neuroscience Department, I had the opportunity to present this paper. 
I reproduced the simulations of the IK(Izhikevich) neuron network and the derived mean-field model described in the following paper:
R.Gast, S.A. Solla,& A. Kennedy, "Neural heterogeneity controls computations in spiking neural networks," Proc. Natl. Acad. Sci. U.S.A. 121 (3) e2311885121, https://doi.org/10.1073/pnas.2311885121 (2024).

This paper utilized the Lorentz distribution, which is equivalent to the Cauchy distribution.
Before directly using the equation in the paper, I mathematically formalized the Cauchy distribution and simulated it to check whether it follows the theoretical
probability density shape. Which is in the "Cauchy–Lorentz distribution.pdf" file.

In the Spiking_Neural_Networks.m file, I made parameters for Regular spiking neurons, Fast spiking neurons, and Low threshold spiking neurons. In each type of spiking neuron, I checked how much the results of the simulation differ between the SNN(Single Population Neuron Network) and the Mean-field model.

As we can confirm from the results of the PDF files, the simulations quite well follow the paper.
