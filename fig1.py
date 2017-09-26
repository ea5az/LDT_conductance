import nest
import pylab as py
import numpy as np 

## NEST configuration
nest.ResetKernel()
nest.SetKernelStatus({"overwrite_files": True,
                      "data_path": "",
                      "data_prefix": "",
                      "resolution":0.01})
## Constants from paper 
# Kuhn, Alexandre, Ad Aertsen, and Stefan Rotter. 
# "Neuronal integration of synaptic input in the fluctuation-driven regime." 
# Journal of Neuroscience 24.10 (2004): 2345-2356.
N = 100
A = np.array([[1178 , 1],[100000,1]])
b = np.array([[0] , [52149]])
I_e = 0.
E_in = -75.
V_th = -50.

## Set up neuron populations
neuron = nest.Create("iaf_cond_alpha",N,params={"I_e":I_e,"E_in":E_in,"V_th":V_th})
noise_ex = nest.Create("poisson_generator")
noise_in = nest.Create("poisson_generator")
syn_dict_ex = {"weight": 7.1}
syn_dict_in = {"weight": -3.7}
nest.Connect(noise_ex, neuron, syn_spec=syn_dict_ex)
nest.Connect(noise_in, neuron, syn_spec=syn_dict_in)

## Determine linear coefficients for constant mean free membrane potential
coeff = np.linalg.solve(A,b)
E_rates = np.logspace(np.log10(1200),np.log10(100000),100)

## Iterate over rates and determine firing frequency
frates = np.zeros(E_rates.shape)
for ii , erate in enumerate(E_rates):
	print(ii)
	nest.SetStatus(noise_ex, {"rate": erate})
	nest.SetStatus(noise_in, {"rate": np.double(coeff[0]*erate + coeff[1])})
	## Reset the spike detector
	spikedetector = nest.Create("spike_detector",
	            	params={"withgid": True, "withtime": True})
	nest.Connect(neuron, spikedetector)
	nest.Simulate(1000.0)

	dSD = nest.GetStatus(spikedetector,keys="events")[0]
	evs = dSD["senders"]
	print(evs)
	frates[ii] = len(evs)*1.0/N

py.semilogx(E_rates , frates)
py.xlabel(r'$\lambda_e (1/s)$')
py.ylabel(r'Firing rate $(s^{-1})$') 
py.grid()
py.savefig('foo.png')
py.show()
