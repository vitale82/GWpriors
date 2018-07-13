Each hdf5 file contains posteriors samples obtained using the 8 priors described in Vitale et al, PRL  119, 251103 (2017)

The hdf5 file contains three groups: priors, posteriors, PSD (that is, power spectral density).

The priors and posteriors groups contain 8 subgroups, one for each of the priors P1.. P8 of Vitale et al.
So for example (using python):

import hdf5
data=h5py.File('GW150914.hdf5','r')
data['posteriors']['P1'] 

gives you access to the posterior distributions obtained using the P1 prior (that is the same prior used by the LIGO and Virgo collaborations).

We release the 15 parameters on which a CBC signal from a precessing binary in quasi-circular orbits depends. 
Those are (see LVC, PRL 116, 241102):

mc (chirp mass, Msun)
q  (mass ratio, in the range [0,1])
a1 (dimensionless spin for the most massive object, in the range [0,0.89])
a2 (dimensionless spin for the least massive object, in the range [0,0.89])
tilt1 (angle between the spin and the orbital angular momentum for the most massive object, at 20Hz, in rads)
tilt2 (angle between the spin and the orbital angular momentum for the least massive object, at 20Hz, in rads)
phi12 (difference of spins azimuthal angles in the plane of the orbit, at 20Hz, in rads)
phi_jl (angle between the total and orbital angular momentum, in rads)

distance (luminosity distance, Mpc)
psi (polarization angle, at 20Hz, rads)
phase (coalescence phase, rads)
time (gps arrival time at the geocenter, secs)
ra (rightascension, rads)
dec (declination, rads)
theta_jn (angle between the total angular momentum and the line of sight, rads)

We also release some useful derived parameters:

chi_eff (effective inspiral spin, in the range [-0.89,0.89])
chi_p (precessing spin, in the range [0,0.89])
m1 (mass of the most massive object, Msun)
m2 (mass of the least massive object, Msun)

The masses above are in the detector frame. We convert them to source-frame masses by dividing by (1+z).
We provide the following supplementary parameters:

m1_source (source-frame m1, Msun)
m2_source (source-frame m2, Msun)
mc_source (source-frame chirp mass, Msun)
redshift  

The redshift is obtained from the gravitational-wave luminosity distance, assuming a Plack 2015 cosmology (Planck collaboration, A&A 594, A13)
The values of the cosmological parameters are stored as attributes of the 'posteriors' group:

for (k,v) in data['posteriors'].attrs.iteritems():
    print k,v

omega_matter 0.3065
omega_lambda 0.6935
hnot 67.9
omega_k 0.0
w2 0.0
w1 0.0
w0 -1.0


Finally, we release the log likelihood and the logprior

logl (natural logarithm of the likelihood)
logp (natural logarithm of the prior)

The units of each parameter are given as attribute of the corresponding dataset. For example, the following command produces a plot of the luminosity distance posterior of GW150914 obtained with the prior P3

data=h5py.File('GW150914.hdf5','r')
dist=data['posteriors']['P3']['distance']
units=dist.attrs['units']
hist(dist,50)
xlabel('Distance [%s]'%units)
ylabel('PDF')

These posteriors were obtained using the nested sampling flavor of LALInference (Veitch et al, PRD 91, 042003).

The nested sampling computes the evidence of the data under the signal hypothesis, that is the likelihood marginalized over all parameters.
We report the natural log of the evidence of the signal hypothesis under all priors as an attribute of the posterior/PX group, e.g.

lognat_signal_evidence=data['posteriors']['P3'].attrs['signal_evidence']

and the natural log of the noise hypothesis

lognat_noise_evidence=data['posteriors']['P3'].attrs['noise_evidence']

The (natural log of the) Bayes Factor signal vs noise (e.g LVC, PRL 116, 241102, Table I) can thus be calculated as

BayesSN=lognat_signal_evidence-lognat_noise_evidence

Finally, the natural log of the Bayes factor between two different priors can be otained as e.g.

p1_signev=data['posteriors']['P1'].attrs['signal_evidence']
p5_signev=data['posteriors']['P5'].attrs['signal_evidence']
BayesP1P5=p1_signev-p5_signev # = 3.84

which shows P1 is favoured by the data

=======

The "priors" group contains samples drawn from the priors P1..P8 (each in the corresponding group).

To access e.g. the chi_eff P8 prior for GW150914

chi_eff_prior=data['priors']['P8']['chi_eff']


======


Finally, the 'PSD' group contains one group for each of the interferometers used in the analysis, H and L.
The PSD were generated with the BayesWave algorithm (Cornish and Littenberg, CQG 32, 13).

There are two datasets for each interferometer, 'frequency' and 'power_spectrum'
Those will be different for different events, but the same for all priors.

Thus, to plot the power spectrum used in the analysis of GW150914:

h1psd=data['PSD']['H']
l1psd=data['PSD']['L']
loglog(h1psd['frequency'],h1psd['power_spectrum'],'r',label='H')
loglog(l1psd['frequency'],l1psd['power_spectrum'],'b',label='L')



