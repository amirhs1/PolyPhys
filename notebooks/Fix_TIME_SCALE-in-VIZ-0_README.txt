
1. File naming convention in "animation" directory

The file naming convention in a file name like the one given below

	N200epshm29nh37ac1nc47747l25dt0.005ndump2000adump5000ens1.ring-frame0_500-speed0.5

is this:

N200: 200 DNA monomers with size (diameter) sigma=1.0
epshm29: 29k_bT is the LJ interaction strength between proteins poles and DNA monomers.
nh37: 37 HNS proteins.
ac1: 1.0sigma is the crowders size.
nc47747: 47747 is the number of crowders.
l25: 25sigma is half of the side of a cubic simulation box.
dt0.005: 0.005 is the integration time step in the MD algorithm in the production run.
ndump2000: 2000 is the frequency of dumping nucleoid particles (monomer+protein) configurations.
adump5000: 5000 is the frequency of dumping all particles (monomer+protein+crowder) configurations.
ens1: 1 is the id of the first independent run of this system (set of input parameters)
ring: Topology of the polymer
frame0_500: 0 to 500 is the range of frames (configurations converted to an animation.
speed0.5: 0.5 per second is the playback speed; the frame rate is still 32.


Note:

1 time step=dt=0.005tau_LJ where tau_LJ=1.0 is the LJ time unit
1 dumping time : 2000*dt=2000*0.005tau_LJ (nucleoid) or 5000*dt=5000*0.005tau_LJ (all) = 10tau_LJ or 25tau_LJ: In the plots (see below), tau=10tau_LJ.
Total simulation = 2*10^8dt=10^6tau_LJ=


Observation:

The animations/videos with frame0_500 show how the initial time interval, in which the proteins find and bind the polymer, differs from one system to another.

The animations/videos with frame21500_20000 show how the some time interval, in the middle of simulation, in which the proteins have already bridged different parts of the polymer and the polymer is intra-looped, differs from one system to another. 



2. File naming convention in "tseries-whole_range" and "tseries-2times10e7_first_timesteps" directories

The time evaluation of three physical properties of the ring polymer are shown, namely radius of gyration (R_g), shape parameters (S), and asphiricity (Delta). In my lingo, "space" means a group of simulations in which N(number of DNA monomers), epshm (the LJ interaction strength between proteins poles and DNA monomers), nh (number of HNS proteins) and ac (size of crowders) are fixed. Additionally, "ensemble" means N, epshm, nh, ac, and phic (volume fraction of crowders) are fixed.


These two directories are different in the length of a plotted time series. The whole simulation time (2*10^8) is depicted in the  "tseries-whole_range" directory while the the first 10^7 time step (1/20 of whole run) is plotted.


Each directory has 3 child directories:

	"per_space-per_ensemble-per_property_plot": Each file is for a given pairs of  space and physical property. In a given file, each subplot is for an ensemble.
	
	"per_space-per_property_plot": Each file is is for a given physical property with fixed ac=2.0. In a given file, each subplot is for a fixed number of HNS proteins (a space). Line
 	colors show different ensembles.

	"per_space-per_property_plot-z_score_normalization": Similar to "per_space-per_property_plot" directory. However, here, I normalized each physical property by the z-score; that is, I
	subtracted the mean from the time series of a given property and divided it by the variance. This way, I wanted to make the variations due to bridging or crowding more pronounced.

Observation:

3. "equilibriumProperties-mean" and "equilibriumProperties-norm" files

These files show the radius of gyration (R_g) shape parameters (S) and asphiricity (Delta) of the chain as a functions of the volume fraction of crowders. 


