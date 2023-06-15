#Physical constants
ETA = 3.6e-3 #Blood viscosity
BE = 7/3 #Bifurcation exponent H-K model
#Model specific constants(boundary conditions)
Q_PERF = 8.33e-6 #Perfusion flow rate
P_PERF = 1.33e4 #Perfusion pressure
P_TERM = 8.38e3 #Terminal pressure
R_PERF = 0.05 # Radius of perfusion
N_TERM =  250 #Number of terminal segments, decides the resolution of the tree
SEGMENT_ASPECT_RATIO = 2

WALL_THICKNESS_LUMEN_RATIO = 0.1
EPSILON = 1e-3

Q_TERM = Q_PERF /N_TERM

OPENSCAD_EXPORT_SCALING_FACTOR = 1000

CORE_COUNT=1
