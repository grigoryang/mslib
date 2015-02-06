MYSOURCE  = 

MYPROGS   = peptideGrapheneAssembly1.noMPI.forRoy assemblyOptimizer assemblyOptimizer4Craig generateCrystalLatticeForC60 wipeChainID C60-to-C60-gap confind

MYHEADERS = 

MYDIR=myProgs/gevorg

#CCOPTIM = g++ -Wall -Wno-sign-compare -O3 -msse3 -mfpmath=sse -funroll-loops -fopenmp -g
#CCOPTIM = g++ -Wall -Wno-sign-compare -O3 -msse3 -mfpmath=sse -funroll-loops -fopenmp
#CCOPTIM = mpiCC -Wall -Wno-sign-compare -O3 -msse3 -mfpmath=sse -funroll-loops -fopenmp
#CCOPTIM = mpiCC -Wall -Wno-sign-compare -O3 -msse3 -mfpmath=sse -funroll-loops -fopenmp -pg -g


# this is for the graphene/peptide stuff
#MSL_EXTERNAL_LIB_DIR=/export/apps/contrib/lib
#MSL_EXTERNAL_INCLUDE_DIR=/export/apps/contrib/include

#MSL_EXTERNAL_LIB_DIR=/share/apps/contrib/lib
#MSL_EXTERNAL_INCLUDE_DIR=/share/apps/contrib/include

MSL_EXTERNAL_LIB_DIR=/usr/local/lib
MSL_EXTERNAL_INCLUDE_DIR=/usr/local/include

#CCOPTIM = mpiCC -Wall -Wno-sign-compare -O3 -msse3 -mfpmath=sse -funroll-loops -fopenmp
#CCDEBUG = mpiCC -Wall -Wno-sign-compare -msse3 -mfpmath=sse -funroll-loops -fopenmp -g
