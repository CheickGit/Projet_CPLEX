
SYSTEM = x86-64_linux
LIBFORMAT = static_pic
CONCERTDIR = /opt/ibm/ILOG/CPLEX_Studio2212/concert
CPLEXDIR = 	/opt/ibm/ILOG/CPLEX_Studio2212/cplex
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

#LIBTORCH = /home/cheick/libtorch

#CXXFLAGS += -I$(LIBTORCH)/include -I$(LIBTORCH)/include/torch/csrc/api/include
#LDFLAGS += -L$(LIBTORCH)/lib -Wl,-rpath,$(LIBTORCH)/lib -ltorch -lc10


#SYSTEM = x86-64_osx
#LIBFORMAT = static_pic
#CONCERTDIR = /Applications/CPLEX_Studio128/concert
#CPLEXDIR = 	/Applications/CPLEX_Studio128/cplex
#CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
#CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

# ---------------------------------------------------------------------
# Flags
# ---------------------------------------------------------------------
CPLEXFLAGS = -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include 


CC = g++
CCFLAGS = -O3 -DIL_STD
OBJECTS = main_modele1_lineaire.o modele1_lineaire.o

#---------------------------------------------------------
# Files
#---------------------------------------------------------
exec : $(OBJECTS)
	   $(CC) $(CCFLAGS) -o mcnd_mip $(OBJECTS) -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -ldl -lpthread
	   rm -f $(OBJECTS)
	   
.cpp.o : 
		$(CC) $(CCFLAGS) $(CPLEXFLAGS) -c $< -o $@	
		

		
	

