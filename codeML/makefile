# ---------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------
SYSTEM = x86-64_linux
LIBFORMAT = static_pic
CPLEX_STUDIO = /opt/ibm/ILOG/CPLEX_Studio2212
CPLEXDIR = $(CPLEX_STUDIO)/cplex
CONCERTDIR = $(CPLEX_STUDIO)/concert
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

PYTHON = python3
VENV = .venv
VENV_PYTHON = $(VENV)/bin/python
PYTHON_INCLUDE = $(shell $(VENV_PYTHON) -c "from sysconfig import get_paths; print(get_paths()['include'])")
# Utiliser python3.11-config pour obtenir tous les indicateurs de liaison nécessaires pour l'embarquement
# Cela suppose que python3.11-config est dans le PATH et correspond à la version de Python utilisée.
# $(VENV_PYTHON)-config serait plus spécifique au venv, mais si le venv utilise le python système,
# python3.11-config (système) est généralement correct.
PYTHON_EMBED_LDFLAGS = $(shell python3.11-config --ldflags --embed)
# ---------------------------------------------------------------------
# Compilation flags
# ---------------------------------------------------------------------
CPLEXFLAGS = -I$(CPLEXDIR)/include -I$(CONCERTDIR)/include
CC = g++
CCFLAGS = -O3 -std=c++17 -DIL_STD -I$(PYTHON_INCLUDE)

# ---------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------
SRC = main_modele1_lineaire.cpp modele1_lineaire.cpp
OBJ = $(SRC:.cpp=.o)
TARGET = mcnd_mip

# ---------------------------------------------------------------------
# Règles
# ---------------------------------------------------------------------
all: $(TARGET)

$(TARGET): $(OBJ)
	$(CC) $(CCFLAGS) -o $@ $(OBJ) \
	-L$(CPLEXLIBDIR) -lilocplex -lcplex \
	-L$(CONCERTLIBDIR) -lconcert \
	$(PYTHON_EMBED_LDFLAGS) # Nouvelle ligne utilisant python3.11-config
	
%.o: %.cpp
	$(CC) $(CCFLAGS) $(CPLEXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJ) $(TARGET) branching_data.csv ml_model.pkl regression_weights.txt fileout

# Crée l'environnement virtuel s'il n'existe pas
venv:
	test -d $(VENV) || $(PYTHON) -m venv $(VENV)
	$(VENV)/bin/pip install --upgrade pip
	$(VENV)/bin/pip install pandas scikit-learn joblib numpy

# Entraîne le modèle ML dans l’environnement virtuel
train: venv
	$(VENV)/bin/python entraine.py

.PHONY: all clean train venv
