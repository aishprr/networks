
Large Scale Biological Networks Attack:

By Ryan Park, Aishwarya Prem Renu

This simulation simulates the following organisms and cellular structures

1. Bacterial Cells (Upto 5 different species)
2. Auto Inducers associated with each species
3. Macrophages
4. Helper T-cells
5. Killer T-cells

High level working of the simulation is with every time step, certain operations
take place. A cycle repeats every STEP_MULTIPLE time steps. 
Some of time steps out of STEP_MULTIPLE time steps is dedicated to reproduction
of bacterial cells, production of autoinducers, and the production of killer 
T-cells. Some of the time steps are dedicated to movement of all of the cells 
based on biological contraints. The movement of bacteria is based on the 
concentration of auto-inducers in the immediate environment, and the concentration
of autoinducers is in turn dependent on the concentration of the bacteria since
they are produced by the bacterial cells. Each bacterial cell produces AI_PER_BAC
number of AI molecules in every time step.
When the bacterial cells are closer than BACT_DIS_THRESH from each other, they
make connections with one another if both the cells involved have less than 
BACT_DEG_THRESH degree.
When macrophages come closeer to bacterial cells than MACRO_EAT_DIS, if they 
have less than MACRO_MAX_BACT_EAT that they have eaten, but not ingested yet, 
they will eat the bactera they are close to up to that limit.
After the bacteria inside the macrophage mature up BACT_IN_MACRO_REPRO_AGE number
of time steps, if they are still surviving they product BACT_IN_MACRO_REPR 
children each. Once there are more than MACRO_BACT_TO_DIE inside the macrophage, 
the bacteria take over the macrophage, and kill it from within.
Once macrophages with bacterial cells inside come closer than HELPER_MACRO_DIS
to helper T-cells, this will cause the production of KILL_PER_MAC number of 
killer T-cells. Once these killer T-cells come close than KILLTCELL_EAT_DIS
to a bacterial cell, the bacterial cell gets eaten by the kill T-cell, and each
killer T-cell can only eat one bacteria per time step. 
Strength of each bacterial species starts off at the BACT_STRENGTH vector 
initially, and this strength vector changes with time based on the connectedness
of the previous species in the cascade, and the cascade order is in increasing
order of BactType number. There exist BACT_IN_MACRO_REPRO_THRESH and 
BACT_IN_MACRO_KILL_THRESH which determine the thresholds of strength for the 
bacterial species to have special powers to reproduce inside the macrophage and
kill the macrophage respectively.
Initial bacterial count is determined by BACT_INIT_COUNT, and their counts are
limited by BACT_COUNT_LIMIT. Similar init counts exist for all the cells.

Speeds of movement of all the cells are determined by the speed variables for 
all the cells.
AIC_HCOUNT and AIC_WCOUNT determine the number of grid boxes along the vertical
axis and the horizontal axis respectively. This determines the granularity of 
concentration calculation.



Invocation:

python bact_mov.py

Files: 

bact_mov.py
Can also be found at: https://github.com/aishprr/networks

Libraries required:

1. networkx
2. matplotlib
3. numpy

