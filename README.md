# SiWEcal 

Geant4 simulation of a CALICE Silicon-Tungsten 3 layer ECAL 

More info: 

1. https://twiki.cern.ch/twiki/bin/view/Main/SiWECALAnalysis

2. https://agenda.linearcollider.org/event/6557/session/0/contribution/164

Geometry implementation is instantiated by detector versions in an enum: src/DetectorConstruction.cc and src/SamplingSection.cc

A small ntuple is stored with the energy deposits: src/EventAction.cc 

When changing the ttree content, adding homemade classes rather than
simple objects, the classes should be added to userlib. The dictionary
for root to understand the classes also need to be remade. Use "make
dictionary" before make, inside of userlib/. 

## Setup the environment (SLC6)

source g4env.sh

## Compile

mkdir -p userlib/{lib,obj,bin} && cd userlib && make dictionary && make -j 5 && cd - && make -j 5

## Run and visualization

For visualization go to the executable folder (in lxplus $HOME/geant4_workdir/bin/Linux-g++/) and run with vis.mac as an argument 

SiWEcal vis.mac


