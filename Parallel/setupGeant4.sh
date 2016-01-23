#!/bin/bash

echo "The script starts now."

echo "System: "
uname -a

source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh

source /afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/g4env.sh

export PWD=`pwd`

cp /afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/runparallel.mac .

sed -e "s/PARTICLEGUN/TYPEGUN/g" runparallel.mac > voodoo
sed -e "s/ENERGYOFGUN/TYPEENE/g" voodoo > voodoo1
sed -e "s/SEEDNUM/THESEEDNUM/g" voodoo1 > voodoo2

mv voodoo2 run.mac 

SiWEcal run.mac THEDETECTORCONFIG 

mv SiWEcal.root OUTPUT  

cp OUTPUT /afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/Parallel/DetectorConfigurations/Config_THEDETECTORCONFIG/TYPEGUN/TYPEENE/results 

#/afs/cern.ch/project/eos/installation/cms/bin/eos.select cp $PWD/OUTPUT $EOSMAINPATH/TYPE/PandoraOutput/OUTPUT 


 
