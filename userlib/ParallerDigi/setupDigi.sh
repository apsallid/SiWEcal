#!/bin/bash

echo "The script starts now."

echo "System: "
uname -a

source /afs/cern.ch/cms/LCG/LCG-2/UI/cms_ui_env.sh

source /afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/g4env.sh

export EOSMAINPATH=/eos/cms/store/group/phys_b2g/apsallid/SiWEcal/

export PWD=`pwd`

cp /afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/userlib/bin/digitizer $PWD/. 

$PWD/digitizer 0 root://eoscms/$EOSMAINPATH/DetectorConfigurations/Config_THEDETECTORCONFIG/TYPEGUN/TYPEENE/results/INPUT $PWD '0-11:1' '0-11:0.05556' '0-11:0' 1 12 0 4 1 1 

mv DigiPFcal_withDigiHits_withSimHits.root Digi_INPUT

#cp OUTPUT /afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/Parallel/DetectorConfigurations/Config_THEDETECTORCONFIG/TYPEGUN/TYPEENE/results 

/afs/cern.ch/project/eos/installation/cms/bin/eos.select cp $PWD/Digi_INPUT $EOSMAINPATH/Digi/DetectorConfigurations/Config_THEDETECTORCONFIG/TYPEGUN/TYPEENE/results/Digi_INPUT


 
