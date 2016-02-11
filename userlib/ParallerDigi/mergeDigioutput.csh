#!/bin/tcsh

# Detector Configurations
# For details: https://twiki.cern.ch/twiki/bin/view/Main/SiWECALAnalysisDataStructure
setenv detectorconfigurations "1 2 3 4 5 6 7 8 9 10 11"
# Particle 
#setenv suddirineoslistguns "e+ pi+ mu+"
setenv suddirineoslistguns "mu- mu+"
# Energies in GeV
#setenv suddirineoslistenergies "15 30 50 80 100 150 "
#setenv suddirineoslistenergies "15 30 50"
setenv suddirineoslistenergies "15"
setenv eospath "/eos/cms/store/group/phys_b2g/apsallid/SiWEcal/Digi/DetectorConfigurations"
# Starting the loop through all test beam setups, particles and energies
foreach detconf ($detectorconfigurations)
echo "===================================================================================="
echo "Detector configuration $detconf"
foreach subdirineosguns  ($suddirineoslistguns)
echo "===================================================================================="
foreach subdirineosenergies  ($suddirineoslistenergies)
echo "------------------------"
echo "${subdirineosguns} ${subdirineosenergies} GeV"

#setenv workpath "/afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/Parallel/DetectorConfigurations/Config_${detconf}/$subdirineosguns/$subdirineosenergies"
setenv workpath "$eospath/Config_${detconf}/$subdirineosguns/$subdirineosenergies"

setenv output SiWEcal_${detconf}_${subdirineosguns}_${subdirineosenergies}.root

#mkdir -p /tmp/apsallid/Config_${detconf}/$subdirineosguns/$subdirineosenergies

#/afs/cern.ch/project/eos/installation/cms/bin/eos.select cp -r ${workpath}/ /tmp/apsallid/Config_${detconf}/$subdirineosguns/$subdirineosenergies/

/afs/cern.ch/project/eos/installation/cms/bin/eos.select cp -r ${workpath}/ /tmp/apsallid/Configs/


#hadd $output ${workpath}/SiWEcal_${detconf}_${subdirineosguns}_${subdirineosenergies}GeV_*.root
#echo "/tmp/apsallid/Config_${detconf}/$subdirineosguns/$subdirineosenergies/SiWEcal_${detconf}_${subdirineosguns}_${subdirineosenergies}GeV_*.root"
 
#hadd -f /tmp/apsallid/$output /tmp/apsallid/Config_${detconf}/$subdirineosguns/$subdirineosenergies/results/SiWEcal_${detconf}_${subdirineosguns}_${subdirineosenergies}GeV_*.root

end 

end

end

