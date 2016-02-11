#!/bin/tcsh

# Detector Configurations
# For details: https://twiki.cern.ch/twiki/bin/view/Main/SiWECALAnalysisDataStructure
setenv detectorconfigurations "1 2 3 4 5 6 7 8 9 10 11"
#setenv detectorconfigurations "2"
# Particle 
#setenv suddirineoslistguns "e+ pi+ mu+"
setenv suddirineoslistguns "mu- mu+"
# Energies in GeV
#setenv suddirineoslistenergies "15 30 50 80 100 150 "
setenv suddirineoslistenergies "15 30 50"
# Starting the loop through all test beam setups, particles and energies
foreach detconf ($detectorconfigurations)
echo "===================================================================================="
echo "Detector configuration $detconf"
foreach subdirineosguns  ($suddirineoslistguns)
echo "===================================================================================="
foreach subdirineosenergies  ($suddirineoslistenergies)
echo "------------------------"
echo "${subdirineosguns} ${subdirineosenergies} GeV"

setenv workpath "/afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/userlib/ParallerDigi/DetectorConfigurations/Config_${detconf}/$subdirineosguns/$subdirineosenergies/jobs"

setenv runumberslist ` ls -ltr ${workpath} | grep .job |  awk '{print $9}' `

foreach run  ($runumberslist)

#echo ${run}
chmod 755 ${workpath}/${run}

echo "Sending ${run}"
bsub -q 2nd -o /tmp/junk ${workpath}/${run}

end

end 

end

end



