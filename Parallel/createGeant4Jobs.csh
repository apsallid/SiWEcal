#!/bin/tcsh

# Detector Configurations
# For details: https://twiki.cern.ch/twiki/bin/view/Main/SiWECALAnalysisDataStructure
setenv detectorconfigurations "1 2 3 4 5 6 7 8 9 10 11"
# Particle 
setenv suddirineoslistguns "e+ pi+ mu+"
# Energies in GeV
setenv suddirineoslistenergies "15 30 50 80 100 150 "
# Delete some auxilliary files just in case
rm voodoo voodoo1 voodoo2 voodoo3 voodoo4 voodoo5
# Starting the loop through all test beam setups, particles and energies
foreach detconf ($detectorconfigurations)
echo "===================================================================================="
echo "Detector configuration $detconf" 
foreach subdirineosguns  ($suddirineoslistguns)
foreach subdirineosenergies  ($suddirineoslistenergies)
echo "------------------------"
echo "${subdirineosguns} ${subdirineosenergies} GeV"

foreach num (`seq 1 10`)

set num1=`expr ${num} \* 1257 + 158725983 `

cat setupGeant4.sh > voodoo

setenv file SiWEcal_${detconf}_${subdirineosguns}_${subdirineosenergies}GeV_${num}.root  
setenv workpath "/afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/Parallel/DetectorConfigurations/Config_${detconf}/$subdirineosguns/$subdirineosenergies"

if ($num == 1) then
rm -rf $workpath/jobs $workpath/logfiles $workpath/results
mkdir -p $workpath/jobs $workpath/logfiles $workpath/results
chmod 755 $workpath/jobs $workpath/logfiles $workpath/results
endif

sed -e "s/OUTPUT/$file/g" voodoo > voodoo1
sed -e "s/TYPEGUN/$subdirineosguns/g" voodoo1 > voodoo2
sed -e "s/TYPEENE/$subdirineosenergies/g" voodoo2 > voodoo3
sed -e "s/THESEEDNUM/${num1}/g" voodoo3 > voodoo4
sed -e "s/THEDETECTORCONFIG/${detconf}/g" voodoo4 > voodoo5

mv voodoo5 $workpath/jobs/geant_${num}.job
chmod 755 $workpath/jobs/geant_${num}.job

echo geant_${num}.job

rm voodoo 
rm voodoo1
rm voodoo2
rm voodoo3
rm voodoo4

end

end

end

end
~

