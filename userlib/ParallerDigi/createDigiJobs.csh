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
# Delete some auxilliary files just in case
rm voodoo voodoo1 voodoo2 voodoo3 voodoo4 voodoo5
# Clear eos space for new files
setenv eospath "/eos/cms/store/group/phys_b2g/apsallid/SiWEcal/Digi/DetectorConfigurations"
/afs/cern.ch/project/eos/installation/cms/bin/eos.select rm -r $eospath
# Starting the loop through all test beam setups, particles and energies
foreach detconf ($detectorconfigurations)
echo "===================================================================================="
echo "Detector configuration $detconf" 
foreach subdirineosguns  ($suddirineoslistguns)
foreach subdirineosenergies  ($suddirineoslistenergies)
echo "------------------------"
echo "${subdirineosguns} ${subdirineosenergies} GeV"

foreach num (`seq 1 1`)

#set num1=`expr ${num} \* 1257 + 158725983 `

cat setupDigi.sh > voodoo

setenv file SiWEcal_${detconf}_${subdirineosguns}_${subdirineosenergies}GeV_${num}.root  
setenv workpath "/afs/cern.ch/work/a/apsallid/CMS/Geant4/SiWEcal/userlib/ParallerDigi/DetectorConfigurations/Config_${detconf}/$subdirineosguns/$subdirineosenergies"

if ( $num == 1 ) then
rm -rf $workpath/jobs $workpath/logfiles $workpath/results
mkdir -p $workpath/jobs $workpath/logfiles $workpath/results
chmod 755 $workpath/jobs $workpath/logfiles $workpath/results

/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir $eospath 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir $eospath/Config_${detconf} 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir $eospath/Config_${detconf}/$subdirineosguns 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir $eospath/Config_${detconf}/$subdirineosguns/$subdirineosenergies 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select mkdir $eospath/Config_${detconf}/$subdirineosguns/$subdirineosenergies/results 

/afs/cern.ch/project/eos/installation/cms/bin/eos.select chmod 755 $eospath 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select chmod 755 $eospath/Config_${detconf} 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select chmod 755 $eospath/Config_${detconf}/$subdirineosguns 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select chmod 755 $eospath/Config_${detconf}/$subdirineosguns/$subdirineosenergies 
/afs/cern.ch/project/eos/installation/cms/bin/eos.select chmod 755 $eospath/Config_${detconf}/$subdirineosguns/$subdirineosenergies/results 
endif

sed -e "s/INPUT/$file/g" voodoo > voodoo1
sed -e "s/TYPEGUN/$subdirineosguns/g" voodoo1 > voodoo2
sed -e "s/TYPEENE/$subdirineosenergies/g" voodoo2 > voodoo3
sed -e "s/THEDETECTORCONFIG/${detconf}/g" voodoo3 > voodoo4

mv voodoo4 $workpath/jobs/digi_${num}.job
chmod 755 $workpath/jobs/digi_${num}.job

echo digi_${num}.job

rm voodoo 
rm voodoo1
rm voodoo2
rm voodoo3

end

end

end

end
~

