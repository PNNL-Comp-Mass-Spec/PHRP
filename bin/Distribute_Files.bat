copy PeptideHitResultsProcessor.dll "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common"
copy PeptideHitResultsProcessor.dll "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin"

rem PeptideHitResultsProcRunner is compiled as AnyCPU
copy *.dll "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP"
copy PeptideHitResultsProcessor.pdb "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP"
copy PeptideHitResultsProcRunner.exe "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP"
copy PeptideHitResultsProcRunner.pdb "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP"

copy *.dll "C:\DMS_Programs\PHRP"
copy PHRPReader.pdb "C:\DMS_Programs\PHRP"
copy PeptideHitResultsProcessor.pdb "C:\DMS_Programs\PHRP"
copy PeptideHitResultsProcRunner.exe "C:\DMS_Programs\PHRP"
copy PeptideHitResultsProcRunner.pdb "C:\DMS_Programs\PHRP"


pause
