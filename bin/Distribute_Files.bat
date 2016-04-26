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

if not exist "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86" mkdir "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86"
if not exist "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64" mkdir "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64"
xcopy "x86\SQLite.Interop.dll" "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86\" /D /Y
xcopy "x64\SQLite.Interop.dll" "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64\" /D /Y

if not exist "C:\DMS_Programs\PHRP\x86" mkdir "C:\DMS_Programs\PHRP\x86"
if not exist "C:\DMS_Programs\PHRP\x64" mkdir "C:\DMS_Programs\PHRP\x64"
xcopy "x86\SQLite.Interop.dll" "C:\DMS_Programs\PHRP\x86\" /D /Y
xcopy "x64\SQLite.Interop.dll" "C:\DMS_Programs\PHRP\x64\" /D /Y

if exist "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\SQLite.Interop.dll" (del "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\SQLite.Interop.dll")
if exist "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin\SQLite.Interop.dll" (del "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin\SQLite.Interop.dll")
if exist "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\SQLite.Interop.dll" (del "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\SQLite.Interop.dll")
if exist "C:\DMS_Programs\PHRP\SQLite.Interop.dll" (del "C:\DMS_Programs\PHRP\SQLite.Interop.dll")

pause
