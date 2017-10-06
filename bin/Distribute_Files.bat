copy Release\PeptideHitResultsProcessor.dll "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common"
copy Release\PeptideHitResultsProcessor.dll "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin"

rem PeptideHitResultsProcRunner is compiled as AnyCPU
copy Release\*.dll "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP"
copy Release\PeptideHitResultsProcessor.pdb "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP"
copy Release\PeptideHitResultsProcRunner.exe "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP"
copy Release\PeptideHitResultsProcRunner.pdb "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP"

copy Release\*.dll "C:\DMS_Programs\PHRP"
copy Release\PHRPReader.pdb "C:\DMS_Programs\PHRP"
copy Release\PeptideHitResultsProcessor.pdb "C:\DMS_Programs\PHRP"
copy Release\PeptideHitResultsProcRunner.exe "C:\DMS_Programs\PHRP"
copy Release\PeptideHitResultsProcRunner.pdb "C:\DMS_Programs\PHRP"

pause

if not exist "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86" mkdir "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86"
if not exist "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64" mkdir "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64"
xcopy "Release\x86\SQLite.Interop.dll" "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86\" /D /Y
xcopy "Release\x64\SQLite.Interop.dll" "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64\" /D /Y

if not exist "C:\DMS_Programs\PHRP\x86" mkdir "C:\DMS_Programs\PHRP\x86"
if not exist "C:\DMS_Programs\PHRP\x64" mkdir "C:\DMS_Programs\PHRP\x64"
xcopy "Release\x86\SQLite.Interop.dll" "C:\DMS_Programs\PHRP\x86\" /D /Y
xcopy "Release\x64\SQLite.Interop.dll" "C:\DMS_Programs\PHRP\x64\" /D /Y

if exist "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\SQLite.Interop.dll" (del "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\SQLite.Interop.dll")
if exist "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin\SQLite.Interop.dll" (del "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin\SQLite.Interop.dll")
if exist "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\SQLite.Interop.dll" (del "F:\My Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\SQLite.Interop.dll")
if exist "C:\DMS_Programs\PHRP\SQLite.Interop.dll" (del "C:\DMS_Programs\PHRP\SQLite.Interop.dll")

pause
