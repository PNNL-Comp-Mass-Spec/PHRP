@echo off

echo.
echo Be sure to build PeptideHitResultsProcessor.dll in Release mode prior to calling this batch file
echo.
pause

@echo on

xcopy Release\PeptideHitResultsProcessor.dll "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common" /D /Y
xcopy Release\PeptideHitResultsProcessor.dll "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin" /D /Y

rem PeptideHitResultsProcRunner is compiled as AnyCPU
xcopy Release\*.dll "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y
xcopy Release\PeptideHitResultsProcessor.pdb "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y
xcopy Release\PeptideHitResultsProcRunner.exe "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y
xcopy Release\PeptideHitResultsProcRunner.pdb "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y
xcopy Release\PHRPReader.pdb "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y
xcopy Release\PHRPReader.xml "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y

xcopy Release\*.dll "C:\DMS_Programs\PHRP" /D /Y
xcopy Release\PHRPReader.pdb "C:\DMS_Programs\PHRP" /D /Y
xcopy Release\PeptideHitResultsProcessor.pdb "C:\DMS_Programs\PHRP" /D /Y
xcopy Release\PeptideHitResultsProcRunner.exe "C:\DMS_Programs\PHRP" /D /Y
xcopy Release\PeptideHitResultsProcRunner.pdb "C:\DMS_Programs\PHRP" /D /Y
xcopy Release\PHRPReader.pdb "C:\DMS_Programs\PHRP" /D /Y
xcopy Release\PHRPReader.xml "C:\DMS_Programs\PHRP" /D /Y

xcopy Release\*.dll \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Release\PHRPReader.pdb \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Release\PeptideHitResultsProcessor.pdb \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Release\PeptideHitResultsProcRunner.exe \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Release\PeptideHitResultsProcRunner.pdb \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Release\PHRPReader.pdb \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Release\PHRPReader.xml \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y

pause

if not exist "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86" mkdir "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86"
if not exist "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64" mkdir "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64"
xcopy "Release\x86\SQLite.Interop.dll" "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86\" /D /Y
xcopy "Release\x64\SQLite.Interop.dll" "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64\" /D /Y

if not exist "C:\DMS_Programs\PHRP\x86" mkdir "C:\DMS_Programs\PHRP\x86"
if not exist "C:\DMS_Programs\PHRP\x64" mkdir "C:\DMS_Programs\PHRP\x64"
xcopy "Release\x86\SQLite.Interop.dll" "C:\DMS_Programs\PHRP\x86\" /D /Y
xcopy "Release\x64\SQLite.Interop.dll" "C:\DMS_Programs\PHRP\x64\" /D /Y

if exist "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\SQLite.Interop.dll" (del "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\SQLite.Interop.dll")
if exist "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin\SQLite.Interop.dll" (del "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin\SQLite.Interop.dll")
if exist "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\SQLite.Interop.dll" (del "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\SQLite.Interop.dll")
if exist "C:\DMS_Programs\PHRP\SQLite.Interop.dll" (del "C:\DMS_Programs\PHRP\SQLite.Interop.dll")

pause
