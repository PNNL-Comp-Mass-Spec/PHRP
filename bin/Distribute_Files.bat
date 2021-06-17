@echo off

echo.
echo Build PeptideHitResultsProcessor.dll in Debug mode prior to calling this batch file
echo.
pause

@echo on

xcopy Debug\PeptideHitResultsProcessor.dll "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common" /D /Y
xcopy Debug\PeptideHitResultsProcessor.dll "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin" /D /Y

xcopy Debug\PHRPReader.* "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common" /D /Y
xcopy Debug\PHRPReader.* "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin" /D /Y

rem PeptideHitResultsProcRunner is compiled as AnyCPU
xcopy Debug\*.dll "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y
xcopy Debug\PeptideHitResultsProcessor.pdb "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y
xcopy Debug\PeptideHitResultsProcRunner.exe "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y
xcopy Debug\PeptideHitResultsProcRunner.pdb "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y
xcopy Debug\PHRPReader.pdb "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y
xcopy Debug\PHRPReader.xml "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP" /D /Y

xcopy Debug\*.dll "C:\DMS_Programs\PHRP" /D /Y
xcopy Debug\PHRPReader.pdb "C:\DMS_Programs\PHRP" /D /Y
xcopy Debug\PeptideHitResultsProcessor.pdb "C:\DMS_Programs\PHRP" /D /Y
xcopy Debug\PeptideHitResultsProcRunner.exe "C:\DMS_Programs\PHRP" /D /Y
xcopy Debug\PeptideHitResultsProcRunner.pdb "C:\DMS_Programs\PHRP" /D /Y
xcopy Debug\PHRPReader.pdb "C:\DMS_Programs\PHRP" /D /Y
xcopy Debug\PHRPReader.xml "C:\DMS_Programs\PHRP" /D /Y

if not exist "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86" mkdir "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86"
if not exist "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64" mkdir "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64"
xcopy Debug\x86\SQLite.Interop.dll "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x86\" /D /Y
xcopy Debug\x64\SQLite.Interop.dll "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\x64\" /D /Y

if not exist C:\DMS_Programs\PHRP\x86 mkdir C:\DMS_Programs\PHRP\x86
if not exist C:\DMS_Programs\PHRP\x64 mkdir C:\DMS_Programs\PHRP\x64
xcopy Debug\x86\SQLite.Interop.dll C:\DMS_Programs\PHRP\x86\ /D /Y
xcopy Debug\x64\SQLite.Interop.dll C:\DMS_Programs\PHRP\x64\ /D /Y

if exist "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\SQLite.Interop.dll" (del "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\SQLite.Interop.dll")
if exist "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin\SQLite.Interop.dll" (del "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Program\bin\SQLite.Interop.dll")
if exist "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\SQLite.Interop.dll" (del "F:\Documents\Projects\DataMining\DMS_Managers\Analysis_Manager\AM_Common\PHRP\SQLite.Interop.dll")
if exist "C:\DMS_Programs\PHRP\SQLite.Interop.dll" (del "C:\DMS_Programs\PHRP\SQLite.Interop.dll")

pushd ..\PHRPReader\bin
call Distribute_Files.bat NoPause
popd

@echo off
echo.
echo.
echo About to copy to \\pnl\projects\OmicsSW\DMS_Programs
echo.
pause
@echo on

xcopy Debug\*.dll \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Debug\PHRPReader.pdb \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Debug\PeptideHitResultsProcessor.pdb \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Debug\PeptideHitResultsProcRunner.exe \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Debug\PeptideHitResultsProcRunner.pdb \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Debug\PHRPReader.pdb \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y
xcopy Debug\PHRPReader.xml \\pnl\projects\OmicsSW\DMS_Programs\AnalysisToolManagerDistribution\PHRP\ /D /Y

pause
