@echo off

set ExePath=PeptideHitResultsProcRunner.exe

if exist %ExePath% goto DoWork
if exist ..\%ExePath% set ExePath=..\%ExePath% && goto DoWork
if exist ..\bin\Debug\%ExePath% set ExePath=..\bin\Debug\%ExePath% && goto DoWork

echo Executable not found: %ExePath%
goto Done

:DoWork
echo.
echo Processing with %ExePath%
echo.

xcopy \\gigasax\DMS_FASTA_File_Archive\Dynamic\Forward\ID_008358_7EC82878.fasta TestData\FASTA\ /D

del TestData\DNN202303271534_Auto2168355\PHRP_LogFile.txt
%ExePath% TestData\DNN202303271534_Auto2168355\report.tsv /AltBasePath:..\..\Data /O:TestData\DNN202303271534_Auto2168355_Compare\ /M:TestData\DNN202303271534_Auto2168355\DiaNN_Tryp_Dyn_MetOx_Stat_CysAlk_ModDefs.txt /N:TestData\DNN202303271534_Auto2168355\DiaNN_Tryp_Dyn_MetOx_Stat_CysAlk.txt  /DiaNNQValue:0.1 /DiaNNCScore:0.25 /L:TestData\DNN202303271534_Auto2168355_Compare\PHRP_LogFile.txt /F:TestData\FASTA\ID_008358_7EC82878.fasta

:Done
