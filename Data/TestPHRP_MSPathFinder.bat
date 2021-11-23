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

xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_006407_8F27399B.fasta TestData\FASTA\ /D
del TestData\MSP201804280928_Auto1578533_Compare\PHRP_LogFile.txt
%ExePath% TestData\MSP201804280928_Auto1578533\QC_Shew_17_01_4_27Apr18_Merry_18-02-04_IcTda.tsv  /O:TestData\MSP201804280928_Auto1578533_Compare\ /M:TestData\MSP201804280928_Auto1578533\MSPF_MetOx_CysDehydro_NTermAcet_SingleInternalCleavage_ModDefs.txt /N:TestData\MSP201804280928_Auto1578533\MSPF_MetOx_CysDehydro_NTermAcet_SingleInternalCleavage.txt /SynPvalue:0.2  /SynProb:0.05  /L:TestData\MSP201804280928_Auto1578533_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_006407_8F27399B.fasta

:Done
