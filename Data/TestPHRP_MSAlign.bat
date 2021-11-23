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

xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_003679_A3EA5E14.fasta TestData\FASTA\ /D
del TestData\MSA201211061153_Auto892640_Compare\PHRP_LogFile.txt
%ExePath% TestData\MSA201211061153_Auto892640\FAB_reduced_2_CID_MSAlign_ResultTable.txt  /O:TestData\MSA201211061153_Auto892640_Compare\ /M:TestData\MSA201211061153_Auto892640\MSAlign_CID_15ppm_2011-10-13_ModDefs.txt /T:TestData\MSA201211061153_Auto892640\Mass_Correction_Tags.txt /N:TestData\MSA201211061153_Auto892640\MSAlign_CID_15ppm_2011-10-13.txt /SynPvalue:0.2  /SynProb:0.05  /L:TestData\MSA201211061153_Auto892640_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_003679_A3EA5E14.fasta

xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_007359_B85F4A48.fasta TestData\FASTA\ /D
del TestData\MSA201807121159_Auto1609998_Compare\PHRP_LogFile.txt
%ExePath% TestData\MSA201807121159_Auto1609998\20180709PDX_13_MSAlign_ResultTable.txt  /O:TestData\MSA201807121159_Auto1609998_Compare\ /T:TestData\MSA201211061153_Auto892640\Mass_Correction_Tags.txt /N:TestData\MSA201807121159_Auto1609998\MSAlign_15ppm_0pt01_FDR_2012-01-03.txt /SynPvalue:0.2  /SynProb:0.05  /L:TestData\MSA201807121159_Auto1609998_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_007359_B85F4A48.fasta

:Done
