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

xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_003456_9B916A8B.fasta TestData\FASTA\ /D
del TestData\MXQ202103181341_Auto1878805_Compare\PHRP_LogFile.txt
%ExePath% TestData\MXQ202103181341_Auto1878805\msms.txt /AltBasePath:..\..\Data /O:TestData\MXQ202103181341_Auto1878805_Compare\ /M:TestData\MXQ202103181341_Auto1878805\MaxQuant_Tryp_Dyn_MetOx_NTermAcet_20ppmParTol_ModDefs.txt /T:TestData\MXQ202103181341_Auto1878805\Mass_Correction_Tags.txt /N:TestData\MXQ202103181341_Auto1878805\MaxQuant_Tryp_Dyn_MetOx_NTermAcet_20ppmParTol.xml /SynPvalue:0.2  /SynProb:0.05  /L:TestData\MXQ202103181341_Auto1878805_Compare\PHRP_LogFile.txt /ProteinMods  /F:TestData\FASTA\ID_003456_9B916A8B.fasta /FHT:False /Syn:True

xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_003456_9B916A8B.fasta TestData\FASTA\ /D
del TestData\MXQ202105151319_Auto1899827_Compare\PHRP_LogFile.txt
%ExePath% TestData\MXQ202105151319_Auto1899827\txt\msms.txt /AltBasePath:..\..\Data /O:TestData\MXQ202105151319_Auto1899827_Compare\ /T:TestData\MXQ202105151319_Auto1899827\Mass_Correction_Tags.txt /N:TestData\MXQ202105151319_Auto1899827\MaxQuant_Tryp_Dyn_MetOx_NTermAcet_Stat_CysAlk_TMT_6Plex_20ppmParTol.xml /SynPvalue:0.2  /SynProb:0.05  /L:TestData\MXQ202105151319_Auto1899827_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_003456_9B916A8B.fasta

:Done
