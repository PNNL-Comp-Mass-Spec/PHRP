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
del TestData\MSG201205081554_Auto838496_Compare\PHRP_LogFile.txt
%ExePath% TestData\MSG201205081554_Auto838496\QC_Shew_11_06_Run-03_7May12_Roc_12-04-08_msgfdb.txt /AltBasePath:..\..\Data /O:TestData\MSG201205081554_Auto838496_Compare\ /M:TestData\MSG201205081554_Auto838496\MSGFPlus_PartTryp_MetOx_50ppmParTol_NoDecoy_EstimateFDR_ModDefs.txt /T:TestData\MSG201205081554_Auto838496\Mass_Correction_Tags.txt /N:TestData\MSG201205081554_Auto838496\MSGFPlus_PartTryp_MetOx_50ppmParTol_NoDecoy_EstimateFDR.txt /SynPvalue:0.2  /SynProb:0.05  /L:TestData\MSG201205081554_Auto838496_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_003456_9B916A8B.fasta

xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_004208_295531A4.fasta TestData\FASTA\ /D
del TestData\MSG201802011337_Auto1547784_Compare\PHRP_LogFile.txt
%ExePath% TestData\MSG201802011337_Auto1547784\MCF-7_pSTY_2_OldGr_15Mar17_Merry__IntEmtr_msgfplus.tsv /AltBasePath:..\..\Data /O:TestData\MSG201802011337_Auto1547784_Compare\ /M:TestData\MSG201802011337_Auto1547784\MSGFPlus_Tryp_DynSTYPhos_Stat_CysAlk_20ppmParTol_ModDefs.txt /T:TestData\MSG201802011337_Auto1547784\Mass_Correction_Tags.txt /N:TestData\MSG201802011337_Auto1547784\MSGFPlus_Tryp_DynSTYPhos_Stat_CysAlk_20ppmParTol.txt /SynPvalue:0.2  /SynProb:0.05  /L:TestData\MSG201802011337_Auto1547784_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_004208_295531A4_decoy.fasta

:Done
