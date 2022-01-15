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

xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_007741_1516CCE3.fasta TestData\FASTA\ /D
del TestData\MSF202008051507_Auto1822588_Compare\PHRP_LogFile.txt
%ExePath% TestData\MSF202008051507_Auto1822588\MCF7Cell_O-GlcNAc_Petyuk_R1_15May20_Rage_Rep-20-05-01.tsv /AltBasePath:..\..\Data /O:TestData\MSF202008051507_Auto1822588_Compare\ /M:TestData\MSF202008051507_Auto1822588\MSFragger_Tryp_Dyn_MetOx_ProtNTermAcet_StatCysAlk_20ppmParTol_ModDefs.txt /T:TestData\MSF202008051507_Auto1822588\Mass_Correction_Tags.txt /N:TestData\MSF202008051507_Auto1822588\MSFragger_Tryp_Dyn_MetOx_ProtNTermAcet_StatCysAlk_20ppmParTol.params /SynPvalue:0.2 /SynProb:0.05 /L:TestData\MSF202008051507_Auto1822588_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_007741_1516CCE3_decoy.fasta /FHT:False /Syn:True

del TestData\MSF202008051507_Auto1822588_Compare_UsePSM\PHRP_LogFile.txt
%ExePath% TestData\MSF202008051507_Auto1822588\MCF7Cell_O-GlcNAc_Petyuk_R1_15May20_Rage_Rep-20-05-01_psm.tsv /AltBasePath:..\..\Data /O:TestData\MSF202008051507_Auto1822588_Compare_UsePSM\ /M:TestData\MSF202008051507_Auto1822588\MSFragger_Tryp_Dyn_MetOx_ProtNTermAcet_StatCysAlk_20ppmParTol_ModDefs.txt /T:TestData\MSF202008051507_Auto1822588\Mass_Correction_Tags.txt /N:TestData\MSF202008051507_Auto1822588\MSFragger_Tryp_Dyn_MetOx_ProtNTermAcet_StatCysAlk_20ppmParTol.params /SynPvalue:0.2 /SynProb:0.05 /L:TestData\MSF202008051507_Auto1822588_Compare_UsePSM\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_007741_1516CCE3_decoy.fasta /FHT:False /Syn:True

xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_008026_7A1842EC.fasta TestData\FASTA\ /D
del TestData\MSF202110151330_Auto1966198_Compare\PHRP_LogFile.txt
%ExePath% TestData\MSF202110151330_Auto1966198\Aggregation_psm.tsv /AltBasePath:..\..\Data /O:TestData\MSF202110151330_Auto1966198_Compare\ /T:TestData\MSF202110151330_Auto1966198\Mass_Correction_Tags.txt /N:TestData\MSF202110151330_Auto1966198\MSFragger_Tryp_Dyn_MetOx_ProtNTermAcet_StatCysAlk_20ppmParTol.params /SynPvalue:0.2 /SynProb:0.05 /L:TestData\MSF202110151330_Auto1966198_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_008026_7A1842EC_decoy.fasta

:Done
