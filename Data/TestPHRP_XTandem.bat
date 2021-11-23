
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
del TestData\XTM201305070949_Auto940873_Compare\PHRP_LogFile.txt
%ExePath% TestData\XTM201305070949_Auto940873\*_xt.xml  /O:TestData\XTM201305070949_Auto940873_Compare\ /M:TestData\XTM201305070949_Auto940873\xtandem_Rnd1PartTryp_Rnd2DynMetOx_ModDefs.txt /T:TestData\XTM201305070949_Auto940873\Mass_Correction_Tags.txt /N:TestData\XTM201305070949_Auto940873\xtandem_Rnd1PartTryp_Rnd2DynMetOx.xml /SynPvalue:0.2  /SynProb:0.05  /L:TestData\XTM201305070949_Auto940873_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_003456_9B916A8B.fasta

xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_004208_295531A4.fasta TestData\FASTA\ /D
del TestData\XTM201702142141_Auto1408699_Compare\PHRP_LogFile.txt
%ExePath% TestData\XTM201702142141_Auto1408699\*_xt.xml  /O:TestData\XTM201702142141_Auto1408699_Compare\ /M:TestData\XTM201702142141_Auto1408699\xtandem_ETD_Rnd1Tryp_StatCysAlk_STYPhos_NoRefinement_20ppmParent_0pt5DaFrag_ModDefs.txt /T:TestData\XTM201702142141_Auto1408699\Mass_Correction_Tags.txt /N:TestData\XTM201702142141_Auto1408699\xtandem_ETD_Rnd1Tryp_StatCysAlk_STYPhos_NoRefinement_20ppmParent_0pt5DaFrag.xml /SynPvalue:0.2  /SynProb:0.05  /L:TestData\XTM201702142141_Auto1408699_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_004208_295531A4.fasta

:Done
