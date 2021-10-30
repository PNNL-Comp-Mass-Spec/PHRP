@echo off

set ExePath=PeptideHitResultsProcRunner.exe

if exist %ExePath% goto DoWork
if exist ..\%ExePath% set ExePath=..\%ExePath% && goto DoWork
if exist ..\..\Bin\Debug\%ExePath% set ExePath=..\..\Bin\Debug\%ExePath% && goto DoWork

echo Executable not found: %ExePath%
goto Done

:DoWork
echo.
echo Processing with %ExePath%
echo.

%ExePath% /i:msms.txt /n:MaxQuant_Tryp_Stat_CysAlk_Dyn_MetOx_NTermAcet_20ppmParTol.xml /t:Mass_Correction_Tags.txt /L /ProteinMods /F:M_musculus_UniProt_SPROT_2013_09_2013-09-18_Tryp_Pig_Bov.fasta > PHRP_ConsoleOutput.txt

:Done

pause
