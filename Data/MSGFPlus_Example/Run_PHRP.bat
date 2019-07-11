@echo off

set ExePath=PeptideHitResultsProcRunner.exe

if exist %ExePath% goto DoWork
if exist ..\%ExePath% set ExePath=..\%ExePath% && goto DoWork
if exist ..\..\Bin\Debug\%ExePath% set ExePath=..\..\Bin\Debug\%ExePath% && goto DoWork

echo Executable not found: %ExePath%
goto Done

:DoWork
echo.
echo Procesing with %ExePath%
echo.

%ExePath% /i:QC_Shew_13_05b_HCD_500ng_24Mar14_Tiger_14-03-04_msgfplus.tsv /m:MSGFDB_PartTryp_MetOx_20ppmParTol_ModDefs.txt /n:MSGFDB_PartTryp_MetOx_20ppmParTol.txt /t:Mass_Correction_Tags.txt /L /ProteinMods /F:Shewanella_oneidensis_MR1_2010-04-22_Tryp_Pig_Bov.fasta > PHRP_ConsoleOutput.txt

:Done

pause
