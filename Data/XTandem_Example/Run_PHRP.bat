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

%ExePath% -i:QC_Shew_12_02_pt5_2c_20Dec12_Leopard_12-11-10_xt.xml -m:xtandem_Rnd1PartTryp_Rnd2DynMetOx_ModDefs.txt -n:xtandem_Rnd1PartTryp_Rnd2DynMetOx.xml -t:Mass_Correction_Tags.txt -L -ProteinMods -F:Shewanella_oneidensis_MR1_2010-04-22_Tryp_Pig_Bov.fasta > PHRP_ConsoleOutput.txt

:Done

pause
