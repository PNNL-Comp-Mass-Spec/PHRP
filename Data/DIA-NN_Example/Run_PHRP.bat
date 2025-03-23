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

%ExePath% /i:QC_Hela_23_01_DIA_15min_5mz_10IT_2ng_Monty_PepMap_75x50_ES903_Trap_r3_report.parquet /m:DiaNN_Tryp_Dyn_MetOx_NTermAcet_Stat_CysAlk_Precursor400-900_ModDefs.txt /n:DiaNN_Tryp_Dyn_MetOx_NTermAcet_Stat_CysAlk_Precursor400-900.txt /t:Mass_Correction_Tags.txt /L /ProteinMods /F:H_sapiens_UniProt_SPROT_2023-09-01_Tryp_Pig_Bov.fasta > PHRP_ConsoleOutput.txt

:Done

pause
