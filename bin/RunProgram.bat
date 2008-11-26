@echo off
rem PeptideHitResultsProcRunner.exe /i:EIF_NAF_10NPI_1_15Jan06_Doc_05-11-10_xt.xml /m:xtandem_Rnd1Tryp_Rnd2DynMetOxCysIodoNTermAcet_ModDefs.txt
rem PeptideHitResultsProcRunner.exe /i:EIF_NAF_11CPI_1_14Jan06_Doc_05-11-11_syn.txt /m:sequest_N14_NE_Stat_C_Iodoacetimide_ModDefs.txt

PeptideHitResultsProcRunner.exe /i:QC_Shew_07_02_pt5_27Aug07_Sphinx_07-04-11_inspect.txt /m:Inspect_NoMods_FTHybrid_ModDefs.txt /n:Inspect_NoMods_FTHybrid.txt /t:Mass_Correction_Tags.txt
