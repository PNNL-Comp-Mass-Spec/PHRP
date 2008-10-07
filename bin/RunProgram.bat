@echo off
rem PeptideHitResultsProcRunner.exe /i:EIF_NAF_10NPI_1_15Jan06_Doc_05-11-10_xt.xml /m:xtandem_Rnd1Tryp_Rnd2DynMetOxCysIodoNTermAcet_ModDefs.txt
rem PeptideHitResultsProcRunner.exe /i:EIF_NAF_11CPI_1_14Jan06_Doc_05-11-11_syn.txt /m:sequest_N14_NE_Stat_C_Iodoacetimide_ModDefs.txt

PeptideHitResultsProcRunner.exe /i:Lowdose_control_IMAC_top10c_inspect.txt /m:Inspect_Dyn_STYPhosCtermMet_Stat_DEMet_FTHybrid_ModDefs.txt /n:Inspect_Dyn_STYPhosCtermMet_Stat_DEMet_FTHybrid.txt /t:Mass_Correction_Tags.txt
