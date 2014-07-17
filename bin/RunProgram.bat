@echo off
rem PeptideHitResultsProcRunner.exe /i:EIF_NAF_10NPI_1_15Jan06_Doc_05-11-10_xt.xml /m:xtandem_Rnd1Tryp_Rnd2DynMetOxCysIodoNTermAcet_ModDefs.txt
rem PeptideHitResultsProcRunner.exe /i:EIF_NAF_11CPI_1_14Jan06_Doc_05-11-11_syn.txt /m:sequest_N14_NE_Stat_C_Iodoacetimide_ModDefs.txt

rem PeptideHitResultsProcRunner.exe /i:QC_Shew_07_02_pt5_27Aug07_Sphinx_07-04-11_inspect.txt /m:Inspect_NoMods_FTHybrid_ModDefs.txt /n:Inspect_NoMods_FTHybrid.txt /t:Mass_Correction_Tags.txt


rem PeptideHitResultsProcRunner.exe /i:INS200901081452_Auto358179\QC_Shew_07_02_pt5_27Aug07_Sphinx_07-04-11_inspect.txt /m:INS200901081452_Auto358179\Inspect_CysAlk_ModDefs.txt /n:INS200901081452_Auto358179\Inspect_CysAlk.txt /t:INS200901081452_Auto358179\Mass_Correction_Tags.txt /L

rem PeptideHitResultsProcRunner.exe /i:INS200902131714_Auto365855\Bats_0012g_OrbiA_19Dec07_Draco_07-09-13_inspect.txt /m:INS200902131714_Auto365855\Inspect_MetOx_FTHybrid_50ppm_ModDefs.txt /n:INS200902131714_Auto365855\Inspect_MetOx_FTHybrid_50ppm.txt /t:INS200902131714_Auto365855\Mass_Correction_Tags.txt /L

PeptideHitResultsProcRunner.exe /i:MSGFPlus_Example\QC_Shew_13_05b_HCD_500ng_24Mar14_Tiger_14-03-04_msgfplus.tsv /m:MSGFPlus_Example\MSGFDB_PartTryp_MetOx_20ppmParTol_ModDefs.txt /n:MSGFPlus_Example\MSGFDB_PartTryp_MetOx_20ppmParTol.txt /t:MSGFPlus_Example\Mass_Correction_Tags.txt /L /ProteinMods /F:MSGFPlus_Example\Shewanella_oneidensis_MR1_2010-04-22_Tryp_Pig_Bov.revcat.fasta

rem cd x64
rem PeptideHitResultsProcRunner.exe /i:..\MSGFPlus_Example\QC_Shew_13_05b_HCD_500ng_24Mar14_Tiger_14-03-04_msgfplus.tsv /m:..\MSGFPlus_Example\MSGFDB_PartTryp_MetOx_20ppmParTol_ModDefs.txt /n:..\MSGFPlus_Example\MSGFDB_PartTryp_MetOx_20ppmParTol.txt /t:..\MSGFPlus_Example\Mass_Correction_Tags.txt /L /ProteinMods /F:..\MSGFPlus_Example\Shewanella_oneidensis_MR1_2010-04-22_Tryp_Pig_Bov.revcat.fasta


pause
