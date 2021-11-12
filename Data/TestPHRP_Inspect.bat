xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_002614_23305E80.fasta TestData\FASTA\ /D
del TestData\INS201108111004_Auto731708_Compare\PHRP_LogFile.txt
..\bin\Debug\PeptideHitResultsProcRunner.exe TestData\INS201108111004_Auto731708\QC_Shew_11_03_0pt5_a_10Aug11_Cougar_11-01-17_inspect.txt /O:TestData\INS201108111004_Auto731708_Compare\ /M:TestData\INS201108111004_Auto731708\Inspect_MetOx_FTHybrid_20ppm_ModDefs.txt /N:TestData\INS201108111004_Auto731708\Inspect_MetOx_FTHybrid_20ppm.txt /SynPvalue:0.2  /SynProb:0.05  /L:TestData\INS201108111004_Auto731708_Compare\PHRP_LogFile.txt /F:TestData\FASTA\ID_002614_23305E80.fasta

