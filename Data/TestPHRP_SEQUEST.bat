xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\decoy\ID_003797_2DFD81D6.fasta TestData\FASTA\ /D
..\bin\Debug\PeptideHitResultsProcRunner.exe TestData\Seq201211261855_Auto896700\Sor_Lysate_HR_30Mar11_Draco_11-03-10_syn.txt  /O:TestData\Seq201211261855_Auto896700_Compare\ /M:TestData\Seq201211261855_Auto896700\sequest_N14_NE_StatC_Iodo_DynK_Guanid_ModDefs.txt /T:TestData\Seq201211261855_Auto896700\Mass_Correction_Tags.txt /N:TestData\Seq201211261855_Auto896700\sequest_N14_NE_StatC_Iodo_DynK_Guanid.params /SynPvalue:0.2  /SynProb:0.05  /L:PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_003797_2DFD81D6.fasta
pause

rem Note: This FASTA file has nucleic acids
xcopy \\gigasax\DMS_Organism_Files\None\TestData\FASTA\NA2B2_Ds_EST.fasta TestData\FASTA\ /D
..\bin\Debug\PeptideHitResultsProcRunner.exe TestData\Seq201212131618_Auto901984\121130_Polle_6_agi_150_syn.txt  /O:TestData\Seq201212131618_Auto901984_Compare\ /M:TestData\Seq201212131618_Auto901984\sequest_DNA_N14_NE_Dyn_Met_Ox_Stat_Cys_Iodo_ModDefs.txt /T:TestData\Seq201212131618_Auto901984\Mass_Correction_Tags.txt /N:TestData\Seq201212131618_Auto901984\sequest_DNA_N14_NE_Dyn_Met_Ox_Stat_Cys_Iodo.params /SynPvalue:0.2  /SynProb:0.05  /L:PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\NA2B2_Ds_EST.fasta
