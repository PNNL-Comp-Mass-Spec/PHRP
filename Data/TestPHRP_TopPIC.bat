xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_005649_5ECA3D02.fasta TestData\FASTA\ /D
xcopy \\gigasax\DMS_FASTA_File_Archive\dynamic\forward\ID_006407_8F27399B.fasta TestData\FASTA\ /D

del TestData\TPC201808171745_Auto1624200_Compare\PHRP_LogFile.txt
..\bin\Debug\PeptideHitResultsProcRunner.exe TestData\TPC201808171745_Auto1624200\MZ20170525_Mnx_PFA_TopPIC_PrSMs.txt  /O:TestData\TPC201808171745_Auto1624200_Compare\ /N:TestData\TPC201808171745_Auto1624200\TopPIC_15ppmParTol_NumShift1_2018-08-16.txt /SynPvalue:0.2  /SynProb:0.05  /L:TestData\TPC201808171745_Auto1624200_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_005649_5ECA3D02.fasta

del TestData\TPC202102011550_Auto1869702_Compare\PHRP_LogFile.txt
..\bin\Debug\PeptideHitResultsProcRunner.exe TestData\TPC202102011550_Auto1869702\Hubmap_nanoPOTs_top_down_QC_intact_20ng_FAIMS_LowNCE_r3_TopPIC_PrSMs.txt  /O:TestData\TPC202102011550_Auto1869702_Compare\ /N:TestData\TPC202102011550_Auto1869702\TopPIC_15ppmParTol_NumShift1_2020-03-02.txt /SynPvalue:0.2  /SynProb:0.05  /L:TestData\TPC202102011550_Auto1869702_Compare\PHRP_LogFile.txt /ProteinMods /F:TestData\FASTA\ID_006407_8F27399B.fasta
