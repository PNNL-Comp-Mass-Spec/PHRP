@echo off
echo Option 1: use the Mzid-To-Tsv-Converter standalone application, available at https://github.com/PNNL-Comp-Mass-Spec/Mzid-To-Tsv-Converter/releases

echo Option 2: use MSGFPlus.jar

@echo on
java.exe -Xmx1000M -cp C:\DMS_Programs\MSGFPlus\MSGFPlus.jar edu.ucsd.msjava.ui.MzIDToTsv -i MSGFPlus_Job1217459\48488_Aniger_BU_F9_C_8Jun15_Tiger_15-03-05_msgfplus.mzid -o MSGFPlus_Job1217459\48488_Aniger_BU_F9_C_8Jun15_Tiger_15-03-05_msgfplus.tsv -showQValue 1 -showDecoy 1 -unroll 1

pause
