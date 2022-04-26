@echo off

if "%1"=="NoPause" Goto StartCopy
echo Be sure to build PHRPReader.dll in Debug mode
if not "%1"=="NoPause" pause

:StartCopy

@echo on
call Distribute_Files_Work.bat "..\..\..\CodeTestCS\lib"
call Distribute_Files_Work.bat "..\..\..\CodeTestCS\bin"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\AM_Common"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\AM_Common\PHRP"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\AM_Shared\bin\Debug"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\AM_Program\bin"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_Extraction_PlugIn\bin"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_IDPicker_PlugIn\bin"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_Mage_PlugIn\bin\Debug"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_MSAlign_Plugin\bin"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_MSGF_PlugIn\bin"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_MSGF_PlugIn\MSGFResultsSummarizerDLL\bin\Debug"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_MSGF_PlugIn\MSGFResultsSummarizerExe\bin\Debug"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_MSGFDB_PlugIn\bin"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_MSXML_Gen_PlugIn\bin"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_Phospho_FDR_Aggregator_PlugIn\bin"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_PRIDE_Converter_PlugIn\bin"
call Distribute_Files_Work.bat "..\..\..\DMS_Managers\Analysis_Manager\Plugins\AM_PRIDE_MzXML_PlugIn\bin"

call Distribute_Files_Work.bat "..\..\..\MASICResultsMerger\bin"
call Distribute_Files_Work.bat "..\..\..\MASICResultsMerger\Lib"
call Distribute_Files_Work.bat "..\..\..\MTDB_Creator\libs"
call Distribute_Files_Work.bat "..\..\..\PeptideHitResultsProcessor\bin\Debug"
call Distribute_Files_Work.bat "..\..\..\PeptideHitResultsProcessor\PeptideHitResultsProcessor\bin"
call Distribute_Files_Work.bat "..\..\..\PeptideHitResultsProcessor\CreateMSGFPlusResultsFileFromPHRP\bin\"
call Distribute_Files_Work.bat "..\..\..\PeptideHitResultsProcessor\Test_PHRPReader\bin"
call Distribute_Files_Work.bat "..\..\..\PeptideListToXML\bin"
call Distribute_Files_Work.bat "..\..\..\PeptideListToXML\Lib"
call Distribute_Files_Work.bat "..\..\..\Protein_Coverage_Summarizer\PeptideToProteinMapper\bin\Debug"
call Distribute_Files_Work.bat "..\..\..\Protein_Coverage_Summarizer\PeptideToProteinMapper\bin\Release"
call Distribute_Files_Work.bat "..\..\..\Protein_Coverage_Summarizer\PeptideToProteinMapper\PeptideToProteinMapEngine\Lib"
call Distribute_Files_Work.bat "..\..\..\Protein_Coverage_Summarizer\PeptideToProteinMapper\PeptideToProteinMapEngine\bin\Debug"
call Distribute_Files_Work.bat "..\..\..\Protein_Coverage_Summarizer\PeptideToProteinMapper\PeptideToProteinMapEngine\bin\Release"

call Distribute_Files_Work.bat "..\..\..\Protein_Coverage_Summarizer\PeptideToProteinMapper\PeptideToProteinMapEngine\Lib"
call Distribute_Files_Work.bat "..\..\..\SMAQC\SMAQC\dll"
call Distribute_Files_Work.bat "..\..\..\SMAQC\SMAQC\SMAQC\bin\Debug"

@echo off
echo.
if not exist ..\..\..\..\JoshAldrich echo Directory not found: skip copying to AScore
if not exist ..\..\..\..\JoshAldrich echo.
if not exist ..\..\..\..\JoshAldrich goto SkipAscore

@echo on
call Distribute_Files_Work.bat "..\..\..\..\JoshAldrich\AScore\AScore_DLL\lib"
call Distribute_Files_Work.bat "..\..\..\..\JoshAldrich\AScore\AScore_DLL\bin\AnyCPU\Release"
call Distribute_Files_Work.bat "..\..\..\..\JoshAldrich\AScore\AScore_DLL\bin\AnyCPU\Debug"
call Distribute_Files_Work.bat "..\..\..\..\JoshAldrich\AScore\AScore_Console\bin\Debug"
call Distribute_Files_Work.bat "..\..\..\..\JoshAldrich\AScore\AScore_Console\bin\Release"

@echo off
:SkipAscore

if "%1"=="NoPause" Goto Exit
if not "%1"=="NoPause" pause

:Exit
