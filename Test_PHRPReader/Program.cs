using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Text.RegularExpressions;
using PHRPReader;

namespace Test_PHRPReader
{
    static class Program
    {
        private static clsPHRPReader withEventsField_mPHRPReader;
        private static clsPHRPReader mPHRPReader {
            get { return withEventsField_mPHRPReader; }
            set {
                if (withEventsField_mPHRPReader != null) {
                    withEventsField_mPHRPReader.ErrorEvent -= mPHRPReader_ErrorEvent;
                    withEventsField_mPHRPReader.MessageEvent -= mPHRPReader_MessageEvent;
                    withEventsField_mPHRPReader.WarningEvent -= mPHRPReader_WarningEvent;
                }
                withEventsField_mPHRPReader = value;
                if (withEventsField_mPHRPReader != null) {
                    withEventsField_mPHRPReader.ErrorEvent += mPHRPReader_ErrorEvent;
                    withEventsField_mPHRPReader.MessageEvent += mPHRPReader_MessageEvent;
                    withEventsField_mPHRPReader.WarningEvent += mPHRPReader_WarningEvent;
                }
            }
        }

        public static void Main()
        {
            //const string strSequestSynFilePath = @"Seq201304121552_Auto934225\Firestone_Soil_07_18_05APR13_Frodo_12-12-04_syn.txt";
            //const string strMSGFPlusSynFilePath = @"MSG201304261714_Auto938181\FSFA-299b_25Apr13_Methow_13-02-13_msgfdb_syn.txt";

            const string strSequestFolder = @"\\proto-7\VOrbi05\2013_2\Firestone_Soil_07_18_05APR13_Frodo_12-12-04\Seq201304121552_Auto934225";
            // const string strMSGFPlusFolder = "MSG201304261714_Auto938181"

            // const string strMSGFPlusFolder = "\\proto-7\VOrbiETD03\2015_1\proteogeomics_32_crude_heavy_peptides_200f_25Feb15_Tiger_15-01-26\MSG201503091410_Auto1169297"
            const string strMSGFPlusFolder = @"C:\DMS_WorkDir";

            //const string strXTandemFolder = @"\\proto-7\VOrbiETD01\2013_3\QC_Shew_13_04_pt1_1_2_27Jun13_Leopard_13-05-20\XTM201307011524_Auto958319"

            //const string strMSAlignFolder = @"\\proto-9\VOrbiETD02\2014_1\Synocho_D2_2\MSA201402281500_Auto1030272"

            string strSynOrFHTFile = null;
            clsPHRPReader.ePeptideHitResultType eMatchedResultType = default(clsPHRPReader.ePeptideHitResultType);

            if (false) {
                Console.WriteLine();
                strSynOrFHTFile = clsPHRPReader.AutoDetermineBestInputFile(strSequestFolder, out eMatchedResultType);
                if (!string.IsNullOrEmpty(strSynOrFHTFile) && eMatchedResultType != clsPHRPReader.ePeptideHitResultType.Unknown) {
                    TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates: true);
                }
            }

            var diMSGFPlusFolder = new DirectoryInfo(strMSGFPlusFolder);
            if (!diMSGFPlusFolder.Exists) {
                Console.WriteLine("Warning, Folder not found: " + strMSGFPlusFolder);
            }

            if (true & diMSGFPlusFolder.Exists) {
                Console.WriteLine();
                strSynOrFHTFile = clsPHRPReader.AutoDetermineBestInputFile(strMSGFPlusFolder, out eMatchedResultType);
                if (!string.IsNullOrEmpty(strSynOrFHTFile) && eMatchedResultType != clsPHRPReader.ePeptideHitResultType.Unknown) {
                    TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates: false);
                }
            }

            if (true & diMSGFPlusFolder.Exists) {
                Console.WriteLine();
                strSynOrFHTFile = clsPHRPReader.AutoDetermineBestInputFile(strMSGFPlusFolder, out eMatchedResultType);
                if (!string.IsNullOrEmpty(strSynOrFHTFile) && eMatchedResultType != clsPHRPReader.ePeptideHitResultType.Unknown) {
                    TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates: true);
                }
            }

            if (true & diMSGFPlusFolder.Exists) {
                // Look for an MSGF+ parameter file to parse
                var lstFiles = diMSGFPlusFolder.GetFiles("MSGFDB*.txt");
                if (lstFiles.Length > 0) {
                    TestMSGFPlusParamFileParsing(lstFiles.First().FullName);
                }
            }

            var diSequestFolder = new DirectoryInfo(strSequestFolder);
            if (!diSequestFolder.Exists) {
                Console.WriteLine("Warning, Folder not found: " + strSequestFolder);
            }

            if (true | !diSequestFolder.Exists)
                return;

            DateTime dtStartTimeNoSkipDup = default(DateTime);
            DateTime dtEndTimeNoSkipDup = default(DateTime);

            DateTime dtStartTimeSkipDup = default(DateTime);
            DateTime dtEndTimeSkipDup = default(DateTime);

            Console.WriteLine();
            strSynOrFHTFile = clsPHRPReader.AutoDetermineBestInputFile(strSequestFolder);
            if (!string.IsNullOrEmpty(strSynOrFHTFile) && eMatchedResultType != clsPHRPReader.ePeptideHitResultType.Unknown) {
                dtStartTimeNoSkipDup = DateTime.UtcNow;
                TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates: false);
                dtEndTimeNoSkipDup = DateTime.UtcNow;
            }

            Console.WriteLine();
            strSynOrFHTFile = clsPHRPReader.AutoDetermineBestInputFile(strSequestFolder, out eMatchedResultType);
            if (!string.IsNullOrEmpty(strSynOrFHTFile) && eMatchedResultType != clsPHRPReader.ePeptideHitResultType.Unknown) {
                dtStartTimeSkipDup = DateTime.UtcNow;
                TestPHRPReader(strSynOrFHTFile, blnSkipDuplicates: true);
                dtEndTimeSkipDup = DateTime.UtcNow;
            }

            Console.WriteLine();

            Console.WriteLine("Elapsed time (Keep Duplicates): " + dtEndTimeNoSkipDup.Subtract(dtStartTimeNoSkipDup).TotalSeconds.ToString("0.0") + " seconds");
            Console.WriteLine("Elapsed time (Skip Duplicates): " + dtEndTimeSkipDup.Subtract(dtStartTimeSkipDup).TotalSeconds.ToString("0.0") + " seconds");
        }

        private static void TestMSGFPlusParamFileParsing(string msgfPlusParamFilePath)
        {
            string localErrorMsg = string.Empty;
            var modFileProcessor = new clsMSGFPlusParamFileModExtractor("MSGF+");

            modFileProcessor.ErrorOccurred += ModExtractorErrorHandler;
            modFileProcessor.WarningMessageEvent += ModExtractorWarningHandler;

            var peptideMassCalculator = new clsPeptideMassCalculator();
            var success = clsPHRPParserMSGFDB.UpdateMassCalculatorMasses(msgfPlusParamFilePath, modFileProcessor, peptideMassCalculator, out localErrorMsg);

            var udtModInfo = new List<clsPeptideMassCalculator.udtPeptideSequenceModInfoType>();

            var monoMass = peptideMassCalculator.ComputeSequenceMass("PEPTIDE", udtModInfo);

            Console.WriteLine("Mono mass of PEPTIDE: " + monoMass.ToString("0.0000"));

            var udtModifiedResidue = new clsPeptideMassCalculator.udtPeptideSequenceModInfoType();
            udtModifiedResidue.ResidueLocInPeptide = 4;
            udtModifiedResidue.ModificationMass = 79.966;
            udtModInfo.Add(udtModifiedResidue);

            var monoMassModified = peptideMassCalculator.ComputeSequenceMass("PEPTIDE", udtModInfo);

            Console.WriteLine("Mono mass of PEPT*IDE: " + monoMassModified.ToString("0.0000"));

            Debug.Assert(Math.Abs(monoMass - 799.359926865) < 1E-07);

            Debug.Assert(Math.Abs(monoMassModified - 879.325926865) < 1E-07);
        }

        private static void TestPHRPReader(string strSynOrFHTFile, bool blnSkipDuplicates)
        {
            FileInfo fiInputFile = default(FileInfo);
            fiInputFile = new FileInfo(strSynOrFHTFile);

            Console.WriteLine("Instantiating reader");
            var oStartupOptions = new clsPHRPStartupOptions();
            oStartupOptions.LoadModsAndSeqInfo = true;
            oStartupOptions.LoadMSGFResults = true;
            oStartupOptions.LoadScanStatsData = false;
            oStartupOptions.MaxProteinsPerPSM = 100;

            mPHRPReader = new clsPHRPReader(fiInputFile.FullName, clsPHRPReader.ePeptideHitResultType.Unknown, oStartupOptions);
            mPHRPReader.EchoMessagesToConsole = false;
            mPHRPReader.SkipDuplicatePSMs = blnSkipDuplicates;

            // Check for any load errors
            if (mPHRPReader.ErrorMessages.Count > 0) {
                Console.WriteLine("Error(s) instantiating the reader:");
                foreach (var errorMessage in mPHRPReader.ErrorMessages) {
                    Console.WriteLine("  " + errorMessage);
                }
            }

            const bool fastReadEnabled = true;
            mPHRPReader.FastReadMode = fastReadEnabled;

            var oMassCalculator = new clsPeptideMassCalculator();

            if (!mPHRPReader.CanRead) {
                Console.WriteLine("Aborting since PHRPReader is not ready: " + mPHRPReader.ErrorMessage);
                return;
            }

            var lstValues = new List<string>();
            int intIsotopeErrorComputed = 0;
            string strMassErrorPPM = null;

            var intPSMsRead = 0;
            var intModifiedPSMsRead = 0;

            var dctCachedValues = new Dictionary<int, clsPSM>();

            Console.WriteLine("Reading data");

            while (mPHRPReader.MoveNext()) {
                clsPSM oPsm = mPHRPReader.CurrentPSM;

                intPSMsRead += 1;
                lstValues.Clear();

                mPHRPReader.FinalizeCurrentPSM();

                var primarySequence = string.Empty;
                var prefix = string.Empty;
                var suffix = string.Empty;
                clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(oPsm.Peptide, out primarySequence, out prefix, out suffix);

                intIsotopeErrorComputed = 0;
                strMassErrorPPM = GetCorrectedMassErrorPPM(oPsm, ref intIsotopeErrorComputed);

                lstValues.Add(mPHRPReader.DatasetName + "_dta.txt");                                             // #SpecFile
                lstValues.Add("index=" + intPSMsRead);                                                           // SpecID
                lstValues.Add(oPsm.ScanNumber.ToString());                                                       // ScanNum
                lstValues.Add(oPsm.CollisionMode);                                                               // FragMethod
                lstValues.Add(oMassCalculator.ConvoluteMass(oPsm.PrecursorNeutralMass, 0, oPsm.Charge).ToString());      // Precursor m/z

                lstValues.Add(strMassErrorPPM);                                                                  // PrecursorError(ppm)
                lstValues.Add(oPsm.Charge.ToString());                                                           // Charge
                lstValues.Add(oPsm.NumTrypticTerminii.ToString());                                               // Tryptic state (0, 1, or 2)
                lstValues.Add(CleanupPeptide(oPsm.PeptideWithNumericMods));                                      // Peptide

                if (oPsm.SeqID <= 0) {
                    lstValues.Add("**" + oPsm.SeqID + "**");                                                     // SeqID is undefined
                } else {
                    lstValues.Add(oPsm.SeqID.ToString());                                                        // SeqID
                }

                lstValues.Add(oPsm.ProteinFirst);                                                                // Protein First

                if (oPsm.ProteinDetails.Count > 0) {
                    var oFirstProteinDetail = oPsm.ProteinDetails.First();                                       // Protein Details first

                    if (!string.Equals(oPsm.ProteinFirst, oFirstProteinDetail.Key)) {
                        lstValues.Add(oFirstProteinDetail.Key);
                    } else {
                        lstValues.Add("<Match>");
                    }
                    lstValues.Add(oFirstProteinDetail.Value.ResidueStart.ToString());
                    lstValues.Add(oFirstProteinDetail.Value.ResidueEnd.ToString());
                }

                var strXCorr = GetScore(oPsm, clsPHRPParserSequest.DATA_COLUMN_XCorr, "0");
                lstValues.Add(strXCorr);                                                              // XCorr

                lstValues.Add(GetScore(oPsm, clsPHRPParserSequest.DATA_COLUMN_Sp, "0"));              // SP
                lstValues.Add(oPsm.MSGFSpecProb);                // MSGF SpecProb
                lstValues.Add(GetScore(oPsm, clsPHRPParserSequest.DATA_COLUMN_DelCn2, "0"));          // DelCn2

                lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_PValue, "0"));           // PValue
                lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_EValue, "0"));           // EValue
                lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_Rank_MSGFPlus_SpecEValue, "0"));           // SpecEValue
                lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFDB.DATA_COLUMN_FDR, "1"));              // FDR

                if (oPsm.PeptideCleanSequence == "QQIEESTSDYDKEK") {
                    Console.WriteLine(oPsm.Peptide + " in scan " + oPsm.ScanNumber);

                    var parentIonMZ = oMassCalculator.ConvoluteMass(oPsm.PrecursorNeutralMass, 0, oPsm.Charge);

                    Console.WriteLine("ParentIonMZ   = " + parentIonMZ);
                    Console.WriteLine("PeptideWithNumericMods   = " + oPsm.PeptideWithNumericMods);
                }

                if (oPsm.ModifiedResidues.Count > 0) {
                    intModifiedPSMsRead += 1;

                    if (intModifiedPSMsRead % 500 == 0) {
                        Console.WriteLine("PeptideWithNumericMods   = " + oPsm.PeptideWithNumericMods);
                        foreach (var modifiedResidue in oPsm.ModifiedResidues) {
                            Console.WriteLine("  " + modifiedResidue.Residue + modifiedResidue.EndResidueLocInPeptide + ": " + modifiedResidue.ModDefinition.ModificationMassAsText);
                        }
                    }

                    var dblPeptideMassRecomputed = oMassCalculator.ComputeSequenceMassNumericMods(oPsm.PeptideWithNumericMods);
                    if (Math.Abs(oPsm.PeptideMonoisotopicMass - dblPeptideMassRecomputed) > 0.1) {
                        Console.WriteLine("  Peptide mass disagreement: " + (oPsm.PeptideMonoisotopicMass - dblPeptideMassRecomputed).ToString("0.0000000"));
                    }
                }

                var strFlattened = FlattenList(lstValues);

                if (intPSMsRead % 1000 == 0) {
                    //Console.WriteLine(intPSMsRead.ToString().PadRight(8) & " " & oPsm.Peptide.PadRight(40) & "   " & strXCorr)
                    Console.WriteLine(strFlattened);
                }

                dctCachedValues.Add(intPSMsRead, oPsm);
            }
        }

        private static readonly Regex RegexFindItraq = new Regex(@"^([A-Z][^A-Z]*)(\+144\.\d+)(.+)", RegexOptions.Compiled | RegexOptions.IgnoreCase);
        private static string CleanupPeptide(string strPeptide)
        {
            var strPrimarySequence = string.Empty;
            var strPrefix = string.Empty;
            var strSuffix = string.Empty;

            Match reMatch = default(Match);

            if (clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strPeptide, out strPrimarySequence, out strPrefix, out strSuffix)) {
                // Look for an N-terminal iTraq mod
                reMatch = RegexFindItraq.Match(strPrimarySequence);

                if (reMatch.Success) {
                    strPeptide = strPrefix + "." + reMatch.Groups[2].Value + reMatch.Groups[1].Value + reMatch.Groups[3].Value + "." + strSuffix;
                }
            }

            return strPeptide;
        }

        private static string FlattenList(List<string> lstValues)
        {
            return FlattenList(lstValues, '\t');
        }

        private static string FlattenList(List<string> lstValues, char chSepChar)
        {
            var sbOutline = new StringBuilder();

            for (var intIndex = 0; intIndex <= lstValues.Count - 1; intIndex++) {
                if (intIndex > 0) {
                    sbOutline.Append(chSepChar);
                }
                sbOutline.Append(lstValues[intIndex]);
            }

            return sbOutline.ToString();
        }

        private static string GetCorrectedMassErrorPPM(clsPSM oPsm, ref int intIsotopeError)
        {
            const double MASS_C13 = 1.00335483;

            double dblDelM = 0;
            double dblMassErrorPPM = 0;
            intIsotopeError = 0;

            if (double.TryParse(oPsm.MassErrorDa, out dblDelM)) {
                // Examine dblDelM to determine which isotope was chosen
                if (dblDelM >= -0.5) {
                    // This is the typical case
                    while (dblDelM > 0.5) {
                        dblDelM -= MASS_C13;
                        intIsotopeError += 1;
                    }
                } else {
                    // This happens less often; but we'll still account for it
                    // In this case, intCorrectionCount will be negative
                    while (dblDelM < -0.5) {
                        dblDelM += MASS_C13;
                        intIsotopeError -= 1;
                    }
                }

                dblMassErrorPPM = clsPeptideMassCalculator.MassToPPM(dblDelM, oPsm.PrecursorNeutralMass);
            }

            return dblMassErrorPPM.ToString("0.0000");
        }

        private static string GetScore(clsPSM oPsm, string strScoreName, string strValueIfMissing)
        {
            var strScoreValue = string.Empty;

            if (!oPsm.TryGetScore(strScoreName, out strScoreValue)) {
                strScoreValue = strValueIfMissing;
            }

            return strScoreValue;
        }

        private static void mPHRPReader_ErrorEvent(string strErrorMessage)
        {
            Console.WriteLine("Error: " + strErrorMessage);
        }

        private static void mPHRPReader_MessageEvent(string strMessage)
        {
            Console.WriteLine(strMessage);
        }

        private static void mPHRPReader_WarningEvent(string strWarningMessage)
        {
            Console.WriteLine("Warning: " + strWarningMessage);
        }

        #region "Event Handlers"
        private static void ModExtractorErrorHandler(string errMsg)
        {
            Console.WriteLine("Error: " + errMsg);
        }

        private static void ModExtractorWarningHandler(string warningMsg)
        {
            Console.WriteLine("Warning: " + warningMsg);
        }
        #endregion
    }
}
