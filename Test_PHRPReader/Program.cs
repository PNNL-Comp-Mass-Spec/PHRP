using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using PHRPReader;

namespace Test_PHRPReader
{
    static class Program
    {
        public static void Main()
        {
            //const string synFileSequest = @"Seq201304121552_Auto934225\Firestone_Soil_07_18_05APR13_Frodo_12-12-04_syn.txt";
            //const string synFileMSGFPlus = @"MSG201304261714_Auto938181\FSFA-299b_25Apr13_Methow_13-02-13_msgfdb_syn.txt";

            const string resultsPathSequest = @"\\proto-7\VOrbi05\2013_2\Firestone_Soil_07_18_05APR13_Frodo_12-12-04\Seq201304121552_Auto934225";
            // const string resultsPathMSGFPlus = "MSG201304261714_Auto938181"

            // const string resultsPathMSGFPlus = "\\proto-7\VOrbiETD03\2015_1\proteogeomics_32_crude_heavy_peptides_200f_25Feb15_Tiger_15-01-26\MSG201503091410_Auto1169297"
            const string resultsPathMSGFPlus = @"C:\DMS_WorkDir";

            //const string resultsXTandem = @"\\proto-7\VOrbiETD01\2013_3\QC_Shew_13_04_pt1_1_2_27Jun13_Leopard_13-05-20\XTM201307011524_Auto958319"

            //const string resultsMSAlign = @"\\proto-9\VOrbiETD02\2014_1\Synocho_D2_2\MSA201402281500_Auto1030272"

            string synOrFhtFile;
            var matchedResultType = clsPHRPReader.ePeptideHitResultType.Unknown;

            if (false) {
                Console.WriteLine();
                synOrFhtFile = clsPHRPReader.AutoDetermineBestInputFile(resultsPathSequest, out matchedResultType);
                if (!string.IsNullOrEmpty(synOrFhtFile) && matchedResultType != clsPHRPReader.ePeptideHitResultType.Unknown) {
                    TestPHRPReader(synOrFhtFile, blnSkipDuplicates: true);
                }
            }

            var msgfPlusDirectory = new DirectoryInfo(resultsPathMSGFPlus);
            if (!msgfPlusDirectory.Exists) {
                Console.WriteLine("Warning, Folder not found: " + resultsPathMSGFPlus);
            }

            if (true & msgfPlusDirectory.Exists) {
                Console.WriteLine();
                synOrFhtFile = clsPHRPReader.AutoDetermineBestInputFile(resultsPathMSGFPlus, out matchedResultType);
                if (!string.IsNullOrEmpty(synOrFhtFile) && matchedResultType != clsPHRPReader.ePeptideHitResultType.Unknown) {
                    TestPHRPReader(synOrFhtFile, blnSkipDuplicates: false);
                }
            }

            if (true & msgfPlusDirectory.Exists) {
                Console.WriteLine();
                synOrFhtFile = clsPHRPReader.AutoDetermineBestInputFile(resultsPathMSGFPlus, out matchedResultType);
                if (!string.IsNullOrEmpty(synOrFhtFile) && matchedResultType != clsPHRPReader.ePeptideHitResultType.Unknown) {
                    TestPHRPReader(synOrFhtFile, blnSkipDuplicates: true);
                }
            }

            if (true & msgfPlusDirectory.Exists) {
                // Look for an MSGF+ parameter file to parse
                var lstFiles = msgfPlusDirectory.GetFiles("MSGFDB*.txt");
                if (lstFiles.Length > 0) {
                    TestMSGFPlusParamFileParsing(lstFiles.First().FullName);
                }
            }

            var sequestDirectory = new DirectoryInfo(resultsPathSequest);
            if (!sequestDirectory.Exists) {
                Console.WriteLine("Warning, Folder not found: " + resultsPathSequest);
            }

            if (!sequestDirectory.Exists)
                return;

            var startTimeNoSkipDup = default(DateTime);
            var endTimeNoSkipDup = default(DateTime);

            var startTimeSkipDup = default(DateTime);
            var endTimeSkipDup = default(DateTime);

            Console.WriteLine();
            synOrFhtFile = clsPHRPReader.AutoDetermineBestInputFile(resultsPathSequest);
            if (!string.IsNullOrEmpty(synOrFhtFile) && matchedResultType != clsPHRPReader.ePeptideHitResultType.Unknown) {
                startTimeNoSkipDup = DateTime.UtcNow;
                TestPHRPReader(synOrFhtFile, blnSkipDuplicates: false);
                endTimeNoSkipDup = DateTime.UtcNow;
            }

            Console.WriteLine();
            synOrFhtFile = clsPHRPReader.AutoDetermineBestInputFile(resultsPathSequest, out matchedResultType);
            if (!string.IsNullOrEmpty(synOrFhtFile) && matchedResultType != clsPHRPReader.ePeptideHitResultType.Unknown) {
                startTimeSkipDup = DateTime.UtcNow;
                TestPHRPReader(synOrFhtFile, blnSkipDuplicates: true);
                endTimeSkipDup = DateTime.UtcNow;
            }

            Console.WriteLine();

            Console.WriteLine("Elapsed time (Keep Duplicates): " + endTimeNoSkipDup.Subtract(startTimeNoSkipDup).TotalSeconds.ToString("0.0") + " seconds");
            Console.WriteLine("Elapsed time (Skip Duplicates): " + endTimeSkipDup.Subtract(startTimeSkipDup).TotalSeconds.ToString("0.0") + " seconds");
            Console.WriteLine();
        }

        private static void TestMSGFPlusParamFileParsing(string msgfPlusParamFilePath)
        {
            var modFileProcessor = new clsMSGFPlusParamFileModExtractor("MSGF+");

            modFileProcessor.ErrorEvent += ErrorEventHandler;
            modFileProcessor.WarningEvent += WarningEventHandler;

            var peptideMassCalculator = new clsPeptideMassCalculator();
            clsPHRPParserMSGFDB.UpdateMassCalculatorMasses(msgfPlusParamFilePath, modFileProcessor, peptideMassCalculator, out _);

            var udtModInfo = new List<clsPeptideMassCalculator.udtPeptideSequenceModInfoType>();

            var monoMass = peptideMassCalculator.ComputeSequenceMass("PEPTIDE", udtModInfo);

            Console.WriteLine("Mono mass of PEPTIDE: " + monoMass.ToString("0.0000"));

            var udtModifiedResidue = new clsPeptideMassCalculator.udtPeptideSequenceModInfoType
            {
                ResidueLocInPeptide = 4,
                ModificationMass = 79.966
            };
            udtModInfo.Add(udtModifiedResidue);

            var monoMassModified = peptideMassCalculator.ComputeSequenceMass("PEPTIDE", udtModInfo);

            Console.WriteLine("Mono mass of PEPT*IDE: " + monoMassModified.ToString("0.0000"));

            Debug.Assert(Math.Abs(monoMass - 799.359926865) < 1E-07);

            Debug.Assert(Math.Abs(monoMassModified - 879.325926865) < 1E-07);
        }

        private static void TestPHRPReader(string synOrFhtFile, bool blnSkipDuplicates)
        {
            var fiInputFile = new FileInfo(synOrFhtFile);

            Console.WriteLine("Instantiating reader");
            var oStartupOptions = new clsPHRPStartupOptions
            {
                LoadModsAndSeqInfo = true,
                LoadMSGFResults = true,
                LoadScanStatsData = false,
                MaxProteinsPerPSM = 100
            };

            var phrpReader =
                new clsPHRPReader(fiInputFile.FullName, clsPHRPReader.ePeptideHitResultType.Unknown, oStartupOptions)
                {
                    EchoMessagesToConsole = false,
                    SkipDuplicatePSMs = blnSkipDuplicates
                };

            // Check for any load errors
            if (phrpReader.ErrorMessages.Count > 0) {
                Console.WriteLine("Error(s) instantiating the reader:");
                foreach (var errorMessage in phrpReader.ErrorMessages) {
                    Console.WriteLine("  " + errorMessage);
                }
            }

            phrpReader.ErrorEvent += ErrorEventHandler;
            phrpReader.StatusEvent += MessageEventHandler;
            phrpReader.WarningEvent += WarningEventHandler;

            const bool fastReadEnabled = true;
            phrpReader.FastReadMode = fastReadEnabled;

            var oMassCalculator = new clsPeptideMassCalculator();

            if (!phrpReader.CanRead) {
                Console.WriteLine("Aborting since PHRPReader is not ready: " + phrpReader.ErrorMessage);
                return;
            }

            var lstValues = new List<string>();

            var intPSMsRead = 0;
            var intModifiedPSMsRead = 0;

            var dctCachedValues = new Dictionary<int, clsPSM>();

            Console.WriteLine("Reading data");

            while (phrpReader.MoveNext()) {
                var oPsm = phrpReader.CurrentPSM;

                intPSMsRead += 1;
                lstValues.Clear();

                phrpReader.FinalizeCurrentPSM();

                clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(oPsm.Peptide, out _, out _, out _);

                var strMassErrorPPM = GetCorrectedMassErrorPPM(oPsm, out _);

                lstValues.Add(phrpReader.DatasetName + "_dta.txt");                                         // #SpecFile
                lstValues.Add("index=" + intPSMsRead);                                                      // SpecID
                lstValues.Add(oPsm.ScanNumber.ToString());                                                      // ScanNum
                lstValues.Add(oPsm.CollisionMode);                                                              // FragMethod
                lstValues.Add(oMassCalculator.ConvoluteMass(oPsm.PrecursorNeutralMass, 0, oPsm.Charge).ToString(CultureInfo.InvariantCulture));      // Precursor m/z

                lstValues.Add(strMassErrorPPM);                                                                 // PrecursorError(ppm)
                lstValues.Add(oPsm.Charge.ToString());                                                          // Charge
                lstValues.Add(oPsm.NumTrypticTerminii.ToString());                                              // Tryptic state (0, 1, or 2)
                lstValues.Add(CleanupPeptide(oPsm.PeptideWithNumericMods));                                     // Peptide

                if (oPsm.SeqID <= 0) {
                    lstValues.Add("**" + oPsm.SeqID + "**");                                                // SeqID is undefined
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
                lstValues.Add(oPsm.MSGFSpecEValue);                                                   // MSGF SpecEValue
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
            if (clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strPeptide, out var strPrimarySequence, out var strPrefix, out var strSuffix))
            {
                // Look for an N-terminal iTraq mod
                var reMatch = RegexFindItraq.Match(strPrimarySequence);

                if (reMatch.Success) {
                    strPeptide = strPrefix + "." + reMatch.Groups[2].Value + reMatch.Groups[1].Value + reMatch.Groups[3].Value + "." + strSuffix;
                }
            }

            return strPeptide;
        }

        private static string FlattenList(IReadOnlyList<string> lstValues, string chSepChar = "\t")
        {
            return string.Join(chSepChar, lstValues);
        }

        private static string GetCorrectedMassErrorPPM(clsPSM oPsm, out int intIsotopeError)
        {
            const double MASS_C13 = 1.00335483;

            double dblMassErrorPPM = 0;
            intIsotopeError = 0;

            if (double.TryParse(oPsm.MassErrorDa, out var dblDelM)) {
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
            if (!oPsm.TryGetScore(strScoreName, out var strScoreValue)) {
                strScoreValue = strValueIfMissing;
            }

            return strScoreValue;
        }

        private static void ErrorEventHandler(string message, Exception ex)
        {
            PRISM.ConsoleMsgUtils.ShowError(message);
        }

        private static void MessageEventHandler(string message)
        {
            Console.WriteLine(message);
        }

        private static void WarningEventHandler(string message)
        {
            PRISM.ConsoleMsgUtils.ShowWarning(message);
        }

    }
}
