using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.RegularExpressions;
using PHRPReader;
using PHRPReader.Data;
using PHRPReader.Reader;

namespace Test_PHRPReader
{
    static class Program
    {
        public static void Main()
        {
            //const string synFileSequest = @"Seq201304121552_Auto934225\Firestone_Soil_07_18_05APR13_Frodo_12-12-04_syn.txt";
            //const string synFileMSGFPlus = @"MSG201304261714_Auto938181\FSFA-299b_25Apr13_Methow_13-02-13_msgfplus_syn.txt";

            const string resultsPathSequest = @"\\proto-7\VOrbi05\2013_2\Firestone_Soil_07_18_05APR13_Frodo_12-12-04\Seq201304121552_Auto934225";
            // const string resultsPathMSGFPlus = "MSG201304261714_Auto938181"

            // const string resultsPathMSGFPlus = "\\proto-7\VOrbiETD03\2015_1\proteogeomics_32_crude_heavy_peptides_200f_25Feb15_Tiger_15-01-26\MSG201503091410_Auto1169297"
            const string resultsPathMSGFPlus = @"C:\DMS_WorkDir";

            //const string resultsXTandem = @"\\proto-7\VOrbiETD01\2013_3\QC_Shew_13_04_pt1_1_2_27Jun13_Leopard_13-05-20\XTM201307011524_Auto958319"

            //const string resultsMSAlign = @"\\proto-9\VOrbiETD02\2014_1\Synocho_D2_2\MSA201402281500_Auto1030272"

            string synOrFhtFile;
            var matchedResultType = PHRPReader.PeptideHitResultTypes.Unknown;

            if (false)
            {
                Console.WriteLine();
                synOrFhtFile = ReaderFactory.AutoDetermineBestInputFile(resultsPathSequest, out matchedResultType);
                if (!string.IsNullOrEmpty(synOrFhtFile) && matchedResultType != PHRPReader.PeptideHitResultTypes.Unknown)
                {
                    TestPHRPReader(synOrFhtFile, blnSkipDuplicates: true);
                }
            }

            var msgfPlusDirectory = new DirectoryInfo(resultsPathMSGFPlus);
            if (!msgfPlusDirectory.Exists)
            {
                Console.WriteLine("Warning, Folder not found: " + resultsPathMSGFPlus);
            }

            if (msgfPlusDirectory.Exists)
            {
                Console.WriteLine();
                synOrFhtFile = ReaderFactory.AutoDetermineBestInputFile(resultsPathMSGFPlus, out matchedResultType);
                if (!string.IsNullOrEmpty(synOrFhtFile) && matchedResultType != PHRPReader.PeptideHitResultTypes.Unknown)
                {
                    TestPHRPReader(synOrFhtFile, blnSkipDuplicates: false);
                }
            }

            if (msgfPlusDirectory.Exists)
            {
                Console.WriteLine();
                synOrFhtFile = ReaderFactory.AutoDetermineBestInputFile(resultsPathMSGFPlus, out matchedResultType);
                if (!string.IsNullOrEmpty(synOrFhtFile) && matchedResultType != PHRPReader.PeptideHitResultTypes.Unknown)
                {
                    TestPHRPReader(synOrFhtFile, blnSkipDuplicates: true);
                }
            }

            if (msgfPlusDirectory.Exists)
            {
                // Look for an MSGF+ parameter file to parse
                var lstFiles = msgfPlusDirectory.GetFiles("MSGFDB*.txt");
                if (lstFiles.Length > 0)
                {
                    TestMSGFPlusParamFileParsing(lstFiles.First().FullName);
                }
            }

            var sequestDirectory = new DirectoryInfo(resultsPathSequest);
            if (!sequestDirectory.Exists)
            {
                Console.WriteLine("Warning, Folder not found: " + resultsPathSequest);
            }

            if (!sequestDirectory.Exists)
                return;

            var startTimeNoSkipDup = default(DateTime);
            var endTimeNoSkipDup = default(DateTime);

            var startTimeSkipDup = default(DateTime);
            var endTimeSkipDup = default(DateTime);

            Console.WriteLine();
            synOrFhtFile = ReaderFactory.AutoDetermineBestInputFile(resultsPathSequest);
            if (!string.IsNullOrEmpty(synOrFhtFile) && matchedResultType != PHRPReader.PeptideHitResultTypes.Unknown)
            {
                startTimeNoSkipDup = DateTime.UtcNow;
                TestPHRPReader(synOrFhtFile, blnSkipDuplicates: false);
                endTimeNoSkipDup = DateTime.UtcNow;
            }

            Console.WriteLine();
            synOrFhtFile = ReaderFactory.AutoDetermineBestInputFile(resultsPathSequest, out matchedResultType);
            if (!string.IsNullOrEmpty(synOrFhtFile) && matchedResultType != PHRPReader.PeptideHitResultTypes.Unknown)
            {
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
            var modFileProcessor = new MSGFPlusParamFileModExtractor("MSGF+");

            modFileProcessor.ErrorEvent += ErrorEventHandler;
            modFileProcessor.WarningEvent += WarningEventHandler;

            var peptideMassCalculator = new PeptideMassCalculator();
            MSGFPlusSynFileReader.UpdateMassCalculatorMasses(msgfPlusParamFilePath, modFileProcessor, peptideMassCalculator, out _);

            var modifiedResidues = new List<PeptideMassCalculator.PeptideSequenceModInfo>();

            var monoMass = peptideMassCalculator.ComputeSequenceMass("PEPTIDE", modifiedResidues);

            Console.WriteLine("Mono mass of PEPTIDE: " + monoMass.ToString("0.0000"));

            var modifiedResidue = new PeptideMassCalculator.PeptideSequenceModInfo
            {
                ResidueLocInPeptide = 4,
                ModificationMass = 79.966
            };
            modifiedResidues.Add(modifiedResidue);

            var monoMassModified = peptideMassCalculator.ComputeSequenceMass("PEPTIDE", modifiedResidues);

            Console.WriteLine("Mono mass of PEPT*IDE: " + monoMassModified.ToString("0.0000"));

            const double MONO_MASS_EXPECTED = 799.359926865;
            const double MODIFIED_MONO_MASS_EXPECTED = 879.325926865;

            if (Math.Abs(monoMass - MONO_MASS_EXPECTED) > 1E-07)
            {
                PRISM.ConsoleMsgUtils.ShowWarning("Computed mass of {0:F6} does not match the expected mass, {1:F6}", monoMass, MONO_MASS_EXPECTED);
            }

            if (Math.Abs(monoMassModified - MODIFIED_MONO_MASS_EXPECTED) > 1E-07)
            {
                PRISM.ConsoleMsgUtils.ShowWarning("Modified computed mass of {0:F6} does not match the expected mass, {1:F6}", monoMassModified, MODIFIED_MONO_MASS_EXPECTED);
            }
        }

        private static void TestPHRPReader(string synOrFhtFile, bool blnSkipDuplicates)
        {
            var inputFile = new FileInfo(synOrFhtFile);

            Console.WriteLine("Instantiating reader");
            var startupOptions = new StartupOptions
            {
                LoadModsAndSeqInfo = true,
                LoadMSGFResults = true,
                LoadScanStatsData = false,
                MaxProteinsPerPSM = 100
            };

            var phrpReader =
                new ReaderFactory(inputFile.FullName, PeptideHitResultTypes.Unknown, startupOptions)
                {
                    EchoMessagesToConsole = false,
                    SkipDuplicatePSMs = blnSkipDuplicates
                };

            // Check for any load errors
            if (phrpReader.ErrorMessages.Count > 0)
            {
                Console.WriteLine("Error(s) instantiating the reader:");
                foreach (var errorMessage in phrpReader.ErrorMessages)
                {
                    Console.WriteLine("  " + errorMessage);
                }
            }

            phrpReader.ErrorEvent += ErrorEventHandler;
            phrpReader.StatusEvent += MessageEventHandler;
            phrpReader.WarningEvent += WarningEventHandler;

            const bool fastReadEnabled = true;
            phrpReader.FastReadMode = fastReadEnabled;

            var massCalculator = new PeptideMassCalculator();

            if (!phrpReader.CanRead)
            {
                Console.WriteLine("Aborting since PHRPReader is not ready: " + phrpReader.ErrorMessage);
                return;
            }

            var lstValues = new List<string>();

            var intPSMsRead = 0;
            var intModifiedPSMsRead = 0;

            // ReSharper disable once CollectionNeverQueried.Local
            var dctCachedValues = new Dictionary<int, PSM>();

            Console.WriteLine("Reading data");

            while (phrpReader.MoveNext())
            {
                var psm = phrpReader.CurrentPSM;

                intPSMsRead += 1;
                lstValues.Clear();

                phrpReader.FinalizeCurrentPSM();

                PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(psm.Peptide, out _, out _, out _);

                var strMassErrorPPM = GetCorrectedMassErrorPPM(psm, out _);

                lstValues.Add(phrpReader.DatasetName + "_dta.txt");         // #SpecFile
                lstValues.Add("index=" + intPSMsRead);                      // SpecID
                lstValues.Add(psm.ScanNumber.ToString());                   // ScanNum
                lstValues.Add(psm.CollisionMode);                           // FragMethod
                lstValues.Add(massCalculator.ConvoluteMass(psm.PrecursorNeutralMass, 0, psm.Charge).ToString(CultureInfo.InvariantCulture));      // Precursor m/z

                lstValues.Add(strMassErrorPPM);                             // PrecursorError(ppm)
                lstValues.Add(psm.Charge.ToString());                       // Charge
                lstValues.Add(psm.NumTrypticTermini.ToString());            // Tryptic state (0, 1, or 2)
                lstValues.Add(CleanupPeptide(psm.PeptideWithNumericMods));  // Peptide

                if (psm.SeqID <= 0)
                {
                    lstValues.Add("**" + psm.SeqID + "**");                 // SeqID is undefined
                }
                else
                {
                    lstValues.Add(psm.SeqID.ToString());                    // SeqID
                }

                lstValues.Add(psm.ProteinFirst);

                if (psm.ProteinDetails.Count > 0)
                {
                    var firstProteinDetail = psm.ProteinDetails.First();

                    if (!string.Equals(psm.ProteinFirst, firstProteinDetail.Key))
                    {
                        lstValues.Add(firstProteinDetail.Key);
                    }
                    else
                    {
                        lstValues.Add("<Match>");
                    }
                    lstValues.Add(firstProteinDetail.Value.ResidueStart.ToString());
                    lstValues.Add(firstProteinDetail.Value.ResidueEnd.ToString());
                }

                var strXCorr = GetScore(psm, SequestSynFileReader.GetColumnNameByID(SequestSynopsisFileColumns.XCorr), "0");
                lstValues.Add(strXCorr);

                lstValues.Add(GetScore(psm, SequestSynFileReader.GetColumnNameByID(SequestSynopsisFileColumns.Sp), "0"));
                lstValues.Add(psm.MSGFSpecEValue);
                lstValues.Add(GetScore(psm, SequestSynFileReader.GetColumnNameByID(SequestSynopsisFileColumns.DeltaCn2), "0"));

                lstValues.Add(GetScore(psm, MSGFPlusSynFileReader.GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.PValue), "0"));
                lstValues.Add(GetScore(psm, MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.EValue), "0"));
                lstValues.Add(GetScore(psm, MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.RankSpecEValue), "0"));
                lstValues.Add(GetScore(psm, MSGFPlusSynFileReader.GetMSGFDBColumnNameByID(MSGFDBSynFileColumns.FDR), "1"));
                lstValues.Add(GetScore(psm, MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.QValue), "0"));
                lstValues.Add(GetScore(psm, MSGFPlusSynFileReader.GetColumnNameByID(MSGFPlusSynFileColumns.PepQValue), "0"));


                if (psm.PeptideCleanSequence == "QQIEESTSDYDKEK")
                {
                    Console.WriteLine(psm.Peptide + " in scan " + psm.ScanNumber);

                    var parentIonMZ = massCalculator.ConvoluteMass(psm.PrecursorNeutralMass, 0, psm.Charge);

                    Console.WriteLine("ParentIonMZ   = " + parentIonMZ);
                    Console.WriteLine("PeptideWithNumericMods   = " + psm.PeptideWithNumericMods);
                }

                if (psm.ModifiedResidues.Count > 0)
                {
                    intModifiedPSMsRead += 1;

                    if (intModifiedPSMsRead % 500 == 0)
                    {
                        Console.WriteLine("PeptideWithNumericMods   = " + psm.PeptideWithNumericMods);
                        foreach (var modifiedResidue in psm.ModifiedResidues)
                        {
                            Console.WriteLine("  " + modifiedResidue.Residue + modifiedResidue.EndResidueLocInPeptide + ": " + modifiedResidue.ModDefinition.ModificationMassAsText);
                        }
                    }

                    var dblPeptideMassRecomputed = massCalculator.ComputeSequenceMassNumericMods(psm.PeptideWithNumericMods);
                    if (Math.Abs(psm.PeptideMonoisotopicMass - dblPeptideMassRecomputed) > 0.1)
                    {
                        Console.WriteLine("  Peptide mass disagreement: " + (psm.PeptideMonoisotopicMass - dblPeptideMassRecomputed).ToString("0.0000000"));
                    }
                }

                var strFlattened = FlattenList(lstValues);

                if (intPSMsRead % 10000 == 0)
                {
                    Console.WriteLine(strFlattened);
                }

                dctCachedValues.Add(intPSMsRead, psm);
            }
        }

        private static readonly Regex RegexFindItraq = new Regex(@"^([A-Z][^A-Z]*)(\+144\.\d+)(.+)", RegexOptions.Compiled | RegexOptions.IgnoreCase);

        private static string CleanupPeptide(string strPeptide)
        {
            if (PeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strPeptide, out var strPrimarySequence, out var strPrefix, out var strSuffix))
            {
                // Look for an N-terminal iTraq mod
                var reMatch = RegexFindItraq.Match(strPrimarySequence);

                if (reMatch.Success)
                {
                    strPeptide = strPrefix + "." + reMatch.Groups[2].Value + reMatch.Groups[1].Value + reMatch.Groups[3].Value + "." + strSuffix;
                }
            }

            return strPeptide;
        }

        private static string FlattenList(IEnumerable<string> lstValues, string chSepChar = "\t")
        {
            return string.Join(chSepChar, lstValues);
        }

        private static string GetCorrectedMassErrorPPM(PSM psm, out int intIsotopeError)
        {
            const double MASS_C13 = 1.00335483;

            double dblMassErrorPPM = 0;
            intIsotopeError = 0;

            if (double.TryParse(psm.MassErrorDa, out var dblDelM))
            {
                // Examine dblDelM to determine which isotope was chosen
                if (dblDelM >= -0.5)
                {
                    // This is the typical case
                    while (dblDelM > 0.5)
                    {
                        dblDelM -= MASS_C13;
                        intIsotopeError += 1;
                    }
                }
                else
                {
                    // This happens less often; but we'll still account for it
                    // In this case, intCorrectionCount will be negative
                    while (dblDelM < -0.5)
                    {
                        dblDelM += MASS_C13;
                        intIsotopeError -= 1;
                    }
                }

                dblMassErrorPPM = PeptideMassCalculator.MassToPPM(dblDelM, psm.PrecursorNeutralMass);
            }

            return dblMassErrorPPM.ToString("0.0000");
        }

        private static string GetScore(PSM psm, string strScoreName, string strValueIfMissing)
        {
            if (!psm.TryGetScore(strScoreName, out var strScoreValue))
            {
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
