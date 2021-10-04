using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Text.RegularExpressions;
using System.Threading;
using PHRPReader;
using PRISM;

namespace CreateMSGFPlusResultsFileFromPHRP
{
    // This program reads a PHRP-compatible _msgfplus_fht.txt file and creates the
    // equivalent tab-delimited _msgfplus.tsv file that would have been created
    // by MSGFPlus when converting the .mzIdentML file to a .tsv file
    //
    // -------------------------------------------------------------------------------
    // Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
    // Program started March 14, 2013
    //
    // E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
    // Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics
    // -------------------------------------------------------------------------------

    internal static class Program
    {
        // Ignore Spelling: Traq, msgfdb, fht, Frag, Novo, delM

        public const string PROGRAM_DATE = "February 23, 2021";
        private const double MASS_C13 = 1.00335483d;

        // ReSharper disable once InconsistentNaming
        private static readonly Regex mItraqMatcher = new(@"^([A-Z][^A-Z]*)(\+144\.\d+)(.+)", RegexOptions.Compiled | RegexOptions.IgnoreCase);

        private static string mInputFilePath;
        private static string mOutputFilePath;

        public static int Main()
        {
            var commandLineParser = new clsParseCommandLine();

            mInputFilePath = string.Empty;
            mOutputFilePath = string.Empty;

            try
            {
                var proceed = false;
                if (commandLineParser.ParseCommandLine())
                {
                    if (SetOptionsUsingCommandLineParameters(commandLineParser))
                    {
                        proceed = true;
                    }
                }

                if (!proceed ||
                    commandLineParser.NeedToShowHelp ||
                    commandLineParser.ParameterCount + commandLineParser.NonSwitchParameterCount == 0 ||
                    mInputFilePath.Length == 0)
                {
                    ShowProgramHelp();
                    return -1;
                }

                var success = ConvertFile();

                return success ? 0 : -1;
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error occurred in Program->Main", ex);
                return -1;
            }
        }

        private static string CleanupPeptide(string peptide)
        {
            var prefixAndSuffixFound = clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(
                peptide, out var primarySequence, out var prefix, out var suffix);

            if (!prefixAndSuffixFound)
            {
                return peptide;
            }

            // Look for an N-terminal iTraq mod
            var match = mItraqMatcher.Match(primarySequence);
            if (match.Success)
            {
                peptide = prefix + "." + match.Groups[2].Value + match.Groups[1].Value + match.Groups[3].Value + "." + suffix;
            }

            return peptide;
        }

        private static bool ConvertFile()
        {
            try
            {
                var inputFile = new FileInfo(mInputFilePath);

                if (!inputFile.Exists)
                {
                    ShowErrorMessage("Input file not found: " + inputFile.FullName);
                    return false;
                }

                var inputFileDirectoryPath = inputFile.Directory != null ? inputFile.Directory.FullName : string.Empty;

                var datasetName = GetDatasetName(inputFile.Name);

                if (string.IsNullOrEmpty(mOutputFilePath))
                {
                    // Auto-define the output file
                    mOutputFilePath = Path.Combine(inputFileDirectoryPath, datasetName + "_msgfplus.tsv");
                }

                var scanStatsFilePath = Path.Combine(inputFileDirectoryPath, datasetName + "_ScanStats.txt");
                var scanStatsFile = new FileInfo(scanStatsFilePath);

                var phrpOptions = new clsPHRPStartupOptions
                {
                    LoadModsAndSeqInfo = true,
                    LoadMSGFResults = false,
                    LoadScanStatsData = scanStatsFile.Exists
                };

                Console.WriteLine("Opening " + PathUtils.CompactPathString(inputFile.FullName, 100));
                Console.WriteLine();

                var phrpReader = new clsPHRPReader(inputFile.FullName, phrpOptions)
                {
                    EchoMessagesToConsole = false,
                    SkipDuplicatePSMs = false
                };

                phrpReader.StatusEvent += PHRPReader_StatusEvent;
                phrpReader.ErrorEvent += PHRPReader_ErrorEvent;
                phrpReader.WarningEvent += PHRPReader_WarningEvent;

                if (!phrpReader.CanRead)
                {
                    ShowErrorMessage("Aborting since PHRPReader is not ready: " + phrpReader.ErrorMessage);
                    return false;
                }

                var massCalculator = new clsPeptideMassCalculator();

                using var writer = new StreamWriter(new FileStream(mOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read));

                var values = new List<string>();
                var itemsRead = 0;

                // Write the header line
                var headerNames = new List<string> {
                    "#SpecFile", "SpecId", "ScanNum", "ScanTime(Min)", "FragMethod",
                    "Precursor", "IsotopeError", "PrecursorError(ppm)", "Charge",
                    "Peptide", "Protein", "DeNovoScore", "MSGFScore",
                    "SpecEValue", "EValue", "QValue", "PepQValue" };

                writer.WriteLine(FlattenList(headerNames));

                while (phrpReader.MoveNext())
                {
                    var psm = phrpReader.CurrentPSM;
                    itemsRead++;
                    values.Clear();

                    // ReSharper disable once InconsistentNaming
                    var massErrorPPM = GetCorrectedMassErrorPPM(psm, out var isotopeErrorComputed);

                    var elutionTimeMinutes = psm.ElutionTimeMinutes > 0 ? psm.ElutionTimeMinutes.ToString("0.00") : "0";

                    values.Add(phrpReader.DatasetName + ".mzML");                                               // #SpecFile
                    values.Add("index=" + itemsRead);                                                           // SpecID
                    values.Add(psm.ScanNumber.ToString());                                                      // ScanNum
                    values.Add(elutionTimeMinutes);                                                             // ScanTime (unknown)
                    values.Add(psm.CollisionMode);                                                              // FragMethod
                    values.Add(GetPrecursorMZ(massCalculator, psm));                                            // Precursor m/z

                    var isotopeError = GetScore(psm, clsPHRPParserMSGFPlus.DATA_COLUMN_Isotope_Error, "0");
                    if (isotopeError == "0" && isotopeErrorComputed != 0)
                    {
                        isotopeError = isotopeErrorComputed.ToString();
                    }

                    values.Add(isotopeError);                                                                   // IsotopeError
                    values.Add(massErrorPPM);                                                                   // PrecursorError(ppm)
                    values.Add(psm.Charge.ToString());                                                          // Charge
                    values.Add(CleanupPeptide(psm.PeptideWithNumericMods));                                     // Peptide
                    values.Add(psm.ProteinFirst);                                                               // Protein
                    values.Add(GetScore(psm, clsPHRPParserMSGFPlus.DATA_COLUMN_DeNovoScore, "0"));              // DeNovoScore
                    values.Add(GetScore(psm, clsPHRPParserMSGFPlus.DATA_COLUMN_MSGFScore, "0"));                // MSGFScore
                    values.Add(GetScore(psm, clsPHRPParserMSGFPlus.DATA_COLUMN_MSGFPlus_SpecEValue, "0"));      // SpecEValue
                    values.Add(GetScore(psm, clsPHRPParserMSGFPlus.DATA_COLUMN_EValue, "0"));                   // EValue
                    values.Add(GetScore(psm, clsPHRPParserMSGFPlus.DATA_COLUMN_QValue, "0"));                   // QValue
                    values.Add(GetScore(psm, clsPHRPParserMSGFPlus.DATA_COLUMN_PepQValue, "0"));                // PepQValue
                    writer.WriteLine(FlattenList(values));
                }

                Console.WriteLine("Read {0:N0} PSMs from the input file", itemsRead);
                Console.WriteLine();

                Console.WriteLine("Created file " + PathUtils.CompactPathString(mOutputFilePath, 100));
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error occurred in Program->ConvertFile", ex);
                return false;
            }

            return true;
        }

        private static string FlattenList(IReadOnlyList<string> values, char sepChar = '\t')
        {
            var outline = new StringBuilder();
            for (int index = 0, loopTo = values.Count - 1; index <= loopTo; index++)
            {
                if (index > 0)
                {
                    outline.Append(sepChar);
                }

                outline.Append(values[index]);
            }

            return outline.ToString();
        }

        private static string GetAppVersion()
        {
            return Assembly.GetExecutingAssembly().GetName().Version + " (" + PROGRAM_DATE + ")";
        }

        private static string GetCorrectedMassErrorPPM(clsPSM psm, out int isotopeError)
        {
            isotopeError = 0;

            if (!double.TryParse(psm.MassErrorDa, out var delM))
            {
                return "0.0000";
            }

            // Examine delM to determine which isotope was chosen
            if (delM >= -0.5d)
            {
                // This is the typical case
                while (delM > 0.5d)
                {
                    delM -= MASS_C13;
                    isotopeError++;
                }
            }
            else
            {
                // This happens less often; but we'll still account for it
                // In this case, correctionCount will be negative
                while (delM < -0.5d)
                {
                    delM += MASS_C13;
                    isotopeError--;
                }
            }

            // ReSharper disable once InconsistentNaming
            var massErrorPPM = clsPeptideMassCalculator.MassToPPM(delM, psm.PrecursorNeutralMass);

            return massErrorPPM.ToString("0.0000");
        }

        private static string GetDatasetName(string inputFileName)
        {
            var datasetName = Path.GetFileNameWithoutExtension(inputFileName);

            if (datasetName.EndsWith("_msgfdb_fht", StringComparison.OrdinalIgnoreCase) ||
                datasetName.EndsWith("_msgfdb_syn", StringComparison.OrdinalIgnoreCase))
            {
                return datasetName.Substring(0, datasetName.Length - 11);
            }

            if (datasetName.EndsWith("_msgfplus_fht", StringComparison.OrdinalIgnoreCase) ||
                datasetName.EndsWith("_msgfplus_syn", StringComparison.OrdinalIgnoreCase))
            {
                return datasetName.Substring(0, datasetName.Length - 13);
            }

            return datasetName;
        }

        private static string GetPrecursorMZ(clsPeptideMassCalculator massCalculator, clsPSM psm)
        {
            return massCalculator.ConvoluteMass(psm.PrecursorNeutralMass, 0, psm.Charge).ToString(CultureInfo.InvariantCulture);
        }

        private static string GetScore(clsPSM psm, string scoreName, string valueIfMissing)
        {
            if (!psm.TryGetScore(scoreName, out var scoreValue))
            {
                scoreValue = valueIfMissing;
            }

            return scoreValue;
        }

        private static bool SetOptionsUsingCommandLineParameters(clsParseCommandLine commandLineParser)
        {
            // Returns True if no problems; otherwise, returns false

            var validParameters = new List<string> { "I", "O" };

            try
            {
                // Make sure no invalid parameters are present
                if (commandLineParser.InvalidParametersPresent(validParameters))
                {
                    ShowErrorMessage("Invalid command line parameters",
                        (from item in commandLineParser.InvalidParameters(validParameters) select ("/" + item)).ToList());
                    return false;
                }

                // Query commandLineParser to see if various parameters are present
                if (commandLineParser.RetrieveValueForParameter("I", out var inputFilePath))
                {
                    mInputFilePath = string.Copy(inputFilePath);
                }
                else if (commandLineParser.NonSwitchParameterCount > 0)
                {
                    mInputFilePath = commandLineParser.RetrieveNonSwitchParameter(0);
                }

                if (commandLineParser.RetrieveValueForParameter("O", out var outputFilePath))
                {
                    mOutputFilePath = string.Copy(outputFilePath);
                }

                return true;
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error parsing the command line parameters", ex);
            }

            return false;
        }

        private static void ShowErrorMessage(string message, Exception ex = null)
        {
            ConsoleMsgUtils.ShowError(message, ex);
        }

        private static void ShowErrorMessage(string title, IEnumerable<string> errorMessages)
        {
            ConsoleMsgUtils.ShowErrors(title, errorMessages);
        }

        private static void ShowProgramHelp()
        {
            try
            {
                Console.WriteLine(ConsoleMsgUtils.WrapParagraph(
                    "This program reads a PHRP-compatible _msgfplus_fht.txt file and creates the equivalent " +
                    "tab-delimited _msgfplus.tsv file that would have been created by MSGFPlus " +
                    "when converting the .mzIdentML file to a .tsv file"));
                Console.WriteLine();
                Console.WriteLine("Program syntax:");
                Console.WriteLine(Path.GetFileName(Assembly.GetExecutingAssembly().Location) + " InputFilePath [/O:OutputFilePath]");
                Console.WriteLine();
                Console.WriteLine("Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2006");
                Console.WriteLine("Version: " + GetAppVersion());
                Console.WriteLine();
                Console.WriteLine("E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov");
                Console.WriteLine("Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics");

                // Delay for 750 msec in case the user double clicked this file from within Windows Explorer (or started the program via a shortcut)
                Thread.Sleep(750);
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error displaying the program syntax", ex);
            }
        }

        private static void PHRPReader_ErrorEvent(string message, Exception ex)
        {
            ShowErrorMessage(message, ex);
        }

        private static void PHRPReader_StatusEvent(string message)
        {
            Console.WriteLine(message);
        }

        private static void PHRPReader_WarningEvent(string warningMessage)
        {
            Console.WriteLine("Warning: " + warningMessage);
        }
    }
}