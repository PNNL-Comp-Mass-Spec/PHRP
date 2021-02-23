using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.CompilerServices;
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
    // E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
    // Website: http://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/ or http://panomics.pnnl.gov/
    // -------------------------------------------------------------------------------

    static class modMain
    {
        public const string PROGRAM_DATE = "March 14, 2013";
        private const double MASS_C13 = 1.00335483d;
        private static clsPHRPReader _mPHRPReader;

        private static Regex reFindItraq = new (@"^([A-Z][^A-Z]*)(\+144\.\d+)(.+)", RegexOptions.Compiled | RegexOptions.IgnoreCase);

        private static clsPHRPReader mPHRPReader
        {
            [MethodImpl(MethodImplOptions.Synchronized)]
            get
            {
                return _mPHRPReader;
            }

            [MethodImpl(MethodImplOptions.Synchronized)]
            set
            {
                if (_mPHRPReader != null)
                {
                    _mPHRPReader.ErrorEvent -= PHRPReader_ErrorEvent;
                    _mPHRPReader.StatusEvent -= PHRPReader_StatusEvent;
                    _mPHRPReader.WarningEvent -= PHRPReader_WarningEvent;
                }

                _mPHRPReader = value;
                if (_mPHRPReader != null)
                {
                    _mPHRPReader.ErrorEvent += PHRPReader_ErrorEvent;
                    _mPHRPReader.StatusEvent += PHRPReader_StatusEvent;
                    _mPHRPReader.WarningEvent += PHRPReader_WarningEvent;
                }
            }
        }

        private static string mInputFilePath;
        private static string mOutputFilePath;

        public static int Main()
        {
            var commandLineParser = new clsParseCommandLine();
            bool proceed;

            mInputFilePath = string.Empty;
            mOutputFilePath = string.Empty;

            try
            {
                proceed = false;
                if (commandLineParser.ParseCommandLine())
                {
                    if (SetOptionsUsingCommandLineParameters(commandLineParser))
                        proceed = true;
                }

                if (!proceed || commandLineParser.NeedToShowHelp || commandLineParser.ParameterCount + commandLineParser.NonSwitchParameterCount == 0 || mInputFilePath.Length == 0)
                {
                    ShowProgramHelp();
                    return -1;
                }

                ConvertFile();
                return 0;
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error occurred in modMain->Main", ex);
                return -1;
            }
        }

        private static string CleanupPeptide(string peptide)
        {
            if (clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(peptide, out var primarySequence, out var prefix, out var suffix))
            {
                // Look for an N-terminal iTraq mod
                var match = reFindItraq.Match(primarySequence);
                if (match.Success)
                {
                    peptide = prefix + "." + match.Groups[2].Value + match.Groups[1].Value + match.Groups[3].Value + "." + suffix;
                }
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

                if (string.IsNullOrEmpty(mOutputFilePath))
                {
                    // Auto-define the output file

                    mOutputFilePath = Path.GetFileNameWithoutExtension(inputFile.Name);
                    if (mOutputFilePath.ToLower().EndsWith("_msgfplus_fht") || mOutputFilePath.EndsWith("_msgfplus_syn"))
                    {
                        mOutputFilePath = mOutputFilePath.Substring(0, mOutputFilePath.Length - 11);
                    }

                    mOutputFilePath = Path.Combine(inputFile.Directory.FullName, mOutputFilePath + "_msgfplus.tsv");
                }

                mPHRPReader = new clsPHRPReader(inputFile.FullName, clsPHRPReader.PeptideHitResultTypes.Unknown, true, false, false);
                mPHRPReader.EchoMessagesToConsole = false;
                mPHRPReader.SkipDuplicatePSMs = false;
                if (!mPHRPReader.CanRead)
                {
                    ShowErrorMessage("Aborting since PHRPReader is not ready: " + mPHRPReader.ErrorMessage);
                    return false;
                }

                var massCalculator = new clsPeptideMassCalculator();

                using (var writer = new StreamWriter(new FileStream(mOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    var values = new List<string>();
                    var psmsRead = 0;

                    // Write the header line
                    var headerNames = new List<string> {
                        "#SpecFile", "SpecId", "ScanNum", "ScanTime(Min)", "FragMethod",
                        "Precursor", "IsotopeError", "PrecursorError(ppm)", "Charge",
                        "Peptide", "Protein", "DeNovoScore", "MSGFScore",
                        "SpecEValue", "EValue", "QValue", "PepQValue" };

                    writer.WriteLine(FlattenList(headerNames));

                    string massErrorPPM;
                    int isotopeErrorComputed;
                    string isotopeError;
                    while (mPHRPReader.MoveNext())
                    {
                        var psm = mPHRPReader.CurrentPSM;
                        psmsRead += 1;
                        values.Clear();
                        isotopeErrorComputed = 0;
                        massErrorPPM = GetCorrectedMassErrorPPM(psm, ref isotopeErrorComputed);
                        values.Add(mPHRPReader.DatasetName + "_dta.txt");                                           // #SpecFile
                        values.Add("index=" + psmsRead);                                                            // SpecID
                        values.Add(psm.ScanNumber.ToString());                                                      // ScanNum
                        values.Add("0");                                                                            // ScanTime (unknown)
                        values.Add(psm.CollisionMode);                                                              // FragMethod
                        values.Add(GetPrecursorMZ(massCalculator, psm));                                            // Precursor m/z
                        isotopeError = GetScore(psm, clsPHRPParserMSGFPlus.DATA_COLUMN_Isotope_Error, "0");
                        if (isotopeError == "0" & isotopeErrorComputed != 0)
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
                }

                Console.WriteLine("Created file " + mOutputFilePath);
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error occurred in modMain->ConvertFile", ex);
                return false;
            }

            return true;
        }

        private static string FlattenList(List<string> values, char sepChar = '\t')
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
            // Return System.Windows.Forms.Application.ProductVersion & " (" & PROGRAM_DATE & ")"

            return Assembly.GetExecutingAssembly().GetName().Version.ToString() + " (" + PROGRAM_DATE + ")";
        }

        private static string GetCorrectedMassErrorPPM(clsPSM psm, ref int isotopeError)
        {
            double delM;
            var massErrorPPM = 0d;
            isotopeError = 0;

            if (double.TryParse(psm.MassErrorDa, out delM))
            {

                // Examine delM to determine which isotope was chosen
                if (delM >= -0.5d)
                {
                    // This is the typical case
                    while (delM > 0.5d)
                    {
                        delM -= MASS_C13;
                        isotopeError += 1;
                    }
                }
                else
                {
                    // This happens less often; but we'll still account for it
                    // In this case, correctionCount will be negative
                    while (delM < -0.5d)
                    {
                        delM += MASS_C13;
                        isotopeError -= 1;
                    }
                }

                massErrorPPM = clsPeptideMassCalculator.MassToPPM(delM, psm.PrecursorNeutralMass);
            }

            return massErrorPPM.ToString("0.0000");
        }

        private static string GetPrecursorMZ(clsPeptideMassCalculator massCalculator, clsPSM psm)
        {
            return massCalculator.ConvoluteMass(psm.PrecursorNeutralMass, 0, psm.Charge).ToString();
        }

        private static string GetScore(clsPSM psm, string scoreName, string valueIfMissing)
        {
            var scoreValue = string.Empty;
            if (!psm.TryGetScore(scoreName, out scoreValue))
            {
                scoreValue = valueIfMissing;
            }

            return scoreValue;
        }

        private static bool SetOptionsUsingCommandLineParameters(clsParseCommandLine commandLineParser)
        {
            // Returns True if no problems; otherwise, returns false

            var value = string.Empty;
            var validParameters = new List<string> { "I", "O" };
            try
            {
                // Make sure no invalid parameters are present
                if (commandLineParser.InvalidParametersPresent(validParameters))
                {
                    ShowErrorMessage("Invalid command line parameters", (from item in commandLineParser.InvalidParameters(validParameters)
                                                                         select ("/" + item)).ToList());
                    return false;
                }
                else
                {
                    // Query commandLineParser to see if various parameters are present
                    if (commandLineParser.RetrieveValueForParameter("I", out value))
                    {
                        mInputFilePath = string.Copy(value);
                    }
                    else if (commandLineParser.NonSwitchParameterCount > 0)
                    {
                        mInputFilePath = commandLineParser.RetrieveNonSwitchParameter(0);
                    }

                    if (commandLineParser.RetrieveValueForParameter("O", out value))
                        mOutputFilePath = string.Copy(value);
                    return true;
                }
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

        private static void ShowErrorMessage(string title, List<string> errorMessages)
        {
            ConsoleMsgUtils.ShowErrors(title, errorMessages);
        }

        private static void ShowProgramHelp()
        {
            try
            {
                Console.WriteLine("This program reads a PHRP-compatible _msgfplus_fht.txt file and creates the equivalent tab-delimited _msgfplus.tsv file that would have been created by MSGFPlus when converting the .mzIdentML file to a .tsv file");
                Console.WriteLine();
                Console.WriteLine("Program syntax:" + Environment.NewLine + Path.GetFileName(Assembly.GetExecutingAssembly().Location) + " InputFilePath [/O:OutputFilePath]");
                Console.WriteLine();
                Console.WriteLine("Program written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in 2006");
                Console.WriteLine("Version: " + GetAppVersion());
                Console.WriteLine();
                Console.WriteLine("E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com");
                Console.WriteLine("Website: http://omics.pnl.gov/ or http://panomics.pnnl.gov/");

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