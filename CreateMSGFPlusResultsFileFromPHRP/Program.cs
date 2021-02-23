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
            int intReturnCode;
            var commandLineParser = new clsParseCommandLine();
            bool blnProceed;
            intReturnCode = 0;
            mInputFilePath = string.Empty;
            mOutputFilePath = string.Empty;
            try
            {
                blnProceed = false;
                if (commandLineParser.ParseCommandLine())
                {
                    if (SetOptionsUsingCommandLineParameters(commandLineParser))
                        blnProceed = true;
                }

                if (!blnProceed || commandLineParser.NeedToShowHelp || commandLineParser.ParameterCount + commandLineParser.NonSwitchParameterCount == 0 || mInputFilePath.Length == 0)
                {
                    ShowProgramHelp();
                    intReturnCode = -1;
                }
                else
                {
                    ConvertFile();
                }
            }
            catch (Exception ex)
            {
                ShowErrorMessage("Error occurred in modMain->Main", ex);
                intReturnCode = -1;
            }

            return intReturnCode;
        }

        private static string CleanupPeptide(string strPeptide)
        {

            var strPrimarySequence = string.Empty;
            var strPrefix = string.Empty;
            var strSuffix = string.Empty;
            Match reMatch;
            if (clsPeptideCleavageStateCalculator.SplitPrefixAndSuffixFromSequence(strPeptide, out strPrimarySequence, out strPrefix, out strSuffix))
            {
                // Look for an N-terminal iTraq mod
                reMatch = reFindItraq.Match(strPrimarySequence);
                if (reMatch.Success)
                {
                    strPeptide = strPrefix + "." + reMatch.Groups[2].Value + reMatch.Groups[1].Value + reMatch.Groups[3].Value + "." + strSuffix;
                }
            }

            return strPeptide;
        }

        private static bool ConvertFile()
        {
            try
            {
                var fiInputFile = new FileInfo(mInputFilePath);
                if (!fiInputFile.Exists)
                {
                    ShowErrorMessage("Input file not found: " + fiInputFile.FullName);
                    return false;
                }

                if (string.IsNullOrEmpty(mOutputFilePath))
                {
                    // Auto-define the output file

                    mOutputFilePath = Path.GetFileNameWithoutExtension(fiInputFile.Name);
                    if (mOutputFilePath.ToLower().EndsWith("_msgfplus_fht") || mOutputFilePath.EndsWith("_msgfplus_syn"))
                    {
                        mOutputFilePath = mOutputFilePath.Substring(0, mOutputFilePath.Length - 11);
                    }

                    mOutputFilePath = Path.Combine(fiInputFile.Directory.FullName, mOutputFilePath + "_msgfplus.tsv");
                }

                mPHRPReader = new clsPHRPReader(fiInputFile.FullName, clsPHRPReader.PeptideHitResultTypes.Unknown, true, false, false);
                mPHRPReader.EchoMessagesToConsole = false;
                mPHRPReader.SkipDuplicatePSMs = false;
                if (!mPHRPReader.CanRead)
                {
                    ShowErrorMessage("Aborting since PHRPReader is not ready: " + mPHRPReader.ErrorMessage);
                    return false;
                }

                var oMassCalculator = new clsPeptideMassCalculator();
                using (var swOutFile = new StreamWriter(new FileStream(mOutputFilePath, FileMode.Create, FileAccess.Write, FileShare.Read)))
                {
                    var lstValues = new List<string>();
                    var intPSMsRead = 0;

                    // Write the header line
                    swOutFile.WriteLine(FlattenList(new List<string>() { "#SpecFile", "SpecID", "ScanNum", "FragMethod", "Precursor", "IsotopeError", "PrecursorError(ppm)", "Charge", "Peptide", "Protein", "DeNovoScore", "MSGFScore", "SpecEValue", "EValue", "QValue", "PepQValue" }));
                    string strMassErrorPPM;
                    int intIsotopeErrorComputed;
                    string strIsotopeError;
                    while (mPHRPReader.MoveNext())
                    {
                        var oPsm = mPHRPReader.CurrentPSM;
                        intPSMsRead += 1;
                        lstValues.Clear();
                        intIsotopeErrorComputed = 0;
                        strMassErrorPPM = GetCorrectedMassErrorPPM(oPsm, ref intIsotopeErrorComputed);
                        lstValues.Add(mPHRPReader.DatasetName + "_dta.txt");                                             // #SpecFile
                        lstValues.Add("index=" + intPSMsRead);                                                           // SpecID
                        lstValues.Add(oPsm.ScanNumber.ToString());                                                       // ScanNum
                        lstValues.Add(oPsm.CollisionMode);                                                               // FragMethod
                        lstValues.Add(GetPrecursorMZ(oMassCalculator, oPsm));                                                             // Precursor m/z
                        strIsotopeError = GetScore(oPsm, clsPHRPParserMSGFPlus.DATA_COLUMN_Isotope_Error, "0");
                        if (strIsotopeError == "0" & intIsotopeErrorComputed != 0)
                        {
                            strIsotopeError = intIsotopeErrorComputed.ToString();
                        }

                        lstValues.Add(strIsotopeError);                                                                  // IsotopeError
                        lstValues.Add(strMassErrorPPM);                                                                  // PrecursorError(ppm)
                        lstValues.Add(oPsm.Charge.ToString());                                                           // Charge
                        lstValues.Add(CleanupPeptide(oPsm.PeptideWithNumericMods));                                      // Peptide
                        lstValues.Add(oPsm.ProteinFirst);                                                                // Protein
                        lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFPlus.DATA_COLUMN_DeNovoScore, "0"));                 // DeNovoScore
                        lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFPlus.DATA_COLUMN_MSGFScore, "0"));                   // MSGFScore
                        lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFPlus.DATA_COLUMN_MSGFPlus_SpecEValue, "0"));           // SpecEValue
                        lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFPlus.DATA_COLUMN_EValue, "0"));                      // EValue
                        lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFPlus.DATA_COLUMN_QValue, "0"));                      // QValue
                        lstValues.Add(GetScore(oPsm, clsPHRPParserMSGFPlus.DATA_COLUMN_PepQValue, "0"));                   // PepQValue
                        swOutFile.WriteLine(FlattenList(lstValues));
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

        private static string FlattenList(List<string> lstValues)
        {
            return FlattenList(lstValues);
        }

        private static string FlattenList(List<string> lstValues, char chSepChar = '\t')
        {
            var sbOutline = new StringBuilder();
            for (int intIndex = 0, loopTo = lstValues.Count - 1; intIndex <= loopTo; intIndex++)
            {
                if (intIndex > 0)
                {
                    sbOutline.Append(chSepChar);
                }

                sbOutline.Append(lstValues[intIndex]);
            }

            return sbOutline.ToString();
        }

        private static string GetAppVersion()
        {
            // Return System.Windows.Forms.Application.ProductVersion & " (" & PROGRAM_DATE & ")"

            return Assembly.GetExecutingAssembly().GetName().Version.ToString() + " (" + PROGRAM_DATE + ")";
        }

        private static string GetCorrectedMassErrorPPM(clsPSM oPsm, ref int intIsotopeError)
        {
            double dblDelM;
            var dblMassErrorPPM = 0d;
            intIsotopeError = 0;
            if (double.TryParse(oPsm.MassErrorDa, out dblDelM))
            {

                // Examine dblDelM to determine which isotope was chosen
                if (dblDelM >= -0.5d)
                {
                    // This is the typical case
                    while (dblDelM > 0.5d)
                    {
                        dblDelM -= MASS_C13;
                        intIsotopeError += 1;
                    }
                }
                else
                {
                    // This happens less often; but we'll still account for it
                    // In this case, intCorrectionCount will be negative
                    while (dblDelM < -0.5d)
                    {
                        dblDelM += MASS_C13;
                        intIsotopeError -= 1;
                    }
                }

                dblMassErrorPPM = clsPeptideMassCalculator.MassToPPM(dblDelM, oPsm.PrecursorNeutralMass);
            }

            return dblMassErrorPPM.ToString("0.0000");
        }

        private static string GetPrecursorMZ(clsPeptideMassCalculator oMassCalculator, clsPSM oPsm)
        {
            return oMassCalculator.ConvoluteMass(oPsm.PrecursorNeutralMass, 0, oPsm.Charge).ToString();
        }

        private static string GetScore(clsPSM oPsm, string strScoreName, string strValueIfMissing)
        {
            var strScoreValue = string.Empty;
            if (!oPsm.TryGetScore(strScoreName, out strScoreValue))
            {
                strScoreValue = strValueIfMissing;
            }

            return strScoreValue;
        }

        private static bool SetOptionsUsingCommandLineParameters(clsParseCommandLine commandLineParser)
        {
            // Returns True if no problems; otherwise, returns false

            var strValue = string.Empty;
            var lstValidParameters = new List<string>() { "I", "O" };
            try
            {
                // Make sure no invalid parameters are present
                if (commandLineParser.InvalidParametersPresent(lstValidParameters))
                {
                    ShowErrorMessage("Invalid command line parameters", (from item in commandLineParser.InvalidParameters(lstValidParameters)
                                                                         select ("/" + item)).ToList());
                    return false;
                }
                else
                {
                    // Query commandLineParser to see if various parameters are present
                    if (commandLineParser.RetrieveValueForParameter("I", out strValue))
                    {
                        mInputFilePath = string.Copy(strValue);
                    }
                    else if (commandLineParser.NonSwitchParameterCount > 0)
                    {
                        mInputFilePath = commandLineParser.RetrieveNonSwitchParameter(0);
                    }

                    if (commandLineParser.RetrieveValueForParameter("O", out strValue))
                        mOutputFilePath = string.Copy(strValue);
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

        private static void PHRPReader_WarningEvent(string strWarningMessage)
        {
            Console.WriteLine("Warning: " + strWarningMessage);
        }
    }
}