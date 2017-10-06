// This class can be used to parse the text following the program name when a
//  program is started from the command line
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started November 8, 2003

// E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com
// Website: http://panomics.pnnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Notice: This computer software was prepared by Battelle Memorial Institute,
// hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the
// Department of Energy (DOE).  All rights in the computer software are reserved
// by DOE on behalf of the United States Government and the Contractor as
// provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY
// WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
// SOFTWARE.  This notice including this sentence must appear on any copies of
// this computer software.

//
// Last modified March 28, 2017

using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading;

namespace PeptideHitResultsProcRunner
{
    public class clsParseCommandLine
    {
        /// <summary>
        /// Default switch char
        /// </summary>
        public const string DEFAULT_SWITCH_CHAR = "/";

        /// <summary>
        /// Alternate switch char
        /// </summary>
        public const string ALTERNATE_SWITCH_CHAR = "-";

        /// <summary>
        /// Default character between the switch name and a value to associate with the parameter
        /// </summary>
        public const string DEFAULT_SWITCH_PARAM_CHAR = ":";

        private readonly Dictionary<string, string> mSwitches = new Dictionary<string, string>();
        private readonly List<string> mNonSwitchParameters = new List<string>();

        private bool mShowHelp = false;

        /// <summary>
        /// If true, we need to show the syntax to the user due to a switch error, invalid switch, or the presence of /? or /help
        /// </summary>
        public bool NeedToShowHelp
        {
            get { return mShowHelp; }
        }

        // ReSharper disable once UnusedMember.Global
        /// <summary>
        /// Number of switches
        /// </summary>
        public int ParameterCount
        {
            get { return mSwitches.Count; }
        }

        // ReSharper disable once UnusedMember.Global
        /// <summary>
        /// Number of parameters that are not preceded by a switch
        /// </summary>
        public int NonSwitchParameterCount
        {
            get { return mNonSwitchParameters.Count; }
        }

        /// <summary>
        /// Set to true to see extra debug information
        /// </summary>
        public bool DebugMode { get; }

        public clsParseCommandLine(bool blnDebugMode = false)
        {
            DebugMode = blnDebugMode;
        }

        /// <summary>
        /// Compares the parameter names in objParameterList with the parameters at the command line
        /// </summary>
        /// <param name="parameterList">Parameter list</param>
        /// <returns>True if any of the parameters are not present in parameterList()</returns>
        public bool InvalidParametersPresent(List<string> parameterList)
        {
            const bool caseSensitive = false;
            return InvalidParametersPresent(parameterList, caseSensitive);
        }

        // ReSharper disable once UnusedMember.Global
        /// <summary>
        /// Compares the parameter names in parameterList with the parameters at the command line
        /// </summary>
        /// <param name="parameterList">Parameter list</param>
        /// <returns>True if any of the parameters are not present in parameterList()</returns>
        public bool InvalidParametersPresent(string[] parameterList)
        {
            const bool caseSensitive = false;
            return InvalidParametersPresent(parameterList, caseSensitive);
        }

        /// <summary>
        /// Compares the parameter names in parameterList with the parameters at the command line
        /// </summary>
        /// <param name="parameterList">Parameter list</param>
        /// <param name="caseSensitive">True to perform case-sensitive matching of the parameter name</param>
        /// <returns>True if any of the parameters are not present in parameterList()</returns>
        public bool InvalidParametersPresent(string[] parameterList, bool caseSensitive)
        {
            return InvalidParameters(parameterList.ToList(), caseSensitive).Count > 0;
        }

        /// <summary>
        /// Validate that the user-provided parameters are in the validParameters list
        /// </summary>
        /// <param name="validParameters"></param>
        /// <param name="caseSensitive"></param>
        /// <returns></returns>
        public bool InvalidParametersPresent(List<string> validParameters, bool caseSensitive)
        {
            return InvalidParameters(validParameters, caseSensitive).Count > 0;
        }

        /// <summary>
        /// Retrieve a list of the user-provided parameters that are not in validParameters
        /// </summary>
        /// <param name="validParameters"></param>
        /// <returns></returns>
        public List<string> InvalidParameters(List<string> validParameters)
        {
            const bool caseSensitive = false;
            return InvalidParameters(validParameters, caseSensitive);
        }

        /// <summary>
        /// Retrieve a list of the user-provided parameters that are not in validParameters
        /// </summary>
        /// <param name="validParameters"></param>
        /// <param name="caseSensitive"></param>
        /// <returns></returns>
        public List<string> InvalidParameters(List<string> validParameters, bool caseSensitive)
        {
            var invalidParams = new List<string>();

            try
            {
                // Find items in mSwitches whose keys are not in validParameters)
                foreach (KeyValuePair<string, string> item in mSwitches)
                {
                    string itemKey = item.Key;
                    int intMatchCount = 0;

                    if (caseSensitive)
                    {
                        intMatchCount = (from validItem in validParameters where validItem == itemKey select validItem).Count();
                    }
                    else
                    {
                        intMatchCount = (from validItem in validParameters where string.Equals(validItem, itemKey, StringComparison.InvariantCultureIgnoreCase) select validItem).Count();
                    }

                    if (intMatchCount == 0)
                    {
                        invalidParams.Add(item.Key);
                    }
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Error in InvalidParameters", ex);
            }

            return invalidParams;
        }

        /// <summary>
        /// Look for parameter on the command line
        /// </summary>
        /// <param name="paramName">Parameter name</param>
        /// <returns>True if present, otherwise false</returns>
        /// <remarks>Does not work for /? or /help -- for those, use .NeedToShowHelp</remarks>
        public bool IsParameterPresent(string paramName)
        {
            string paramValue = string.Empty;
            const bool caseSensitive = false;
            return RetrieveValueForParameter(paramName, out paramValue, caseSensitive);
        }

        /// <summary>
        /// Parse the parameters and switches at the command line; uses / for the switch character and : for the switch parameter character
        /// </summary>
        /// <returns>Returns True if any command line parameters were found; otherwise false</returns>
        /// <remarks>If /? or /help is found, then returns False and sets mShowHelp to True</remarks>
        public bool ParseCommandLine()
        {
            return ParseCommandLine(DEFAULT_SWITCH_CHAR, DEFAULT_SWITCH_PARAM_CHAR);
        }

        // ReSharper disable once UnusedMember.Global
        /// <summary>
        /// Parse the parameters and switches at the command line; uses : for the switch parameter character
        /// </summary>
        /// <returns>Returns True if any command line parameters were found; otherwise false</returns>
        /// <remarks>If /? or /help is found, then returns False and sets mShowHelp to True</remarks>
        public bool ParseCommandLine(string switchStartChar)
        {
            return ParseCommandLine(switchStartChar, DEFAULT_SWITCH_PARAM_CHAR);
        }

        /// <summary>
        /// Parse the parameters and switches at the command line
        /// </summary>
        /// <param name="switchStartChar"></param>
        /// <param name="switchParameterChar"></param>
        /// <returns>Returns True if any command line parameters were found; otherwise false</returns>
        /// <remarks>If /? or /help is found, then returns False and sets mShowHelp to True</remarks>
        public bool ParseCommandLine(string switchStartChar, string switchParameterChar)
        {
            // Returns True if any command line parameters were found
            // Otherwise, returns false
            //
            // If /? or /help is found, then returns False and sets mShowHelp to True

            string strCmdLine = null;

            mSwitches.Clear();
            mNonSwitchParameters.Clear();

            try
            {
                try
                {
                    // .CommandLine() returns the full command line
                    strCmdLine = Environment.CommandLine;

                    // .GetCommandLineArgs splits the command line at spaces, though it keeps text between double quotes together
                    // Note that .NET will strip out the starting and ending double quote if the user provides a parameter like this:
                    // MyProgram.exe "C:\Program Files\FileToProcess"
                    //
                    // In this case, paramList(1) will not have a double quote at the start but it will have a double quote at the end:
                    //  paramList(1) = C:\Program Files\FileToProcess"

                    // One very odd feature of Environment.GetCommandLineArgs() is that if the command line looks like this:
                    //    MyProgram.exe "D:\My Folder\Subfolder\" /O:D:\OutputFolder
                    // Then paramList will have:
                    //    paramList(1) = D:\My Folder\Subfolder" /O:D:\OutputFolder
                    //
                    // To avoid this problem instead specify the command line as:
                    //    MyProgram.exe "D:\My Folder\Subfolder" /O:D:\OutputFolder
                    // which gives:
                    //    paramList(1) = D:\My Folder\Subfolder
                    //    paramList(2) = /O:D:\OutputFolder
                    //
                    // Due to the idiosyncrasies of .GetCommandLineArgs, we will instead use SplitCommandLineParams to do the splitting
                    // paramList = Environment.GetCommandLineArgs()
                }
                catch (Exception ex)
                {
                    // In .NET 1.x, programs would fail if called from a network share
                    // This appears to be fixed in .NET 2.0 and above
                    // If an exception does occur here, we'll show the error message at the console, then sleep for 2 seconds

                    Console.WriteLine("------------------------------------------------------------------------------");
                    Console.WriteLine("This program cannot be run from a network share.  Please map a drive to the");
                    Console.WriteLine(" network share you are currently accessing or copy the program files and");
                    Console.WriteLine(" required DLL's to your local computer.");
                    Console.WriteLine(" Exception: " + ex.Message);
                    Console.WriteLine("------------------------------------------------------------------------------");

                    PauseAtConsole(5000, 1000);

                    mShowHelp = true;
                    return false;
                }

                if (DebugMode)
                {
                    Console.WriteLine();
                    Console.WriteLine("Debugging command line parsing");
                    Console.WriteLine();
                }

                var paramList = SplitCommandLineParams(strCmdLine);

                if (DebugMode)
                {
                    Console.WriteLine();
                }

                if (string.IsNullOrWhiteSpace(strCmdLine))
                {
                    return false;
                }
                else if (strCmdLine.IndexOf(switchStartChar + "?", StringComparison.Ordinal) > 0 || strCmdLine.ToLower().IndexOf(switchStartChar + "help", StringComparison.Ordinal) > 0)
                {
                    mShowHelp = true;
                    return false;
                }

                // Parse the command line
                // Note that paramList(0) is the path to the Executable for the calling program

                for (var intIndex = 1; intIndex <= paramList.Length - 1; intIndex++)
                {
                    if (paramList[intIndex].Length > 0)
                    {
                        var paramName = paramList[intIndex].TrimStart(' ');
                        var paramValue = string.Empty;
                        bool isSwitchParam = false;

                        if (paramName.StartsWith(switchStartChar))
                        {
                            isSwitchParam = true;
                        }
                        else if (paramName.StartsWith(ALTERNATE_SWITCH_CHAR) || paramName.StartsWith(DEFAULT_SWITCH_CHAR))
                        {
                            isSwitchParam = true;
                        }
                        else
                        {
                            // Parameter doesn't start with switchStartChar or / or -
                            isSwitchParam = false;
                        }

                        if (isSwitchParam)
                        {
                            // Look for switchParameterChar in paramList(intIndex)
                            var charIndex = paramList[intIndex].IndexOf(switchParameterChar);

                            if (charIndex >= 0)
                            {
                                // Parameter is of the form /I:MyParam or /I:"My Parameter" or -I:"My Parameter" or /MyParam:Setting
                                paramValue = paramName.Substring(charIndex + 1).Trim();

                                // Remove any starting and ending quotation marks
                                paramValue = paramValue.Trim('"');

                                paramName = paramName.Substring(0, charIndex);
                            }
                            else
                            {
                                // Parameter is of the form /S or -S
                            }

                            // Remove the switch character from paramName
                            paramName = paramName.Substring(1).Trim();

                            if (DebugMode)
                            {
                                Console.WriteLine("SwitchParam: " + paramName + "=" + paramValue);
                            }

                            // Note: .Item() will add paramName if it doesn't exist (which is normally the case)
                            mSwitches[paramName] = paramValue;
                        }
                        else
                        {
                            // Non-switch parameter since switchParameterChar was not found and does not start with switchStartChar

                            // Remove any starting and ending quotation marks
                            paramName = paramName.Trim('"');

                            if (DebugMode)
                            {
                                Console.WriteLine("NonSwitchParam " + mNonSwitchParameters.Count + ": " + paramName);
                            }

                            mNonSwitchParameters.Add(paramName);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Error in ParseCommandLine", ex);
            }

            if (DebugMode)
            {
                Console.WriteLine();
                Console.WriteLine("Switch Count = " + mSwitches.Count);
                Console.WriteLine("NonSwitch Count = " + mNonSwitchParameters.Count);
                Console.WriteLine();
            }

            return mSwitches.Count + mNonSwitchParameters.Count > 0;
        }

        /// <summary>
        /// Pause the program for the specified number of milliseconds, displaying a period at a set interval while paused
        /// </summary>
        /// <param name="millisecondsToPause">Milliseconds to pause; default 5 seconds</param>
        /// <param name="millisecondsBetweenDots">Seconds between each period; default 1 second</param>
        public static void PauseAtConsole(int millisecondsToPause, int millisecondsBetweenDots)
        {
            int totalIterations = 0;

            Console.WriteLine();
            Console.Write("Continuing in " + (millisecondsToPause / 1000.0).ToString("0") + " seconds ");

            try
            {
                if (millisecondsBetweenDots == 0)
                    millisecondsBetweenDots = millisecondsToPause;

                totalIterations = Convert.ToInt32(Math.Round(millisecondsToPause / Convert.ToDouble(millisecondsBetweenDots), 0));
            }
            catch (Exception)
            {
                totalIterations = 1;
            }

            var iteration = 0;
            do
            {
                Console.Write('.');

                Thread.Sleep(millisecondsBetweenDots);

                iteration += 1;
            } while (iteration < totalIterations);

            Console.WriteLine();
        }

        // ReSharper disable once UnusedMember.Global
        /// <summary>
        /// Returns the value of the non-switch parameter at the given index
        /// </summary>
        /// <param name="parameterIndex">Parameter index</param>
        /// <returns>The value of the parameter at the given index; empty string if no value or invalid index</returns>
        public string RetrieveNonSwitchParameter(int parameterIndex)
        {
            string paramValue = string.Empty;

            if (parameterIndex < mNonSwitchParameters.Count)
            {
                paramValue = mNonSwitchParameters[parameterIndex];
            }

            if (paramValue == null)
            {
                return string.Empty;
            }

            return paramValue;
        }

        // ReSharper disable once UnusedMember.Global
        /// <summary>
        /// Returns the parameter at the given index
        /// </summary>
        /// <param name="parameterIndex">Parameter index</param>
        /// <param name="paramName">Parameter name (output)</param>
        /// <param name="paramValue">Value associated with the parameter; empty string if no value (output)</param>
        /// <returns></returns>
        public bool RetrieveParameter(int parameterIndex, out string paramName, out string paramValue)
        {
            try
            {
                paramName = string.Empty;
                paramValue = string.Empty;

                if (parameterIndex < mSwitches.Count)
                {
                    Dictionary<string, string>.Enumerator iEnum = mSwitches.GetEnumerator();

                    var switchIndex = 0;
                    while (iEnum.MoveNext())
                    {
                        if (switchIndex == parameterIndex)
                        {
                            paramName = iEnum.Current.Key;
                            paramValue = iEnum.Current.Value;
                            return true;
                        }
                        switchIndex += 1;
                    }
                }
                else
                {
                    return false;
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Error in RetrieveParameter", ex);
            }

            return false;
        }

        // ReSharper disable once UnusedMember.Global
        /// <summary>
        /// Look for parameter on the command line and returns its value in paramValue
        /// </summary>
        /// <param name="paramName">Parameter name</param>
        /// <param name="paramValue">Value associated with the parameter; empty string if no value (output)</param>
        /// <returns>True if present, otherwise false</returns>
        public bool RetrieveValueForParameter(string paramName, out string paramValue)
        {
            return RetrieveValueForParameter(paramName, out paramValue, false);
        }

        /// <summary>
        /// Look for parameter on the command line and returns its value in paramValue
        /// </summary>
        /// <param name="paramName">Parameter name</param>
        /// <param name="paramValue">Value associated with the parameter; empty string if no value (output)</param>
        /// <param name="caseSensitive">True to perform case-sensitive matching of the parameter name</param>
        /// <returns>True if present, otherwise false</returns>
        public bool RetrieveValueForParameter(string paramName, out string paramValue, bool caseSensitive)
        {
            try
            {
                paramValue = string.Empty;

                if (caseSensitive)
                {
                    if (mSwitches.ContainsKey(paramName))
                    {
                        paramValue = Convert.ToString(mSwitches[paramName]);
                        return true;
                    }
                    else
                    {
                        return false;
                    }
                }
                else
                {
                    var query = (from item in mSwitches where string.Equals(item.Key, paramName, StringComparison.InvariantCultureIgnoreCase) select item).ToList();

                    if (query.Count == 0)
                    {
                        return false;
                    }

                    paramValue = query.FirstOrDefault().Value;
                    return true;
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Error in RetrieveValueForParameter", ex);
            }
        }

        private string[] SplitCommandLineParams(string strCmdLine)
        {
            List<string> paramList = new List<string>();

            var indexStart = 0;
            var indexEnd = 0;

            try
            {
                if (!string.IsNullOrEmpty(strCmdLine))
                {
                    // Make sure the command line doesn't have any carriage return or linefeed characters
                    strCmdLine = strCmdLine.Replace("\r\n", " ");
                    strCmdLine = strCmdLine.Replace("\r", " ");
                    strCmdLine = strCmdLine.Replace("\n", " ");

                    var insideDoubleQuotes = false;

                    while (indexStart < strCmdLine.Length)
                    {
                        // Step through the characters to find the next space
                        // However, if we find a double quote, then stop checking for spaces

                        if (strCmdLine[indexEnd] == '"')
                        {
                            insideDoubleQuotes = !insideDoubleQuotes;
                        }

                        if (!insideDoubleQuotes || indexEnd == strCmdLine.Length - 1)
                        {
                            if (strCmdLine[indexEnd] == ' ' || indexEnd == strCmdLine.Length - 1)
                            {
                                // Found the end of a parameter
                                var paramName = strCmdLine.Substring(indexStart, indexEnd - indexStart + 1).TrimEnd(' ');

                                if (paramName.StartsWith("\""))
                                {
                                    paramName = paramName.Substring(1);
                                }

                                if (paramName.EndsWith("\""))
                                {
                                    paramName = paramName.Substring(0, paramName.Length - 1);
                                }

                                if (!string.IsNullOrEmpty(paramName))
                                {
                                    if (DebugMode)
                                    {
                                        Console.WriteLine("Param " + paramList.Count + ": " + paramName);
                                    }
                                    paramList.Add(paramName);
                                }

                                indexStart = indexEnd + 1;
                            }
                        }

                        indexEnd += 1;
                    }
                }
            }
            catch (Exception ex)
            {
                throw new Exception("Error in SplitCommandLineParams", ex);
            }

            return paramList.ToArray();
        }
    }
}
