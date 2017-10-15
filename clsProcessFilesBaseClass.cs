using System;
using System.Collections.Generic;
using System.IO;

namespace PeptideHitResultsProcRunner
{
    /// <summary>
    /// This class can be used as a base class for classes that process a file or files, and create
    /// new output files in an output folder.  Note that this class contains simple error codes that
    /// can be set from any derived classes.  The derived classes can also set their own local error codes
    /// </summary>
    /// <remarks>
    /// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
    /// Copyright 2005, Battelle Memorial Institute.  All Rights Reserved.
    /// Started November 9, 2003
    /// </remarks>
    public abstract class clsProcessFilesBaseClass : clsProcessFilesOrFoldersBase
    {
        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks></remarks>
        protected clsProcessFilesBaseClass()
        {
            mFileDate = "October 15, 2017";
            ErrorCode = eProcessFilesErrorCodes.NoError;
        }

        #region "Constants and Enums"
        public enum eProcessFilesErrorCodes
        {
            NoError = 0,
            InvalidInputFilePath = 1,
            InvalidOutputFolderPath = 2,
            ParameterFileNotFound = 4,
            InvalidParameterFile = 8,
            FilePathError = 16,
            LocalizedError = 32,
            UnspecifiedError = -1
        }

        //' Copy the following to any derived classes
        //'Public Enum eDerivedClassErrorCodes
        //'    NoError = 0
        //'    UnspecifiedError = -1
        //'End Enum
        #endregion

        #region "Classwide Variables"
        //'Private mLocalErrorCode As eDerivedClassErrorCodes

        //'Public ReadOnly Property LocalErrorCode() As eDerivedClassErrorCodes
        //'    Get
        //'        Return mLocalErrorCode
        //'    End Get
        //'End Property


        #endregion

        #region "Interface Functions"

        /// <summary>
        /// This option applies when processing files matched with a wildcard
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool IgnoreErrorsWhenUsingWildcardMatching { get; set; }

        /// <summary>
        /// Error code reflecting processing outcome
        /// </summary>
        public eProcessFilesErrorCodes ErrorCode { get; set; }

        #endregion

        protected override void CleanupPaths(ref string inputFileOrFolderPath, ref string outputFolderPath)
        {
            CleanupFilePaths(ref inputFileOrFolderPath, ref outputFolderPath);
        }

        protected bool CleanupFilePaths(ref string inputFilePath, ref string outputFolderPath)
        {
            // Returns True if success, False if failure

            try
            {
                // Make sure inputFilePath points to a valid file
                var inputfile = new FileInfo(inputFilePath);

                if (!inputfile.Exists)
                {
                    if (ShowMessages)
                    {
                        ShowErrorMessage("Input file not found: " + inputFilePath);
                    }
                    else
                    {
                        LogMessage("Input file not found: " + inputFilePath, eMessageTypeConstants.ErrorMsg);
                    }

                    ErrorCode = eProcessFilesErrorCodes.InvalidInputFilePath;
                    return false;
                }

                if (string.IsNullOrWhiteSpace(outputFolderPath))
                {
                    // Define outputFolderPath based on inputFilePath
                    outputFolderPath = inputfile.DirectoryName;
                }

                // Make sure outputFolderPath points to a folder
                var outputFolder = new DirectoryInfo(outputFolderPath);

                if (!outputFolder.Exists)
                {
                    // outputFolderPath points to a non-existent folder; attempt to create it
                    outputFolder.Create();
                }

                mOutputFolderPath = string.Copy(outputFolder.FullName);

                return true;
            }
            catch (Exception ex)
            {
                HandleException("Error cleaning up the file paths", ex);
                return false;
            }

        }

        protected bool CleanupInputFilePath(ref string inputFilePath)
        {
            // Returns True if success, False if failure

            try
            {
                // Make sure inputFilePath points to a valid file
                var inputfile = new FileInfo(inputFilePath);

                if (!inputfile.Exists)
                {
                    if (ShowMessages)
                    {
                        ShowErrorMessage("Input file not found: " + inputFilePath);
                    }
                    else
                    {
                        LogMessage("Input file not found: " + inputFilePath, eMessageTypeConstants.ErrorMsg);
                    }

                    ErrorCode = eProcessFilesErrorCodes.InvalidInputFilePath;
                    return false;
                }

                return true;
            }
            catch (Exception ex)
            {
                HandleException("Error cleaning up the file paths", ex);
                return false;
            }

        }

        protected string GetBaseClassErrorMessage()
        {
            // Returns String.Empty if no error

            string errorMessage;

            switch (ErrorCode)
            {
                case eProcessFilesErrorCodes.NoError:
                    errorMessage = string.Empty;
                    break;
                case eProcessFilesErrorCodes.InvalidInputFilePath:
                    errorMessage = "Invalid input file path";
                    break;
                case eProcessFilesErrorCodes.InvalidOutputFolderPath:
                    errorMessage = "Invalid output folder path";
                    break;
                case eProcessFilesErrorCodes.ParameterFileNotFound:
                    errorMessage = "Parameter file not found";
                    break;
                case eProcessFilesErrorCodes.InvalidParameterFile:
                    errorMessage = "Invalid parameter file";
                    break;
                case eProcessFilesErrorCodes.FilePathError:
                    errorMessage = "General file path error";
                    break;
                case eProcessFilesErrorCodes.LocalizedError:
                    errorMessage = "Localized error";
                    break;
                case eProcessFilesErrorCodes.UnspecifiedError:
                    errorMessage = "Unspecified error";
                    break;
                default:
                    // This shouldn't happen
                    errorMessage = "Unknown error state";
                    break;
            }

            return errorMessage;
        }

        public virtual string[] GetDefaultExtensionsToParse()
        {
            var extensionsToParse = new string[1];

            extensionsToParse[0] = ".*";

            return extensionsToParse;
        }

        public bool ProcessFilesWildcard(string inputFolderPath)
        {
            return ProcessFilesWildcard(inputFolderPath, string.Empty, string.Empty);
        }

        public bool ProcessFilesWildcard(string inputFilePath, string outputFolderPath)
        {
            return ProcessFilesWildcard(inputFilePath, outputFolderPath, string.Empty);
        }

        public bool ProcessFilesWildcard(string inputFilePath, string outputFolderPath, string parameterFilePath)
        {
            return ProcessFilesWildcard(inputFilePath, outputFolderPath, parameterFilePath, true);
        }

        public bool ProcessFilesWildcard(string inputFilePath, string outputFolderPath, string parameterFilePath, bool resetErrorCode)
        {
            // Returns True if success, False if failure

            AbortProcessing = false;
            var success = true;
            try
            {
                // Possibly reset the error code
                if (resetErrorCode)
                    ErrorCode = eProcessFilesErrorCodes.NoError;

                if (!string.IsNullOrWhiteSpace(outputFolderPath))
                {
                    // Update the cached output folder path
                    mOutputFolderPath = string.Copy(outputFolderPath);
                }

                // See if inputFilePath contains a wildcard (* or ?)
                if (inputFilePath != null && inputFilePath.Contains("*") | inputFilePath.Contains("?"))
                {
                    // Obtain a list of the matching  files

                    // Copy the path into cleanPath and replace any * or ? characters with _
                    var cleanPath = inputFilePath.Replace("*", "_");
                    cleanPath = cleanPath.Replace("?", "_");

                    var cleanFileInfo = new FileInfo(cleanPath);
                    string inputFolderPath;
                    if (cleanFileInfo.Directory.Exists)
                    {
                        inputFolderPath = cleanFileInfo.DirectoryName;
                    }
                    else
                    {
                        // Use the directory that has the .exe file
                        inputFolderPath = GetAppFolderPath();
                    }

                    var inputFolder = new DirectoryInfo(inputFolderPath);

                    // Remove any directory information from inputFilePath
                    inputFilePath = Path.GetFileName(inputFilePath);

                    var matchCount = 0;
                    foreach (var inputfile in inputFolder.GetFiles(inputFilePath))
                    {
                        matchCount += 1;

                        success = ProcessFile(inputfile.FullName, outputFolderPath, parameterFilePath, resetErrorCode);

                        if (AbortProcessing)
                        {
                            break;
                        }

                        if (!success && !IgnoreErrorsWhenUsingWildcardMatching)
                        {
                            break;
                        }
                        if (matchCount % 100 == 0)
                            Console.Write(".");
                    }

                    if (matchCount == 0)
                    {
                        if (ErrorCode == eProcessFilesErrorCodes.NoError)
                        {
                            if (ShowMessages)
                            {
                                ShowErrorMessage("No match was found for the input file path: " + inputFilePath);
                            }
                            else
                            {
                                LogMessage("No match was found for the input file path: " + inputFilePath, eMessageTypeConstants.ErrorMsg);
                            }
                        }
                    }
                    else
                    {
                        Console.WriteLine();
                    }
                }
                else
                {
                    success = ProcessFile(inputFilePath, outputFolderPath, parameterFilePath, resetErrorCode);
                }

                return success;
            }
            catch (Exception ex)
            {
                HandleException("Error in ProcessFilesWildcard", ex);
                return false;
            }

        }

        public bool ProcessFile(string inputFilePath)
        {
            return ProcessFile(inputFilePath, string.Empty, string.Empty);
        }

        public bool ProcessFile(string inputFilePath, string outputFolderPath)
        {
            return ProcessFile(inputFilePath, outputFolderPath, string.Empty);
        }

        public bool ProcessFile(string inputFilePath, string outputFolderPath, string parameterFilePath)
        {
            return ProcessFile(inputFilePath, outputFolderPath, parameterFilePath, true);
        }

        // Main function for processing a single file
        public abstract bool ProcessFile(string inputFilePath, string outputFolderPath, string parameterFilePath, bool resetErrorCode);

        public bool ProcessFilesAndRecurseFolders(string inputFolderPath)
        {
            return ProcessFilesAndRecurseFolders(inputFolderPath, string.Empty, string.Empty);
        }

        public bool ProcessFilesAndRecurseFolders(string inputFilePathOrFolder, string outputFolderName)
        {
            return ProcessFilesAndRecurseFolders(inputFilePathOrFolder, outputFolderName, string.Empty);
        }

        public bool ProcessFilesAndRecurseFolders(string inputFilePathOrFolder, string outputFolderName, string parameterFilePath)
        {
            return ProcessFilesAndRecurseFolders(inputFilePathOrFolder, outputFolderName, string.Empty, false, parameterFilePath);
        }

        public bool ProcessFilesAndRecurseFolders(string inputFilePathOrFolder, string outputFolderName, string parameterFilePath, string[] extensionsToParse)
        {
            return ProcessFilesAndRecurseFolders(inputFilePathOrFolder, outputFolderName, string.Empty, false, parameterFilePath, 0, extensionsToParse);
        }

        public bool ProcessFilesAndRecurseFolders(string inputFilePathOrFolder, string outputFolderName, string outputFolderAlternatePath, bool recreateFolderHierarchyInAlternatePath)
        {
            return ProcessFilesAndRecurseFolders(inputFilePathOrFolder, outputFolderName, outputFolderAlternatePath, recreateFolderHierarchyInAlternatePath, string.Empty);
        }

        public bool ProcessFilesAndRecurseFolders(string inputFilePathOrFolder, string outputFolderName, string outputFolderAlternatePath, bool recreateFolderHierarchyInAlternatePath, string parameterFilePath)
        {
            return ProcessFilesAndRecurseFolders(inputFilePathOrFolder, outputFolderName, outputFolderAlternatePath, recreateFolderHierarchyInAlternatePath, parameterFilePath, 0);
        }

        public bool ProcessFilesAndRecurseFolders(string inputFilePathOrFolder, string outputFolderName, string outputFolderAlternatePath, bool recreateFolderHierarchyInAlternatePath, string parameterFilePath, int recurseFoldersMaxLevels)
        {
            return ProcessFilesAndRecurseFolders(inputFilePathOrFolder, outputFolderName, outputFolderAlternatePath, recreateFolderHierarchyInAlternatePath, parameterFilePath, recurseFoldersMaxLevels, GetDefaultExtensionsToParse());
        }

        // Main function for processing files in a folder (and subfolders)
        public bool ProcessFilesAndRecurseFolders(string inputFilePathOrFolder, string outputFolderName, string outputFolderAlternatePath, bool recreateFolderHierarchyInAlternatePath, string parameterFilePath, int recurseFoldersMaxLevels, string[] extensionsToParse)
        {
            // Calls ProcessFiles for all files in inputFilePathOrFolder and below having an extension listed in extensionsToParse()
            // The extensions should be of the form ".TXT" or ".RAW" (i.e. a period then the extension)
            // If any of the extensions is "*" or ".*" then all files will be processed
            // If inputFilePathOrFolder contains a filename with a wildcard (* or ?), then that information will be
            //  used to filter the files that are processed
            // If recurseFoldersMaxLevels is <=0 then we recurse infinitely

            // Examine inputFilePathOrFolder to see if it contains a filename; if not, assume it points to a folder
            // First, see if it contains a * or ?
            try
            {
                string inputFolderPath;
                if (string.IsNullOrWhiteSpace(inputFilePathOrFolder))
                {
                    inputFolderPath = string.Empty;
                }
                else if (inputFilePathOrFolder.Contains("*") || inputFilePathOrFolder.Contains("?"))
                {
                    // Copy the path into cleanPath and replace any * or ? characters with _
                    var cleanPath = inputFilePathOrFolder.Replace("*", "_");
                    cleanPath = cleanPath.Replace("?", "_");

                    var inputFile = new FileInfo(cleanPath);
                    if (inputFile.Directory.Exists)
                    {
                        inputFolderPath = inputFile.DirectoryName;
                    }
                    else
                    {
                        // Use the directory that has the .exe file
                        inputFolderPath = GetAppFolderPath();
                    }

                    // Remove any directory information from inputFilePath
                    inputFilePathOrFolder = Path.GetFileName(inputFilePathOrFolder);
                }
                else
                {
                    var inputfolder = new DirectoryInfo(inputFilePathOrFolder);
                    if (inputfolder.Exists)
                    {
                        inputFolderPath = inputfolder.FullName;
                        inputFilePathOrFolder = "*";
                    }
                    else
                    {
                        if (inputfolder.Parent != null && inputfolder.Parent.Exists)
                        {
                            inputFolderPath = inputfolder.Parent.FullName;
                            inputFilePathOrFolder = Path.GetFileName(inputFilePathOrFolder);
                        }
                        else
                        {
                            // Unable to determine the input folder path
                            inputFolderPath = string.Empty;
                        }
                    }
                }

                if (string.IsNullOrWhiteSpace(inputFolderPath))
                {
                    ErrorCode = eProcessFilesErrorCodes.InvalidInputFilePath;
                    return false;
                }

                // Validate the output folder path
                if (!string.IsNullOrWhiteSpace(outputFolderAlternatePath))
                {
                    try
                    {
                        var outputFolder = new DirectoryInfo(outputFolderAlternatePath);
                        if (!outputFolder.Exists)
                            outputFolder.Create();
                    }
                    catch (Exception ex)
                    {
                        ErrorCode = eProcessFilesErrorCodes.InvalidOutputFolderPath;
                        ShowErrorMessage("Error validating the alternate output folder path in ProcessFilesAndRecurseFolders:" + ex.Message);
                        return false;
                    }
                }

                // Initialize some parameters
                AbortProcessing = false;
                var fileProcessCount = 0;
                var fileProcessFailCount = 0;

                // Call RecurseFoldersWork
                const int recursionLevel = 1;
                var success = RecurseFoldersWork(inputFolderPath, inputFilePathOrFolder, outputFolderName,
                                             parameterFilePath, outputFolderAlternatePath,
                                             recreateFolderHierarchyInAlternatePath, extensionsToParse,
                                             ref fileProcessCount, ref fileProcessFailCount,
                                             recursionLevel, recurseFoldersMaxLevels);

                return success;
            }
            catch (Exception ex)
            {
                HandleException("Error in ProcessFilesAndRecurseFolders", ex);
                return false;
            }

        }

        private bool RecurseFoldersWork(string inputFolderPath, string fileNameMatch, string outputFolderName,
            string parameterFilePath, string outputFolderAlternatePath,
            bool recreateFolderHierarchyInAlternatePath, IList<string> extensionsToParse,
            ref int fileProcessCount, ref int fileProcessFailCount,
            int recursionLevel, int recurseFoldersMaxLevels)
        {
            // If recurseFoldersMaxLevels is <=0 then we recurse infinitely

            DirectoryInfo inputFolder;

            int extensionIndex;
            var processAllExtensions = false;

            string outputFolderPathToUse;
            bool success;

            try
            {
                inputFolder = new DirectoryInfo(inputFolderPath);
            }
            catch (Exception ex)
            {
                // Input folder path error
                HandleException("Error in RecurseFoldersWork", ex);
                ErrorCode = eProcessFilesErrorCodes.InvalidInputFilePath;
                return false;
            }

            try
            {
                if (!string.IsNullOrWhiteSpace(outputFolderAlternatePath))
                {
                    if (recreateFolderHierarchyInAlternatePath)
                    {
                        outputFolderAlternatePath = Path.Combine(outputFolderAlternatePath, inputFolder.Name);
                    }
                    outputFolderPathToUse = Path.Combine(outputFolderAlternatePath, outputFolderName);
                }
                else
                {
                    outputFolderPathToUse = outputFolderName;
                }
            }
            catch (Exception ex)
            {
                // Output file path error
                HandleException("Error in RecurseFoldersWork", ex);
                ErrorCode = eProcessFilesErrorCodes.InvalidOutputFolderPath;
                return false;
            }

            try
            {
                // Validate extensionsToParse()
                for (extensionIndex = 0; extensionIndex <= extensionsToParse.Count - 1; extensionIndex++)
                {
                    if (extensionsToParse[extensionIndex] == null)
                    {
                        extensionsToParse[extensionIndex] = string.Empty;
                    }
                    else
                    {
                        if (!extensionsToParse[extensionIndex].StartsWith("."))
                        {
                            extensionsToParse[extensionIndex] = "." + extensionsToParse[extensionIndex];
                        }

                        if (extensionsToParse[extensionIndex] == ".*")
                        {
                            processAllExtensions = true;
                            break;
                        }

                        extensionsToParse[extensionIndex] = extensionsToParse[extensionIndex].ToUpper();
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error in RecurseFoldersWork", ex);
                ErrorCode = eProcessFilesErrorCodes.UnspecifiedError;
                return false;
            }

            try
            {
                if (!string.IsNullOrWhiteSpace(outputFolderPathToUse))
                {
                    // Update the cached output folder path
                    mOutputFolderPath = string.Copy(outputFolderPathToUse);
                }

                ShowMessage("Examining " + inputFolderPath);

                // Process any matching files in this folder
                success = true;

                foreach (var inputfile in inputFolder.GetFiles(fileNameMatch))
                {
                    for (extensionIndex = 0; extensionIndex <= extensionsToParse.Count - 1; extensionIndex++)
                    {
                        if (processAllExtensions || inputfile.Extension.ToUpper() == extensionsToParse[extensionIndex])
                        {
                            success = ProcessFile(inputfile.FullName, outputFolderPathToUse, parameterFilePath, true);
                            if (!success)
                            {
                                fileProcessFailCount += 1;
                                success = true;
                            }
                            else
                            {
                                fileProcessCount += 1;
                            }
                            break;
                        }

                        if (AbortProcessing)
                            break;
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error in RecurseFoldersWork", ex);
                ErrorCode = eProcessFilesErrorCodes.InvalidInputFilePath;
                return false;
            }

            if (!AbortProcessing)
            {
                // If recurseFoldersMaxLevels is <=0 then we recurse infinitely
                //  otherwise, compare recursionLevel to recurseFoldersMaxLevels
                if (recurseFoldersMaxLevels <= 0 || recursionLevel <= recurseFoldersMaxLevels)
                {
                    // Call this function for each of the subfolders of inputFolder
                    foreach (var subFolder in inputFolder.GetDirectories())
                    {
                        success = RecurseFoldersWork(subFolder.FullName, fileNameMatch, outputFolderName,
                                                        parameterFilePath, outputFolderAlternatePath,
                                                        recreateFolderHierarchyInAlternatePath, extensionsToParse,
                                                        ref fileProcessCount, ref fileProcessFailCount,
                                                        recursionLevel + 1, recurseFoldersMaxLevels);

                        if (!success)
                            break;
                    }
                }
            }

            return success;
        }

        protected void SetBaseClassErrorCode(eProcessFilesErrorCodes eNewErrorCode)
        {
            ErrorCode = eNewErrorCode;
        }

        // The following functions should be placed in any derived class
        // Cannot define as abstract since it contains a customized enumerated type (eDerivedClassErrorCodes) in the function declaration

        // private void SetLocalErrorCode(eDerivedClassErrorCodes eNewErrorCode)
        // {
        //     SetLocalErrorCode(eNewErrorCode, false);
        // }
        //
        // private void SetLocalErrorCode(eDerivedClassErrorCodes eNewErrorCode, bool leaveExistingErrorCodeUnchanged)
        // {
        //     if (leaveExistingErrorCodeUnchanged && mLocalErrorCode != eDerivedClassErrorCodes.NoError)
        //     {
        //         // An error code is already defined; do not change it
        //     }
        //     else
        //     {
        //         mLocalErrorCode = eNewErrorCode;
        //
        //         if (eNewErrorCode == eDerivedClassErrorCodes.NoError)
        //         {
        //             if (base.ErrorCode == clsProcessFilesBaseClass.eProcessFilesErrorCodes.LocalizedError)
        //             {
        //                 base.SetBaseClassErrorCode(clsProcessFilesBaseClass.eProcessFilesErrorCodes.NoError);
        //             }
        //         }
        //         else
        //         {
        //             base.SetBaseClassErrorCode(clsProcessFilesBaseClass.eProcessFilesErrorCodes.LocalizedError);
        //         }
        //     }
        // }
    }
}
