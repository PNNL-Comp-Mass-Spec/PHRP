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
            mFileDate = "October 17, 2013";
            mErrorCode = eProcessFilesErrorCodes.NoError;
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

        private eProcessFilesErrorCodes mErrorCode;

        protected bool mIgnoreErrorsWhenUsingWildcardMatching;

        #endregion

        #region "Interface Functions"

        /// <summary>
        /// This option applies when processing files matched with a wildcard
        /// </summary>
        /// <value></value>
        /// <returns></returns>
        /// <remarks></remarks>
        public bool IgnoreErrorsWhenUsingWildcardMatching
        {
            get => mIgnoreErrorsWhenUsingWildcardMatching;
            set => mIgnoreErrorsWhenUsingWildcardMatching = value;
        }

        public eProcessFilesErrorCodes ErrorCode => mErrorCode;

        #endregion

        protected override void CleanupPaths(ref string strInputFileOrFolderPath, ref string strOutputFolderPath)
        {
            CleanupFilePaths(ref strInputFileOrFolderPath, ref strOutputFolderPath);
        }

        protected bool CleanupFilePaths(ref string strInputFilePath, ref string strOutputFolderPath)
        {
            // Returns True if success, False if failure

            bool blnSuccess;

            try
            {
                // Make sure strInputFilePath points to a valid file
                var inputfile = new FileInfo(strInputFilePath);

                if (!inputfile.Exists)
                {
                    if (ShowMessages)
                    {
                        ShowErrorMessage("Input file not found: " + strInputFilePath);
                    }
                    else
                    {
                        LogMessage("Input file not found: " + strInputFilePath, eMessageTypeConstants.ErrorMsg);
                    }

                    mErrorCode = eProcessFilesErrorCodes.InvalidInputFilePath;
                    blnSuccess = false;
                }
                else
                {
                    if (string.IsNullOrWhiteSpace(strOutputFolderPath))
                    {
                        // Define strOutputFolderPath based on strInputFilePath
                        strOutputFolderPath = inputfile.DirectoryName;
                    }

                    // Make sure strOutputFolderPath points to a folder
                    var outputFolder = new DirectoryInfo(strOutputFolderPath);

                    if (!outputFolder.Exists)
                    {
                        // strOutputFolderPath points to a non-existent folder; attempt to create it
                        outputFolder.Create();
                    }

                    mOutputFolderPath = string.Copy(outputFolder.FullName);

                    blnSuccess = true;
                }
            }
            catch (Exception ex)
            {
                HandleException("Error cleaning up the file paths", ex);
                return false;
            }

            return blnSuccess;
        }

        protected bool CleanupInputFilePath(ref string strInputFilePath)
        {
            // Returns True if success, False if failure

            bool blnSuccess;

            try
            {
                // Make sure strInputFilePath points to a valid file
                var inputfile = new FileInfo(strInputFilePath);

                if (!inputfile.Exists)
                {
                    if (ShowMessages)
                    {
                        ShowErrorMessage("Input file not found: " + strInputFilePath);
                    }
                    else
                    {
                        LogMessage("Input file not found: " + strInputFilePath, eMessageTypeConstants.ErrorMsg);
                    }

                    mErrorCode = eProcessFilesErrorCodes.InvalidInputFilePath;
                    blnSuccess = false;
                }
                else
                {
                    blnSuccess = true;
                }
            }
            catch (Exception ex)
            {
                HandleException("Error cleaning up the file paths", ex);
                return false;
            }

            return blnSuccess;
        }

        protected string GetBaseClassErrorMessage()
        {
            // Returns String.Empty if no error

            string strErrorMessage;

            switch (ErrorCode)
            {
                case eProcessFilesErrorCodes.NoError:
                    strErrorMessage = string.Empty;
                    break;
                case eProcessFilesErrorCodes.InvalidInputFilePath:
                    strErrorMessage = "Invalid input file path";
                    break;
                case eProcessFilesErrorCodes.InvalidOutputFolderPath:
                    strErrorMessage = "Invalid output folder path";
                    break;
                case eProcessFilesErrorCodes.ParameterFileNotFound:
                    strErrorMessage = "Parameter file not found";
                    break;
                case eProcessFilesErrorCodes.InvalidParameterFile:
                    strErrorMessage = "Invalid parameter file";
                    break;
                case eProcessFilesErrorCodes.FilePathError:
                    strErrorMessage = "General file path error";
                    break;
                case eProcessFilesErrorCodes.LocalizedError:
                    strErrorMessage = "Localized error";
                    break;
                case eProcessFilesErrorCodes.UnspecifiedError:
                    strErrorMessage = "Unspecified error";
                    break;
                default:
                    // This shouldn't happen
                    strErrorMessage = "Unknown error state";
                    break;
            }

            return strErrorMessage;
        }

        public virtual string[] GetDefaultExtensionsToParse()
        {
            var strExtensionsToParse = new string[1];

            strExtensionsToParse[0] = ".*";

            return strExtensionsToParse;
        }

        public bool ProcessFilesWildcard(string strInputFolderPath)
        {
            return ProcessFilesWildcard(strInputFolderPath, string.Empty, string.Empty);
        }

        public bool ProcessFilesWildcard(string strInputFilePath, string strOutputFolderPath)
        {
            return ProcessFilesWildcard(strInputFilePath, strOutputFolderPath, string.Empty);
        }

        public bool ProcessFilesWildcard(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath)
        {
            return ProcessFilesWildcard(strInputFilePath, strOutputFolderPath, strParameterFilePath, true);
        }

        public bool ProcessFilesWildcard(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath, bool blnResetErrorCode)
        {
            // Returns True if success, False if failure

            mAbortProcessing = false;
            var blnSuccess = true;
            try
            {
                // Possibly reset the error code
                if (blnResetErrorCode)
                    mErrorCode = eProcessFilesErrorCodes.NoError;

                if (!string.IsNullOrWhiteSpace(strOutputFolderPath))
                {
                    // Update the cached output folder path
                    mOutputFolderPath = string.Copy(strOutputFolderPath);
                }

                // See if strInputFilePath contains a wildcard (* or ?)
                if (strInputFilePath != null && (strInputFilePath.Contains("*") | strInputFilePath.Contains("?")))
                {
                    // Obtain a list of the matching  files

                    // Copy the path into strCleanPath and replace any * or ? characters with _
                    var strCleanPath = strInputFilePath.Replace("*", "_");
                    strCleanPath = strCleanPath.Replace("?", "_");

                    var cleanFileInfo = new FileInfo(strCleanPath);
                    string strInputFolderPath;
                    if (cleanFileInfo.Directory.Exists)
                    {
                        strInputFolderPath = cleanFileInfo.DirectoryName;
                    }
                    else
                    {
                        // Use the directory that has the .exe file
                        strInputFolderPath = GetAppFolderPath();
                    }

                    var inputFolder = new DirectoryInfo(strInputFolderPath);

                    // Remove any directory information from strInputFilePath
                    strInputFilePath = Path.GetFileName(strInputFilePath);

                    var intMatchCount = 0;
                    foreach (var inputfile in inputFolder.GetFiles(strInputFilePath))
                    {
                        intMatchCount += 1;

                        blnSuccess = ProcessFile(inputfile.FullName, strOutputFolderPath, strParameterFilePath, blnResetErrorCode);

                        if (mAbortProcessing)
                        {
                            break;
                        }

                        if (!blnSuccess && !mIgnoreErrorsWhenUsingWildcardMatching)
                        {
                            break;
                        }
                        if (intMatchCount % 100 == 0)
                            Console.Write(".");
                    }

                    if (intMatchCount == 0)
                    {
                        if (mErrorCode == eProcessFilesErrorCodes.NoError)
                        {
                            if (ShowMessages)
                            {
                                ShowErrorMessage("No match was found for the input file path: " + strInputFilePath);
                            }
                            else
                            {
                                LogMessage("No match was found for the input file path: " + strInputFilePath, eMessageTypeConstants.ErrorMsg);
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
                    blnSuccess = ProcessFile(strInputFilePath, strOutputFolderPath, strParameterFilePath, blnResetErrorCode);
                }
            }
            catch (Exception ex)
            {
                HandleException("Error in ProcessFilesWildcard", ex);
                return false;
            }

            return blnSuccess;
        }

        public bool ProcessFile(string strInputFilePath)
        {
            return ProcessFile(strInputFilePath, string.Empty, string.Empty);
        }

        public bool ProcessFile(string strInputFilePath, string strOutputFolderPath)
        {
            return ProcessFile(strInputFilePath, strOutputFolderPath, string.Empty);
        }

        public bool ProcessFile(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath)
        {
            return ProcessFile(strInputFilePath, strOutputFolderPath, strParameterFilePath, true);
        }

        // Main function for processing a single file
        public abstract bool ProcessFile(string strInputFilePath, string strOutputFolderPath, string strParameterFilePath, bool blnResetErrorCode);

        public bool ProcessFilesAndRecurseFolders(string strInputFolderPath)
        {
            return ProcessFilesAndRecurseFolders(strInputFolderPath, string.Empty, string.Empty);
        }

        public bool ProcessFilesAndRecurseFolders(string strInputFilePathOrFolder, string strOutputFolderName)
        {
            return ProcessFilesAndRecurseFolders(strInputFilePathOrFolder, strOutputFolderName, string.Empty);
        }

        public bool ProcessFilesAndRecurseFolders(string strInputFilePathOrFolder, string strOutputFolderName, string strParameterFilePath)
        {
            return ProcessFilesAndRecurseFolders(strInputFilePathOrFolder, strOutputFolderName, string.Empty, false, strParameterFilePath);
        }

        public bool ProcessFilesAndRecurseFolders(string strInputFilePathOrFolder, string strOutputFolderName, string strParameterFilePath, string[] strExtensionsToParse)
        {
            return ProcessFilesAndRecurseFolders(strInputFilePathOrFolder, strOutputFolderName, string.Empty, false, strParameterFilePath, 0, strExtensionsToParse);
        }

        public bool ProcessFilesAndRecurseFolders(string strInputFilePathOrFolder, string strOutputFolderName, string strOutputFolderAlternatePath, bool blnRecreateFolderHierarchyInAlternatePath)
        {
            return ProcessFilesAndRecurseFolders(strInputFilePathOrFolder, strOutputFolderName, strOutputFolderAlternatePath, blnRecreateFolderHierarchyInAlternatePath, string.Empty);
        }

        public bool ProcessFilesAndRecurseFolders(string strInputFilePathOrFolder, string strOutputFolderName, string strOutputFolderAlternatePath, bool blnRecreateFolderHierarchyInAlternatePath, string strParameterFilePath)
        {
            return ProcessFilesAndRecurseFolders(strInputFilePathOrFolder, strOutputFolderName, strOutputFolderAlternatePath, blnRecreateFolderHierarchyInAlternatePath, strParameterFilePath, 0);
        }

        public bool ProcessFilesAndRecurseFolders(string strInputFilePathOrFolder, string strOutputFolderName, string strOutputFolderAlternatePath, bool blnRecreateFolderHierarchyInAlternatePath, string strParameterFilePath, int intRecurseFoldersMaxLevels)
        {
            return ProcessFilesAndRecurseFolders(strInputFilePathOrFolder, strOutputFolderName, strOutputFolderAlternatePath, blnRecreateFolderHierarchyInAlternatePath, strParameterFilePath, intRecurseFoldersMaxLevels, GetDefaultExtensionsToParse());
        }

        // Main function for processing files in a folder (and subfolders)
        public bool ProcessFilesAndRecurseFolders(string strInputFilePathOrFolder, string strOutputFolderName, string strOutputFolderAlternatePath, bool blnRecreateFolderHierarchyInAlternatePath, string strParameterFilePath, int intRecurseFoldersMaxLevels, string[] strExtensionsToParse)
        {
            // Calls ProcessFiles for all files in strInputFilePathOrFolder and below having an extension listed in strExtensionsToParse()
            // The extensions should be of the form ".TXT" or ".RAW" (i.e. a period then the extension)
            // If any of the extensions is "*" or ".*" then all files will be processed
            // If strInputFilePathOrFolder contains a filename with a wildcard (* or ?), then that information will be
            //  used to filter the files that are processed
            // If intRecurseFoldersMaxLevels is <=0 then we recurse infinitely

            bool blnSuccess;

            // Examine strInputFilePathOrFolder to see if it contains a filename; if not, assume it points to a folder
            // First, see if it contains a * or ?
            try
            {
                string strInputFolderPath;
                if (strInputFilePathOrFolder != null && (strInputFilePathOrFolder.Contains("*") | strInputFilePathOrFolder.Contains("?")))
                {
                    // Copy the path into strCleanPath and replace any * or ? characters with _
                    var strCleanPath = strInputFilePathOrFolder.Replace("*", "_");
                    strCleanPath = strCleanPath.Replace("?", "_");

                    var inputFile = new FileInfo(strCleanPath);
                    if (inputFile.Directory.Exists)
                    {
                        strInputFolderPath = inputFile.DirectoryName;
                    }
                    else
                    {
                        // Use the directory that has the .exe file
                        strInputFolderPath = GetAppFolderPath();
                    }

                    // Remove any directory information from strInputFilePath
                    strInputFilePathOrFolder = Path.GetFileName(strInputFilePathOrFolder);
                }
                else
                {
                    var inputfolder = new DirectoryInfo(strInputFilePathOrFolder);
                    if (inputfolder.Exists)
                    {
                        strInputFolderPath = inputfolder.FullName;
                        strInputFilePathOrFolder = "*";
                    }
                    else
                    {
                        if (inputfolder.Parent.Exists)
                        {
                            strInputFolderPath = inputfolder.Parent.FullName;
                            strInputFilePathOrFolder = Path.GetFileName(strInputFilePathOrFolder);
                        }
                        else
                        {
                            // Unable to determine the input folder path
                            strInputFolderPath = string.Empty;
                        }
                    }
                }

                if (!string.IsNullOrWhiteSpace(strInputFolderPath))
                {
                    // Validate the output folder path
                    if (!string.IsNullOrWhiteSpace(strOutputFolderAlternatePath))
                    {
                        try
                        {
                            var outputFolder = new DirectoryInfo(strOutputFolderAlternatePath);
                            if (!outputFolder.Exists)
                                outputFolder.Create();
                        }
                        catch (Exception ex)
                        {
                            mErrorCode = eProcessFilesErrorCodes.InvalidOutputFolderPath;
                            ShowErrorMessage("Error validating the alternate output folder path in ProcessFilesAndRecurseFolders:" + ex.Message);
                            return false;
                        }
                    }

                    // Initialize some parameters
                    mAbortProcessing = false;
                    var intFileProcessCount = 0;
                    var intFileProcessFailCount = 0;

                    // Call RecurseFoldersWork
                    const int intRecursionLevel = 1;
                    blnSuccess = RecurseFoldersWork(strInputFolderPath, strInputFilePathOrFolder, strOutputFolderName,
                                                    strParameterFilePath, strOutputFolderAlternatePath,
                                                    blnRecreateFolderHierarchyInAlternatePath, strExtensionsToParse,
                                                    ref intFileProcessCount, ref intFileProcessFailCount,
                                                    intRecursionLevel, intRecurseFoldersMaxLevels);
                }
                else
                {
                    mErrorCode = eProcessFilesErrorCodes.InvalidInputFilePath;
                    return false;
                }
            }
            catch (Exception ex)
            {
                HandleException("Error in ProcessFilesAndRecurseFolders", ex);
                return false;
            }

            return blnSuccess;
        }

        private bool RecurseFoldersWork(string strInputFolderPath, string strFileNameMatch, string strOutputFolderName,
            string strParameterFilePath, string strOutputFolderAlternatePath,
            bool blnRecreateFolderHierarchyInAlternatePath, IList<string> strExtensionsToParse,
            ref int intFileProcessCount, ref int intFileProcessFailCount,
            int intRecursionLevel, int intRecurseFoldersMaxLevels)
        {
            // If intRecurseFoldersMaxLevels is <=0 then we recurse infinitely

            DirectoryInfo inputFolder;

            int intExtensionIndex;
            var blnProcessAllExtensions = false;

            string strOutputFolderPathToUse;
            bool blnSuccess;

            try
            {
                inputFolder = new DirectoryInfo(strInputFolderPath);
            }
            catch (Exception ex)
            {
                // Input folder path error
                HandleException("Error in RecurseFoldersWork", ex);
                mErrorCode = eProcessFilesErrorCodes.InvalidInputFilePath;
                return false;
            }

            try
            {
                if (!string.IsNullOrWhiteSpace(strOutputFolderAlternatePath))
                {
                    if (blnRecreateFolderHierarchyInAlternatePath)
                    {
                        strOutputFolderAlternatePath = Path.Combine(strOutputFolderAlternatePath, inputFolder.Name);
                    }
                    strOutputFolderPathToUse = Path.Combine(strOutputFolderAlternatePath, strOutputFolderName);
                }
                else
                {
                    strOutputFolderPathToUse = strOutputFolderName;
                }
            }
            catch (Exception ex)
            {
                // Output file path error
                HandleException("Error in RecurseFoldersWork", ex);
                mErrorCode = eProcessFilesErrorCodes.InvalidOutputFolderPath;
                return false;
            }

            try
            {
                // Validate strExtensionsToParse()
                for (intExtensionIndex = 0; intExtensionIndex <= strExtensionsToParse.Count - 1; intExtensionIndex++)
                {
                    if (strExtensionsToParse[intExtensionIndex] == null)
                    {
                        strExtensionsToParse[intExtensionIndex] = string.Empty;
                    }
                    else
                    {
                        if (!strExtensionsToParse[intExtensionIndex].StartsWith("."))
                        {
                            strExtensionsToParse[intExtensionIndex] = "." + strExtensionsToParse[intExtensionIndex];
                        }

                        if (strExtensionsToParse[intExtensionIndex] == ".*")
                        {
                            blnProcessAllExtensions = true;
                            break;
                        }

                        strExtensionsToParse[intExtensionIndex] = strExtensionsToParse[intExtensionIndex].ToUpper();
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error in RecurseFoldersWork", ex);
                mErrorCode = eProcessFilesErrorCodes.UnspecifiedError;
                return false;
            }

            try
            {
                if (!string.IsNullOrWhiteSpace(strOutputFolderPathToUse))
                {
                    // Update the cached output folder path
                    mOutputFolderPath = string.Copy(strOutputFolderPathToUse);
                }

                ShowMessage("Examining " + strInputFolderPath);

                // Process any matching files in this folder
                blnSuccess = true;

                foreach (var inputfile in inputFolder.GetFiles(strFileNameMatch))
                {
                    for (intExtensionIndex = 0; intExtensionIndex <= strExtensionsToParse.Count - 1; intExtensionIndex++)
                    {
                        if (blnProcessAllExtensions || inputfile.Extension.ToUpper() == strExtensionsToParse[intExtensionIndex])
                        {
                            blnSuccess = ProcessFile(inputfile.FullName, strOutputFolderPathToUse, strParameterFilePath, true);
                            if (!blnSuccess)
                            {
                                intFileProcessFailCount += 1;
                                blnSuccess = true;
                            }
                            else
                            {
                                intFileProcessCount += 1;
                            }
                            break;
                        }

                        if (mAbortProcessing)
                            break;
                    }
                }
            }
            catch (Exception ex)
            {
                HandleException("Error in RecurseFoldersWork", ex);
                mErrorCode = eProcessFilesErrorCodes.InvalidInputFilePath;
                return false;
            }

            if (!mAbortProcessing)
            {
                // If intRecurseFoldersMaxLevels is <=0 then we recurse infinitely
                //  otherwise, compare intRecursionLevel to intRecurseFoldersMaxLevels
                if (intRecurseFoldersMaxLevels <= 0 || intRecursionLevel <= intRecurseFoldersMaxLevels)
                {
                    // Call this function for each of the subfolders of inputFolder
                    foreach (var subFolder in inputFolder.GetDirectories())
                    {
                        blnSuccess = RecurseFoldersWork(subFolder.FullName, strFileNameMatch, strOutputFolderName,
                                                        strParameterFilePath, strOutputFolderAlternatePath,
                                                        blnRecreateFolderHierarchyInAlternatePath, strExtensionsToParse,
                                                        ref intFileProcessCount, ref intFileProcessFailCount,
                                                        intRecursionLevel + 1, intRecurseFoldersMaxLevels);

                        if (!blnSuccess)
                            break;
                    }
                }
            }

            return blnSuccess;
        }

        protected void SetBaseClassErrorCode(eProcessFilesErrorCodes eNewErrorCode)
        {
            mErrorCode = eNewErrorCode;
        }

        // The following functions should be placed in any derived class
        // Cannot define as abstract since it contains a customized enumerated type (eDerivedClassErrorCodes) in the function declaration

        // private void SetLocalErrorCode(eDerivedClassErrorCodes eNewErrorCode)
        // {
        //     SetLocalErrorCode(eNewErrorCode, false);
        // }
        //
        // private void SetLocalErrorCode(eDerivedClassErrorCodes eNewErrorCode, bool blnLeaveExistingErrorCodeUnchanged)
        // {
        //     if (blnLeaveExistingErrorCodeUnchanged && mLocalErrorCode != eDerivedClassErrorCodes.NoError)
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
