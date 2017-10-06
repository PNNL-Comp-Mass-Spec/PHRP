using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text.RegularExpressions;
using System.Threading;

namespace PeptideHitResultsProcRunner
{
    /// <summary>
    /// This class contains functions used by both clsProcessFilesBaseClass and clsProcessFoldersBaseClass
    /// </summary>
    /// <remarks>
    /// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
    /// Created in October 2013
    /// Last updated in January 2017
    /// </remarks>
    public abstract class clsProcessFilesOrFoldersBase
    {
        #region "Constants and Enums"

        protected enum eMessageTypeConstants
        {
            Normal = 0,
            ErrorMsg = 1,
            Warning = 2
        }

        #endregion

        #region "Classwide Variables"
        protected bool mShowMessages = true;

        protected string mFileDate;
        protected bool mAbortProcessing;

        protected bool mLogMessagesToFile;
        protected bool mLogFileUsesDateStamp = true;
        protected string mLogFilePath;
        protected StreamWriter mLogFile;

        // This variable is updated when CleanupFilePaths() is called
        protected string mOutputFolderPath;
        protected string mLogFolderPath;            // If blank, then mOutputFolderPath will be used; if mOutputFolderPath is also blank, then the log is created in the same folder as the executing assembly

        public event ProgressResetEventHandler ProgressReset;
        public delegate void ProgressResetEventHandler();
        public event ProgressChangedEventHandler ProgressChanged;
        public delegate void ProgressChangedEventHandler(string taskDescription, float percentComplete);       // PercentComplete ranges from 0 to 100, but can contain decimal percentage values
        public event ProgressCompleteEventHandler ProgressComplete;
        public delegate void ProgressCompleteEventHandler();

        public event ErrorEventEventHandler ErrorEvent;
        public delegate void ErrorEventEventHandler(string strMessage);
        public event WarningEventEventHandler WarningEvent;
        public delegate void WarningEventEventHandler(string strMessage);
        public event MessageEventEventHandler MessageEvent;
        public delegate void MessageEventEventHandler(string strMessage);

        protected string mProgressStepDescription;
        protected float mProgressPercentComplete;          // Ranges from 0 to 100, but can contain decimal percentage values

        /// <summary>
        /// Keys in this dictionary are the log type and message (separated by an underscore), values are the most recent time the string was logged
        /// </summary>
        /// <remarks></remarks>
        private readonly Dictionary<string, DateTime> mLogDataCache;

        private const int MAX_LOGDATA_CACHE_SIZE = 100000;

        #endregion

        #region "Interface Functions"

        public bool AbortProcessing
        {
            get { return mAbortProcessing; }
            set { mAbortProcessing = value; }
        }

        public string FileVersion
        {
            get { return GetVersionForExecutingAssembly(); }
        }

        public string FileDate
        {
            get { return mFileDate; }
        }

        public string LogFilePath
        {
            get { return mLogFilePath; }
            set
            {
                if (value == null)
                    value = string.Empty;
                mLogFilePath = value;
            }
        }

        public string LogFolderPath
        {
            get { return mLogFolderPath; }
            set { mLogFolderPath = value; }
        }

        public bool LogMessagesToFile
        {
            get { return mLogMessagesToFile; }
            set { mLogMessagesToFile = value; }
        }

        public virtual string ProgressStepDescription
        {
            get { return mProgressStepDescription; }
        }

        // ProgressPercentComplete ranges from 0 to 100, but can contain decimal percentage values
        public float ProgressPercentComplete
        {
            get { return Convert.ToSingle(Math.Round(mProgressPercentComplete, 2)); }
        }

        public bool ShowMessages
        {
            get { return mShowMessages; }
            set { mShowMessages = value; }
        }

        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        /// <remarks></remarks>
        protected clsProcessFilesOrFoldersBase()
        {
            mProgressStepDescription = string.Empty;

            mOutputFolderPath = string.Empty;
            mLogFolderPath = string.Empty;
            mLogFilePath = string.Empty;

            mLogDataCache = new Dictionary<string, DateTime>();
        }

        public virtual void AbortProcessingNow()
        {
            mAbortProcessing = true;
        }

        protected abstract void CleanupPaths(ref string strInputFileOrFolderPath, ref string strOutputFolderPath);

        public void CloseLogFileNow()
        {
            if ((mLogFile != null))
            {
                mLogFile.Close();
                mLogFile = null;

                GarbageCollectNow();
                Thread.Sleep(100);
            }
        }

        /// <summary>
        /// Verifies that the specified .XML settings file exists in the user's local settings folder
        /// </summary>
        /// <param name="strApplicationName">Application name</param>
        /// <param name="strSettingsFileName">Settings file name</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static bool CreateSettingsFileIfMissing(string strApplicationName, string strSettingsFileName)
        {
            string strSettingsFilePathLocal = GetSettingsFilePathLocal(strApplicationName, strSettingsFileName);

            return CreateSettingsFileIfMissing(strSettingsFilePathLocal);
        }

        /// <summary>
        /// Verifies that the specified .XML settings file exists in the user's local settings folder
        /// </summary>
        /// <param name="strSettingsFilePathLocal">Full path to the local settings file, for example C:\Users\username\AppData\Roaming\AppName\SettingsFileName.xml</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static bool CreateSettingsFileIfMissing(string strSettingsFilePathLocal)
        {
            try
            {
                if (!File.Exists(strSettingsFilePathLocal))
                {
                    FileInfo fiMasterSettingsFile = default(FileInfo);
                    fiMasterSettingsFile = new FileInfo(Path.Combine(GetAppFolderPath(), Path.GetFileName(strSettingsFilePathLocal)));

                    if (fiMasterSettingsFile.Exists)
                    {
                        fiMasterSettingsFile.CopyTo(strSettingsFilePathLocal);
                    }
                }
            }
            catch (Exception)
            {
                // Ignore errors, but return false
                return false;
            }

            return true;
        }

        /// <summary>
        /// Perform garbage collection
        /// </summary>
        /// <remarks></remarks>
        public static void GarbageCollectNow()
        {
            const int intMaxWaitTimeMSec = 1000;
            GarbageCollectNow(intMaxWaitTimeMSec);
        }

        /// <summary>
        /// Perform garbage collection
        /// </summary>
        /// <param name="intMaxWaitTimeMSec"></param>
        /// <remarks></remarks>
        public static void GarbageCollectNow(int intMaxWaitTimeMSec)
        {
            const int THREAD_SLEEP_TIME_MSEC = 100;

            int intTotalThreadWaitTimeMsec = 0;
            if (intMaxWaitTimeMSec < 100)
                intMaxWaitTimeMSec = 100;
            if (intMaxWaitTimeMSec > 5000)
                intMaxWaitTimeMSec = 5000;

            Thread.Sleep(100);

            try
            {
                Thread gcThread = new Thread(GarbageCollectWaitForGC);
                gcThread.Start();

                intTotalThreadWaitTimeMsec = 0;
                while (gcThread.IsAlive && intTotalThreadWaitTimeMsec < intMaxWaitTimeMSec)
                {
                    Thread.Sleep(THREAD_SLEEP_TIME_MSEC);
                    intTotalThreadWaitTimeMsec += THREAD_SLEEP_TIME_MSEC;
                }
                if (gcThread.IsAlive)
                    gcThread.Abort();
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        protected static void GarbageCollectWaitForGC()
        {
            try
            {
                GC.Collect();
                GC.WaitForPendingFinalizers();
            }
            catch
            {
                // Ignore errors here
            }
        }

        /// <summary>
        /// Returns the full path to the folder into which this application should read/write settings file information
        /// </summary>
        /// <param name="strAppName"></param>
        /// <returns></returns>
        /// <remarks>For example, C:\Users\username\AppData\Roaming\AppName</remarks>
        public static string GetAppDataFolderPath(string strAppName)
        {
            string strAppDataFolder = null;

            if (string.IsNullOrWhiteSpace(strAppName))
            {
                strAppName = string.Empty;
            }

            try
            {
                strAppDataFolder = Path.Combine(Environment.GetFolderPath(Environment.SpecialFolder.ApplicationData), strAppName);
                if (!Directory.Exists(strAppDataFolder))
                {
                    Directory.CreateDirectory(strAppDataFolder);
                }
            }
            catch (Exception)
            {
                // Error creating the folder, revert to using the system Temp folder
                strAppDataFolder = Path.GetTempPath();
            }

            return strAppDataFolder;
        }

        /// <summary>
        /// Returns the full path to the folder that contains the currently executing .Exe or .Dll
        /// </summary>
        /// <returns></returns>
        /// <remarks></remarks>
        public static string GetAppFolderPath()
        {
            // Could use Application.StartupPath, but .GetExecutingAssembly is better
            return Path.GetDirectoryName(GetAppPath());
        }

        /// <summary>
        /// Returns the full path to the executing .Exe or .Dll
        /// </summary>
        /// <returns>File path</returns>
        /// <remarks></remarks>
        public static string GetAppPath()
        {
            return Assembly.GetExecutingAssembly().Location;
        }

        /// <summary>
        /// Returns the .NET assembly version followed by the program date
        /// </summary>
        /// <param name="strProgramDate"></param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static string GetAppVersion(string strProgramDate)
        {
            return Assembly.GetExecutingAssembly().GetName().Version.ToString() + " (" + strProgramDate + ")";
        }

        public abstract string GetErrorMessage();

        private string GetVersionForExecutingAssembly()
        {
            string strVersion = null;

            try
            {
                strVersion = Assembly.GetExecutingAssembly().GetName().Version.ToString();
            }
            catch (Exception)
            {
                strVersion = "??.??.??.??";
            }

            return strVersion;
        }

        /// <summary>
        /// Returns the full path to this application's local settings file
        /// </summary>
        /// <param name="strApplicationName"></param>
        /// <param name="strSettingsFileName"></param>
        /// <returns></returns>
        /// <remarks>For example, C:\Users\username\AppData\Roaming\AppName\SettingsFileName.xml</remarks>
        public static string GetSettingsFilePathLocal(string strApplicationName, string strSettingsFileName)
        {
            return Path.Combine(GetAppDataFolderPath(strApplicationName), strSettingsFileName);
        }

        protected void HandleException(string strBaseMessage, Exception ex)
        {
            if (string.IsNullOrWhiteSpace(strBaseMessage))
            {
                strBaseMessage = "Error";
            }

            if (ShowMessages)
            {
                // Note that ShowErrorMessage() will call LogMessage()
                ShowErrorMessage(strBaseMessage + ": " + ex.Message, true);
            }
            else
            {
                LogMessage(strBaseMessage + ": " + ex.Message, eMessageTypeConstants.ErrorMsg);
                throw new Exception(strBaseMessage, ex);
            }
        }

        protected void LogMessage(string strMessage)
        {
            LogMessage(strMessage, eMessageTypeConstants.Normal);
        }

        protected void LogMessage(string strMessage, eMessageTypeConstants eMessageType)
        {
            LogMessage(strMessage, eMessageType, intDuplicateHoldoffHours: 0);
        }

        protected void LogMessage(string strMessage, eMessageTypeConstants eMessageType, int intDuplicateHoldoffHours)
        {
            // Note that CleanupPaths() will update mOutputFolderPath, which is used here if mLogFolderPath is blank
            // Thus, be sure to call CleanupPaths (or update mLogFolderPath) before the first call to LogMessage

            string strMessageType = null;
            bool blnOpeningExistingFile = false;

            switch (eMessageType)
            {
                case eMessageTypeConstants.Normal:
                    strMessageType = "Normal";
                    break;
                case eMessageTypeConstants.ErrorMsg:
                    strMessageType = "Error";
                    break;
                case eMessageTypeConstants.Warning:
                    strMessageType = "Warning";
                    break;
                default:
                    strMessageType = "Unknown";
                    break;
            }

            if (mLogFile == null && mLogMessagesToFile)
            {
                try
                {
                    if (string.IsNullOrWhiteSpace(mLogFilePath))
                    {
                        // Auto-name the log file
                        mLogFilePath = Path.GetFileNameWithoutExtension(GetAppPath());
                        mLogFilePath += "_log";

                        if (mLogFileUsesDateStamp)
                        {
                            mLogFilePath += "_" + DateTime.Now.ToString("yyyy-MM-dd") + ".txt";
                        }
                        else
                        {
                            mLogFilePath += ".txt";
                        }
                    }

                    try
                    {
                        if (mLogFolderPath == null)
                            mLogFolderPath = string.Empty;

                        if (string.IsNullOrWhiteSpace(mLogFolderPath))
                        {
                            // Log folder is undefined; use mOutputFolderPath if it is defined
                            if (!string.IsNullOrWhiteSpace(mOutputFolderPath))
                            {
                                mLogFolderPath = string.Copy(mOutputFolderPath);
                            }
                        }

                        if (mLogFolderPath.Length > 0)
                        {
                            // Create the log folder if it doesn't exist
                            if (!Directory.Exists(mLogFolderPath))
                            {
                                Directory.CreateDirectory(mLogFolderPath);
                            }
                        }
                    }
                    catch (Exception)
                    {
                        mLogFolderPath = string.Empty;
                    }

                    if (!Path.IsPathRooted(mLogFilePath) && mLogFolderPath.Length > 0)
                    {
                        mLogFilePath = Path.Combine(mLogFolderPath, mLogFilePath);
                    }

                    blnOpeningExistingFile = File.Exists(mLogFilePath);

                    if ((blnOpeningExistingFile & mLogDataCache.Count == 0))
                    {
                        UpdateLogDataCache(mLogFilePath, DateTime.UtcNow.AddHours(-intDuplicateHoldoffHours));
                    }

                    mLogFile = new StreamWriter(new FileStream(mLogFilePath, FileMode.Append, FileAccess.Write, FileShare.Read));
                    mLogFile.AutoFlush = true;

                    if (!blnOpeningExistingFile)
                    {
                        mLogFile.WriteLine("Date\tType\tMessage");
                    }
                }
                catch (Exception ex)
                {
                    // Error creating the log file; set mLogMessagesToFile to false so we don't repeatedly try to create it
                    mLogMessagesToFile = false;
                    HandleException("Error opening log file", ex);
                    // Note: do not exit this function if an exception occurs
                }
            }

            if ((mLogFile != null))
            {
                var blnWriteToLog = true;

                string strLogKey = strMessageType + "_" + strMessage;
                DateTime dtLastLogTime = default(DateTime);
                bool blnMessageCached = false;

                if (mLogDataCache.TryGetValue(strLogKey, out dtLastLogTime))
                {
                    blnMessageCached = true;
                }
                else
                {
                    blnMessageCached = false;
                    dtLastLogTime = DateTime.UtcNow.AddHours(-(intDuplicateHoldoffHours + 1));
                }

                if (intDuplicateHoldoffHours > 0 && DateTime.UtcNow.Subtract(dtLastLogTime).TotalHours < intDuplicateHoldoffHours)
                {
                    blnWriteToLog = false;
                }

                if (blnWriteToLog)
                {
                    mLogFile.WriteLine(DateTime.Now.ToString("yyyy-MM-dd hh:mm:ss tt") + "\t" +
                                       strMessageType + "\t" +
                                       strMessage);

                    if (blnMessageCached)
                    {
                        mLogDataCache[strLogKey] = DateTime.UtcNow;
                    }
                    else
                    {
                        try
                        {
                            mLogDataCache.Add(strLogKey, DateTime.UtcNow);

                            if (mLogDataCache.Count > MAX_LOGDATA_CACHE_SIZE)
                            {
                                TrimLogDataCache();
                            }
                        }
                        catch (Exception)
                        {
                            // Ignore errors here
                        }
                    }
                }
            }

            RaiseMessageEvent(strMessage, eMessageType);
        }

        private string strLastMessage = "";
        private DateTime dtLastReportTime = DateTime.Now;

        private void RaiseMessageEvent(string strMessage, eMessageTypeConstants eMessageType)
        {
            if (!string.IsNullOrWhiteSpace(strMessage))
            {
                if (string.Equals(strMessage, strLastMessage) && DateTime.UtcNow.Subtract(dtLastReportTime).TotalSeconds < 0.5)
                {
                    // Duplicate message; do not raise any events
                }
                else
                {
                    dtLastReportTime = DateTime.UtcNow;
                    strLastMessage = string.Copy(strMessage);

                    switch (eMessageType)
                    {
                        case eMessageTypeConstants.Normal:
                            if (MessageEvent != null)
                            {
                                MessageEvent(strMessage);
                            }

                            break;
                        case eMessageTypeConstants.Warning:
                            if (WarningEvent != null)
                            {
                                WarningEvent(strMessage);
                            }

                            break;
                        case eMessageTypeConstants.ErrorMsg:
                            if (ErrorEvent != null)
                            {
                                ErrorEvent(strMessage);
                            }

                            break;
                        default:
                            if (MessageEvent != null)
                            {
                                MessageEvent(strMessage);
                            }

                            break;
                    }
                }
            }
        }

        protected void ResetProgress()
        {
            mProgressPercentComplete = 0;
            if (ProgressReset != null)
            {
                ProgressReset();
            }
        }

        protected void ResetProgress(string strProgressStepDescription)
        {
            UpdateProgress(strProgressStepDescription, 0);
            if (ProgressReset != null)
            {
                ProgressReset();
            }
        }

        protected void ShowErrorMessage(string strMessage)
        {
            ShowErrorMessage(strMessage, blnAllowLogToFile: true);
        }

        protected void ShowErrorMessage(string strMessage, bool blnAllowLogToFile)
        {
            ShowErrorMessage(strMessage, blnAllowLogToFile, intDuplicateHoldoffHours: 0);
        }

        protected void ShowErrorMessage(string strMessage, int intDuplicateHoldoffHours)
        {
            ShowErrorMessage(strMessage, blnAllowLogToFile: true, intDuplicateHoldoffHours: intDuplicateHoldoffHours);
        }

        protected void ShowErrorMessage(string strMessage, bool blnAllowLogToFile, int intDuplicateHoldoffHours)
        {
            const string strSeparator = "------------------------------------------------------------------------------";

            Console.WriteLine();
            Console.WriteLine(strSeparator);
            Console.WriteLine(strMessage);
            Console.WriteLine(strSeparator);
            Console.WriteLine();

            if (blnAllowLogToFile)
            {
                // Note that LogMessage will call RaiseMessageEvent
                LogMessage(strMessage, eMessageTypeConstants.ErrorMsg, intDuplicateHoldoffHours);
            }
            else
            {
                RaiseMessageEvent(strMessage, eMessageTypeConstants.ErrorMsg);
            }
        }

        protected void ShowMessage(string strMessage)
        {
            ShowMessage(strMessage, blnAllowLogToFile: true, blnPrecedeWithNewline: false, intDuplicateHoldoffHours: 0);
        }

        protected void ShowMessage(string strMessage, int intDuplicateHoldoffHours)
        {
            ShowMessage(strMessage, blnAllowLogToFile: true, blnPrecedeWithNewline: false, intDuplicateHoldoffHours: intDuplicateHoldoffHours);
        }

        protected void ShowMessage(string strMessage, bool blnAllowLogToFile)
        {
            ShowMessage(strMessage, blnAllowLogToFile, blnPrecedeWithNewline: false, intDuplicateHoldoffHours: 0);
        }

        protected void ShowMessage(string strMessage, bool blnAllowLogToFile, bool blnPrecedeWithNewline)
        {
            ShowMessage(strMessage, blnAllowLogToFile, blnPrecedeWithNewline, intDuplicateHoldoffHours: 0);
        }

        protected void ShowMessage(
            string strMessage,
            bool blnAllowLogToFile,
            bool blnPrecedeWithNewline,
            int intDuplicateHoldoffHours)
        {
            ShowMessage(strMessage, blnAllowLogToFile, blnPrecedeWithNewline, intDuplicateHoldoffHours, eMessageTypeConstants.Normal);
        }

        protected void ShowMessage(
            string strMessage,
            bool blnAllowLogToFile,
            bool blnPrecedeWithNewline,
            int intDuplicateHoldoffHours,
            eMessageTypeConstants eMessageType)
        {
            if (blnPrecedeWithNewline)
            {
                Console.WriteLine();
            }
            Console.WriteLine(strMessage);

            if (blnAllowLogToFile)
            {
                // Note that LogMessage will call RaiseMessageEvent
                LogMessage(strMessage, eMessageType, intDuplicateHoldoffHours);
            }
            else
            {
                RaiseMessageEvent(strMessage, eMessageType);
            }
        }

        protected void ShowWarning(string strMessage)
        {
            ShowMessage(strMessage, blnAllowLogToFile: true, blnPrecedeWithNewline: false, intDuplicateHoldoffHours: 0, eMessageType: eMessageTypeConstants.Warning);
        }

        protected void ShowWarning(string strMessage, int intDuplicateHoldoffHours)
        {
            ShowMessage(strMessage, blnAllowLogToFile: true, blnPrecedeWithNewline: false, intDuplicateHoldoffHours: intDuplicateHoldoffHours, eMessageType: eMessageTypeConstants.Warning);
        }

        protected void ShowWarning(string strMessage, bool blnAllowLogToFile)
        {
            ShowMessage(strMessage, blnAllowLogToFile, blnPrecedeWithNewline: false, intDuplicateHoldoffHours: 0, eMessageType: eMessageTypeConstants.Warning);
        }

        private void TrimLogDataCache()
        {
            if (mLogDataCache.Count < MAX_LOGDATA_CACHE_SIZE)
                return;

            try
            {
                // Remove entries from mLogDataCache so that the list count is 80% of MAX_LOGDATA_CACHE_SIZE

                // First construct a list of dates that we can sort to determine the datetime threshold for removal
                List<System.DateTime> lstDates = (from entry in mLogDataCache select entry.Value).ToList();

                // Sort by date
                lstDates.Sort();

                var intThresholdIndex = Convert.ToInt32(Math.Floor(mLogDataCache.Count - MAX_LOGDATA_CACHE_SIZE * 0.8));
                if (intThresholdIndex < 0)
                    intThresholdIndex = 0;

                var dtThreshold = lstDates[intThresholdIndex];

                // Construct a list of keys to be removed
                List<string> lstKeys = (from entry in mLogDataCache where entry.Value <= dtThreshold select entry.Key).ToList();

                // Remove each of the keys
                foreach (var strKey in lstKeys)
                {
                    mLogDataCache.Remove(strKey);
                }
            }
            catch (Exception)
            {
                // Ignore errors here
            }
        }

        private DateTime dtLastErrorShown = DateTime.MinValue;

        private void UpdateLogDataCache(string strLogFilePath, DateTime dtDateThresholdToStore)
        {
            var reParseLine = new Regex(@"^([^\t]+)\t([^\t]+)\t(.+)", RegexOptions.Compiled);

            try
            {
                mLogDataCache.Clear();

                using (var srLogFile = new StreamReader(new FileStream(strLogFilePath, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                {
                    while (!srLogFile.EndOfStream)
                    {
                        var strLineIn = srLogFile.ReadLine();
                        var reMatch = reParseLine.Match(strLineIn);

                        if (reMatch.Success)
                        {
                            DateTime dtLogTime = default(DateTime);
                            if (DateTime.TryParse(reMatch.Groups[1].Value, out dtLogTime))
                            {
                                dtLogTime = dtLogTime.ToUniversalTime();
                                if (dtLogTime >= dtDateThresholdToStore)
                                {
                                    string strKey = reMatch.Groups[2].Value + "_" + reMatch.Groups[3].Value;

                                    try
                                    {
                                        if (mLogDataCache.ContainsKey(strKey))
                                        {
                                            mLogDataCache[strKey] = dtLogTime;
                                        }
                                        else
                                        {
                                            mLogDataCache.Add(strKey, dtLogTime);
                                        }
                                    }
                                    catch (Exception)
                                    {
                                        // Ignore errors here
                                    }
                                }
                            }
                        }
                    }
                }

                if (mLogDataCache.Count > MAX_LOGDATA_CACHE_SIZE)
                {
                    TrimLogDataCache();
                }
            }
            catch (Exception ex)
            {
                if (DateTime.UtcNow.Subtract(dtLastErrorShown).TotalSeconds > 10)
                {
                    dtLastErrorShown = DateTime.UtcNow;
                    Console.WriteLine("Error caching the log file: " + ex.Message);
                }
            }
        }

        protected void UpdateProgress(string strProgressStepDescription)
        {
            UpdateProgress(strProgressStepDescription, mProgressPercentComplete);
        }

        protected void UpdateProgress(float sngPercentComplete)
        {
            UpdateProgress(ProgressStepDescription, sngPercentComplete);
        }

        protected void UpdateProgress(string strProgressStepDescription, float sngPercentComplete)
        {
            var blnDescriptionChanged = !string.Equals(strProgressStepDescription, mProgressStepDescription);

            mProgressStepDescription = string.Copy(strProgressStepDescription);
            if (sngPercentComplete < 0)
            {
                sngPercentComplete = 0;
            }
            else if (sngPercentComplete > 100)
            {
                sngPercentComplete = 100;
            }
            mProgressPercentComplete = sngPercentComplete;

            if (blnDescriptionChanged)
            {
                if (mProgressPercentComplete < float.Epsilon)
                {
                    LogMessage(mProgressStepDescription.Replace(Environment.NewLine, "; "));
                }
                else
                {
                    LogMessage(mProgressStepDescription + " (" + mProgressPercentComplete.ToString("0.0") + "% complete)".Replace(Environment.NewLine, "; "));
                }
            }

            if (ProgressChanged != null)
            {
                ProgressChanged(ProgressStepDescription, ProgressPercentComplete);
            }
        }

        protected void OperationComplete()
        {
            if (ProgressComplete != null)
            {
                ProgressComplete();
            }
        }
    }
}
