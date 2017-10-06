using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Xml;

namespace PeptideHitResultsProcRunner
{
    /// <summary>
    /// This class is used to read or write settings in an Xml settings file
    /// Based on a class from the DMS Analysis Manager software written by Dave Clark and Gary Kiebel (PNNL, Richland, WA)
    /// Additional features added by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in October 2003
    /// Copyright 2005, Battelle Memorial Institute
    ///
    /// Updated in October 2004 to truly be case-insensitive if IsCaseSensitive = False when calling LoadSettings()
    /// Updated in August 2007 to remove the PRISM.Logging functionality and to include class XMLFileReader inside class XmlSettingsFileAccessor
    /// Updated in December 2010 to rename objects from Ini to XML
    /// </summary>
    public class XmlSettingsFileAccessor
    {
        /// <summary>
        /// Constructor
        /// </summary>
        public XmlSettingsFileAccessor()
        {
            mCaseSensitive = false;
            mSectionNames = new Dictionary<string, string>();

            mCachedSection = new udtRecentSectionType();
            mCachedSection.SectionName = string.Empty;
            mCachedSection.KeyNames = new Dictionary<string, string>();
        }

        private struct udtRecentSectionType
        {
            /// <summary>
            /// Stores the section name whose keys are cached; the section name is capitalized identically to that actually present in the Xml file
            /// </summary>
            public string SectionName;

            /// <summary>
            /// Keys for this section
            /// Keys in KeyNames are the lower case version of the name in the file if mCaseSensitive is true, or the actual version if mCaseSensitive is false
            /// Values in KeyNames are the actual way the key name is capitalized in the Xml file
            /// </summary>
            public Dictionary<string, string> KeyNames;
        }

        // XML file reader
        // Call LoadSettings to initialize, even if simply saving settings
        private string m_XMLFilePath = "";
        private XMLFileReader withEventsField_m_XMLFileAccessor;
        private XMLFileReader m_XMLFileAccessor
        {
            get { return withEventsField_m_XMLFileAccessor; }
            set
            {
                if (withEventsField_m_XMLFileAccessor != null)
                {
                    withEventsField_m_XMLFileAccessor.InformationMessage -= FileAccessorInfoMessageEvent;
                }
                withEventsField_m_XMLFileAccessor = value;
                if (withEventsField_m_XMLFileAccessor != null)
                {
                    withEventsField_m_XMLFileAccessor.InformationMessage += FileAccessorInfoMessageEvent;
                }
            }
        }

        private bool mCaseSensitive;

        // When mCaseSensitive = False, then mSectionNames stores the mapping between lowercase section name and actual section name stored in file
        //   If section is present more than once in file, then only grabs the first occurence of the section
        // When mCaseSensitive = True, then the mappings in mSectionNames are effectively not used
        private readonly Dictionary<string, string> mSectionNames;
        private udtRecentSectionType mCachedSection;
        public event InformationMessageEventHandler InformationMessage;

        public delegate void InformationMessageEventHandler(string msg);

        /// <summary>
        /// Loads the settings for the defined Xml Settings File.  Assumes names are not case sensitive
        /// </summary>
        /// <return>The function returns a boolean that shows if the file was successfully loaded.</return>
        public bool LoadSettings()
        {
            return LoadSettings(m_XMLFilePath, false);
        }

        /// <summary>
        /// Loads the settings for the defined Xml Settings File.   Assumes names are not case sensitive
        /// </summary>
        /// <param name="XmlSettingsFilePath">The path to the XML settings file.</param>
        /// <return>True if the file was successfully loaded (or created)</return>
        /// <remarks>The XML file will be created if it does not exist</remarks>
        public bool LoadSettings(string XmlSettingsFilePath)
        {
            return LoadSettings(XmlSettingsFilePath, false);
        }

        /// <summary>
        /// Loads the settings for the defined Xml Settings File
        /// </summary>
        /// <param name="XmlSettingsFilePath">The path to the XML settings file.</param>
        /// <param name="IsCaseSensitive">Case sensitive names if True.  Non-case sensitive if false.</param>
        /// <return>True if the file was successfully loaded (or created)</return>
        /// <remarks>The XML file will be created if it does not exist</remarks>
        public bool LoadSettings(string XmlSettingsFilePath, bool IsCaseSensitive)
        {
            mCaseSensitive = IsCaseSensitive;

            m_XMLFilePath = XmlSettingsFilePath;

            // Note: Always set IsCaseSensitive = True for XMLFileReader's constructor since this class handles
            //       case sensitivity mapping internally
            m_XMLFileAccessor = new XMLFileReader(m_XMLFilePath, true);
            if (m_XMLFileAccessor == null)
            {
                return false;
            }
            else if (m_XMLFileAccessor.Initialized)
            {
                CacheSectionNames();
                return true;
            }
            else
            {
                return false;
            }
        }

        public bool ManualParseXmlOrIniFile(string strFilePath)
        {
            m_XMLFilePath = strFilePath;

            // Note: Always set IsCaseSensitive = True for XMLFileReader's constructor since this class handles
            //       case sensitivity mapping internally
            m_XMLFileAccessor = new XMLFileReader(string.Empty, true);

            if (m_XMLFileAccessor == null)
            {
                return false;
            }
            else if (m_XMLFileAccessor.ManualParseXmlOrIniFile(strFilePath))
            {
                if (m_XMLFileAccessor.Initialized)
                {
                    CacheSectionNames();
                    return true;
                }
            }

            return false;
        }

        /// <summary>
        /// Saves the settings for the defined Xml Settings File.  Note that you must call LoadSettings to initialize the class prior to setting any values.
        /// </summary>
        /// <return>The function returns a boolean that shows if the file was successfully saved.</return>
        public bool SaveSettings()
        {
            if (m_XMLFileAccessor == null)
            {
                return false;
            }
            else if (m_XMLFileAccessor.Initialized)
            {
                m_XMLFileAccessor.OutputFilename = m_XMLFilePath;
                m_XMLFileAccessor.Save();
                return true;
            }
            else
            {
                return false;
            }
        }

        /// <summary>
        /// Checks if a section is present in the settings file.
        /// </summary>
        /// <param name="sectionName">The name of the section to look for.</param>
        /// <return>The function returns a boolean that shows if the section is present.</return>
        public bool SectionPresent(string sectionName)
        {
            var strSections = m_XMLFileAccessor.AllSections;

            foreach (var candidateSectionName in strSections)
            {
                if (SetNameCase(candidateSectionName) == SetNameCase(sectionName))
                    return true;
            }

            return false;
        }

        private bool CacheKeyNames(string sectionName)
        {
            // Looks up the Key Names for the given section, storing them in mCachedSection
            // This is done so that this class will know the correct capitalization for the key names

            List<string> strKeys = default(List<string>);

            // Lookup the correct capitalization for sectionName (only truly important if mCaseSensitive = False)
            var sectionNameInFile = GetCachedSectionName(sectionName);
            if (string.IsNullOrWhiteSpace(sectionNameInFile))
                return false;

            try
            {
                // Grab the keys for sectionName
                strKeys = m_XMLFileAccessor.AllKeysInSection(sectionNameInFile);
            }
            catch (Exception)
            {
                // Invalid section name; do not update anything
                return false;
            }

            if (strKeys == null)
            {
                return false;
            }

            // Update mCachedSection with the key names for the given section
            mCachedSection.SectionName = sectionNameInFile;
            mCachedSection.KeyNames.Clear();

            foreach (var keyName in strKeys)
            {
                // Change the key name to lowercase if mCaseSensitive is true
                var strKeyNameToStore = SetNameCase(keyName);

                if (!mCachedSection.KeyNames.Keys.Contains(strKeyNameToStore))
                {
                    mCachedSection.KeyNames.Add(strKeyNameToStore, keyName);
                }
                else
                {
                    Console.WriteLine("Note: ignoring duplicate key in the XML file: " + keyName);
                }
            }

            return true;
        }

        private void CacheSectionNames()
        {
            // Looks up the Section Names in the XML file
            // This is done so that this class will know the correct capitalization for the section names

            var strSections = m_XMLFileAccessor.AllSections;

            mSectionNames.Clear();

            foreach (var section in strSections)
            {
                var strSectionNameToStore = SetNameCase(section);

                if (!mSectionNames.ContainsKey(strSectionNameToStore))
                {
                    mSectionNames.Add(strSectionNameToStore, section);
                }
                else
                {
                    Console.WriteLine("Note: ignoring duplicate section in the XML file: " + section);
                }
            }
        }

        private string GetCachedKeyName(string sectionName, string keyName)
        {
            // Looks up the correct capitalization for key keyName in section sectionName
            // Returns String.Empty if not found

            bool blnSuccess = false;
            string sectionNameInFile = null;
            string keyNameToFind = null;

            // Lookup the correct capitalization for sectionName (only truly important if mCaseSensitive = False)
            sectionNameInFile = GetCachedSectionName(sectionName);
            if (string.IsNullOrWhiteSpace(sectionNameInFile))
                return string.Empty;

            if (mCachedSection.SectionName == sectionNameInFile)
            {
                blnSuccess = true;
            }
            else
            {
                // Update the keys for sectionName
                blnSuccess = CacheKeyNames(sectionName);
            }

            if (blnSuccess)
            {
                keyNameToFind = SetNameCase(keyName);
                if (mCachedSection.KeyNames.ContainsKey(keyNameToFind))
                {
                    return mCachedSection.KeyNames[keyNameToFind];
                }
                else
                {
                    return string.Empty;
                }
            }
            else
            {
                return string.Empty;
            }
        }

        private string GetCachedSectionName(string sectionName)
        {
            // Looks up the correct capitalization for sectionName
            // Returns String.Empty if not found

            var sectionNameToFind = SetNameCase(sectionName);

            if (mSectionNames.ContainsKey(sectionNameToFind))
            {
                return mSectionNames[sectionNameToFind];
            }
            else
            {
                return string.Empty;
            }
        }

        /// <summary>
        /// Return sectionName as is if mCaseSensitive is true, or return it lowercase if false
        /// </summary>
        /// <param name="sectionName">The name to be set.</param>
        /// <return>The function returns a string.</return>
        private string SetNameCase(string sectionName)
        {
            if ((mCaseSensitive))
            {
                return sectionName;
            }
            else
            {
                return sectionName.ToLower();
            }
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a String.</return>
        public string GetParam(string sectionName, string keyName, string valueIfMissing)
        {
            return GetParam(sectionName, keyName, valueIfMissing, out _);
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a String.</return>
        public string GetParam(string sectionName, string keyName, string valueIfMissing, out bool valueNotPresent)
        {
            string strResult = string.Empty;
            string sectionNameInFile = null;
            string keyNameInFile = null;
            bool blnValueFound = false;

            if (mCaseSensitive)
            {
                strResult = m_XMLFileAccessor.GetXMLValue(sectionName, keyName);
                if ((strResult != null))
                    blnValueFound = true;
            }
            else
            {
                sectionNameInFile = GetCachedSectionName(sectionName);
                if (sectionNameInFile.Length > 0)
                {
                    keyNameInFile = GetCachedKeyName(sectionName, keyName);
                    if (keyNameInFile.Length > 0)
                    {
                        strResult = m_XMLFileAccessor.GetXMLValue(sectionNameInFile, keyNameInFile);
                        if ((strResult != null))
                            blnValueFound = true;
                    }
                }
            }

            if (strResult == null || !blnValueFound)
            {
                valueNotPresent = true;
                return valueIfMissing;
            }
            else
            {
                valueNotPresent = false;
                return strResult;
            }
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns boolean True if the "value" attribute is "true". Otherwise, returns boolean False.</return>
        public bool GetParam(string sectionName, string keyName, bool valueIfMissing)
        {
            return GetParam(sectionName, keyName, valueIfMissing, out _);
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns boolean True if the "value" attribute is "true". Otherwise, returns boolean False.</return>
        public bool GetParam(string sectionName, string keyName, bool valueIfMissing, out bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), out blnNotFound);
            if (strResult == null || blnNotFound)
            {
                valueNotPresent = true;
                return valueIfMissing;
            }
            else
            {
                valueNotPresent = false;
                if (strResult.ToLower().Equals("true"))
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a Short.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public short GetParam(string sectionName, string keyName, short valueIfMissing)
        {
            return GetParam(sectionName, keyName, valueIfMissing, out _);
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a Short.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public short GetParam(string sectionName, string keyName, short valueIfMissing, out bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;
            short intValue = 0;

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), out blnNotFound);
            if (strResult == null || blnNotFound)
            {
                valueNotPresent = true;
                return valueIfMissing;
            }
            else
            {
                valueNotPresent = false;
                try
                {
                    if (short.TryParse(strResult, out intValue))
                    {
                        return intValue;
                    }
                    else if (strResult.ToLower().Equals("true"))
                    {
                        return -1;
                    }
                    else if (strResult.ToLower().Equals("false"))
                    {
                        return 0;
                    }
                    else
                    {
                        valueNotPresent = true;
                        return valueIfMissing;
                    }
                }
                catch (Exception)
                {
                    valueNotPresent = true;
                    return valueIfMissing;
                }
            }
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as an Integer.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public int GetParam(string sectionName, string keyName, int valueIfMissing)
        {
            return GetParam(sectionName, keyName, valueIfMissing, out _);
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as an Integer. If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public int GetParam(string sectionName, string keyName, int valueIfMissing, out bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;
            int intValue = 0;

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), out blnNotFound);
            if (strResult == null || blnNotFound)
            {
                valueNotPresent = true;
                return valueIfMissing;
            }
            else
            {
                valueNotPresent = false;
                try
                {
                    if (int.TryParse(strResult, out intValue))
                    {
                        return intValue;
                    }
                    else if (strResult.ToLower().Equals("true"))
                    {
                        return -1;
                    }
                    else if (strResult.ToLower().Equals("false"))
                    {
                        return 0;
                    }
                    else
                    {
                        valueNotPresent = true;
                        return valueIfMissing;
                    }
                }
                catch (Exception)
                {
                    valueNotPresent = true;
                    return valueIfMissing;
                }
            }
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a Long. If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public long GetParam(string sectionName, string keyName, long valueIfMissing)
        {
            return GetParam(sectionName, keyName, valueIfMissing, out _);
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a Long. If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public long GetParam(string sectionName, string keyName, long valueIfMissing, out bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;
            Int64 intValue = default(Int64);

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), out blnNotFound);
            if (strResult == null || blnNotFound)
            {
                valueNotPresent = true;
                return valueIfMissing;
            }
            else
            {
                valueNotPresent = false;
                try
                {
                    if (Int64.TryParse(strResult, out intValue))
                    {
                        return intValue;
                    }
                    else if (strResult.ToLower().Equals("true"))
                    {
                        return -1;
                    }
                    else if (strResult.ToLower().Equals("false"))
                    {
                        return 0;
                    }
                    else
                    {
                        valueNotPresent = true;
                        return valueIfMissing;
                    }
                }
                catch (Exception)
                {
                    valueNotPresent = true;
                    return valueIfMissing;
                }
            }
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a Single. If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public float GetParam(string sectionName, string keyName, float valueIfMissing)
        {
            return GetParam(sectionName, keyName, valueIfMissing, out _);
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a Single. If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public float GetParam(string sectionName, string keyName, float valueIfMissing, out bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;
            float sngValue = 0;

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), out blnNotFound);
            if (strResult == null || blnNotFound)
            {
                valueNotPresent = true;
                return valueIfMissing;
            }
            else
            {
                valueNotPresent = false;
                try
                {
                    if (float.TryParse(strResult, out sngValue))
                    {
                        return sngValue;
                    }
                    else if (strResult.ToLower().Equals("true"))
                    {
                        return -1;
                    }
                    else if (strResult.ToLower().Equals("false"))
                    {
                        return 0;
                    }
                    else
                    {
                        valueNotPresent = true;
                        return valueIfMissing;
                    }
                }
                catch (Exception)
                {
                    valueNotPresent = true;
                    return valueIfMissing;
                }
            }
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a Double.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public double GetParam(string sectionName, string keyName, double valueIfMissing)
        {
            return GetParam(sectionName, keyName, valueIfMissing);
        }

        /// <summary>
        /// Gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.  Returned ByRef.</param>
        /// <return>The function returns the name of the "value" attribute as a Double.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public double GetParam(string sectionName, string keyName, double valueIfMissing, out bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;
            double dblValue = 0;

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), out blnNotFound);
            if (strResult == null || blnNotFound)
            {
                valueNotPresent = true;
                return valueIfMissing;
            }
            else
            {
                valueNotPresent = false;
                try
                {
                    if (double.TryParse(strResult, out dblValue))
                    {
                        return dblValue;
                    }
                    else if (strResult.ToLower().Equals("true"))
                    {
                        return -1;
                    }
                    else if (strResult.ToLower().Equals("false"))
                    {
                        return 0;
                    }
                    else
                    {
                        valueNotPresent = true;
                        return valueIfMissing;
                    }
                }
                catch (Exception)
                {
                    valueNotPresent = true;
                    return valueIfMissing;
                }
            }
        }

        /// <summary>
        /// Legacy function name; calls SetXMLFilePath
        /// </summary>
        [Obsolete("Use SetXMLFilePath")]
        public void SetIniFilePath(string XmlSettingsFilePath)
        {
            SetXMLFilePath(XmlSettingsFilePath);
        }

        /// <summary>
        /// Sets the path to the Xml Settings File.
        /// </summary>
        /// <param name="XmlSettingsFilePath">The path to the XML settings file.</param>
        public void SetXMLFilePath(string XmlSettingsFilePath)
        {
            m_XMLFilePath = XmlSettingsFilePath;
        }

        /// <summary>
        /// Sets a new String value for the "value" attribute.
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="newValue">The new value for the "value".</param>
        /// <return>The function returns a boolean that shows if the change was done.</return>
        public bool SetParam(string sectionName, string keyName, string newValue)
        {
            string sectionNameInFile = null;
            string keyNameInFile = null;

            if (!mCaseSensitive)
            {
                sectionNameInFile = GetCachedSectionName(sectionName);
                if (sectionNameInFile.Length > 0)
                {
                    keyNameInFile = GetCachedKeyName(sectionName, keyName);
                    if (keyNameInFile.Length > 0)
                    {
                        // Section and Key are present; update them
                        return m_XMLFileAccessor.SetXMLValue(sectionNameInFile, keyNameInFile, newValue);
                    }
                    else
                    {
                        // Section is present, but the Key isn't; add teh key
                        return m_XMLFileAccessor.SetXMLValue(sectionNameInFile, keyName, newValue);
                    }
                }
            }

            // If we get here, then either mCaseSensitive = True or the section and key weren't found
            return m_XMLFileAccessor.SetXMLValue(sectionName, keyName, newValue);
        }

        /// <summary>
        /// Sets a new Boolean value for the "value" attribute.
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="newValue">The new value for the "value".</param>
        /// <return>The function returns a boolean that shows if the change was done.</return>
        public bool SetParam(string sectionName, string keyName, bool newValue)
        {
            return this.SetParam(sectionName, keyName, Convert.ToString(newValue));
        }

        /// <summary>
        /// Sets a new Short value for the "value" attribute.
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="newValue">The new value for the "value".</param>
        /// <return>The function returns a boolean that shows if the change was done.</return>
        public bool SetParam(string sectionName, string keyName, short newValue)
        {
            return this.SetParam(sectionName, keyName, Convert.ToString(newValue));
        }

        /// <summary>
        /// Sets a new Integer value for the "value" attribute.
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="newValue">The new value for the "value".</param>
        /// <return>The function returns a boolean that shows if the change was done.</return>
        public bool SetParam(string sectionName, string keyName, int newValue)
        {
            return this.SetParam(sectionName, keyName, Convert.ToString(newValue));
        }

        /// <summary>
        /// Sets a new Long value for the "value" attribute.
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="newValue">The new value for the "value".</param>
        /// <return>The function returns a boolean that shows if the change was done.</return>
        public bool SetParam(string sectionName, string keyName, long newValue)
        {
            return this.SetParam(sectionName, keyName, Convert.ToString(newValue));
        }

        /// <summary>
        /// Sets a new Single value for the "value" attribute.
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="newValue">The new value for the "value".</param>
        /// <return>The function returns a boolean that shows if the change was done.</return>
        public bool SetParam(string sectionName, string keyName, float newValue)
        {
            return this.SetParam(sectionName, keyName, Convert.ToString(newValue));
        }

        /// <summary>
        /// Sets a new Double value for the "value" attribute.
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="newValue">The new value for the "value".</param>
        /// <return>The function returns a boolean that shows if the change was done.</return>
        public bool SetParam(string sectionName, string keyName, double newValue)
        {
            return this.SetParam(sectionName, keyName, Convert.ToString(newValue));
        }

        /// <summary>
        /// Renames a section.
        /// </summary>
        /// <param name="sectionNameOld">The name of the old XML section name.</param>
        /// <param name="sectionNameNew">The new name for the XML section.</param>
        /// <return>The function returns a boolean that shows if the change was done.</return>
        public bool RenameSection(string sectionNameOld, string sectionNameNew)
        {
            string strSectionName = null;

            if (!mCaseSensitive)
            {
                strSectionName = GetCachedSectionName(sectionNameOld);
                if (strSectionName.Length > 0)
                {
                    return m_XMLFileAccessor.SetXMLSection(strSectionName, sectionNameNew);
                }
            }

            // If we get here, then either mCaseSensitive = True or the section wasn't found using GetCachedSectionName
            return m_XMLFileAccessor.SetXMLSection(sectionNameOld, sectionNameNew);
        }

        private void FileAccessorInfoMessageEvent(string msg)
        {
            if (InformationMessage != null)
            {
                InformationMessage(msg);
            }
        }

        /// <summary>
        /// Tools to manipulates XML Settings files.
        /// </summary>
        protected class XMLFileReader
        {
            public enum XMLItemTypeEnum
            {
                GetKeys = 0,
                GetValues = 1,
                GetKeysAndValues = 2
            }

            private string m_XmlFilename;
            private XmlDocument m_XmlDoc;

            private List<string> m_sections = new List<string>();
            private bool m_CaseSensitive = false;
            private string m_SaveFilename;
            private bool m_initialized = false;

            public bool NotifyOnEvent;
            public bool NotifyOnException;

            public event InformationMessageEventHandler InformationMessage;
            public delegate void InformationMessageEventHandler(string msg);

            /// <summary>
            /// Constructor: Initializes a new instance of the XMLFileReader (non case-sensitive)
            /// </summary>
            /// <param name="XmlFilename">The name of the XML file.</param>
            public XMLFileReader(string XmlFilename)
            {
                NotifyOnException = false;
                InitXMLFileReader(XmlFilename, false);
            }

            /// <summary>
            /// Constructor: Initializes a new instance of the XMLFileReader.
            /// </summary>
            /// <param name="XmlFilename">The name of the XML file.</param>
            /// <param name="IsCaseSensitive">Case sensitive as boolean.</param>
            /// <remarks>The XML file will be created if it does not exist</remarks>
            public XMLFileReader(string XmlFilename, bool IsCaseSensitive)
            {
                NotifyOnException = true;
                InitXMLFileReader(XmlFilename, IsCaseSensitive);
            }

            /// <summary>
            /// This routine is called by each of the constructors to make the actual assignments.
            /// </summary>
            /// <remarks>The XML file will be created if it does not exist</remarks>
            private void InitXMLFileReader(string strXmlFilename, bool IsCaseSensitive)
            {
                m_CaseSensitive = IsCaseSensitive;
                m_XmlDoc = new XmlDocument();

                if (string.IsNullOrEmpty(strXmlFilename))
                {
                    return;
                }

                // Try to load the file as an XML file
                try
                {
                    if (!File.Exists(strXmlFilename))
                    {
                        ManualParseXmlOrIniFile(strXmlFilename);
                        return;
                    }

                    m_XmlDoc.Load(strXmlFilename);
                    UpdateSections();
                    m_XmlFilename = strXmlFilename;
                    m_initialized = true;
                }
                catch
                {
                    // Exception occurred parsing XmlFilename
                    // Manually parse the file line-by-line
                    ManualParseXmlOrIniFile(strXmlFilename);
                }
            }

            /// <summary>
            /// Legacy property; calls XmlFilename
            /// </summary>
            public string IniFilename
            {
                get { return XmlFilename; }
            }

            /// <summary>
            /// This routine returns the name of the ini file.
            /// </summary>
            /// <return>The function returns the name of ini file.</return>
            public string XmlFilename
            {
                get
                {
                    if (!Initialized)
                    {
                        return string.Empty;
                    }
                    else
                    {
                        return (m_XmlFilename);
                    }
                }
            }

            /// <summary>
            /// Returns a boolean showing if the file was initialized or not.
            /// </summary>
            public bool Initialized
            {
                get { return m_initialized; }
            }

            /// <summary>
            /// Returns a boolean showing if the name is case sensitive or not.
            /// </summary>
            public bool CaseSensitive
            {
                get { return m_CaseSensitive; }
            }

            /// <summary>
            /// Return sectionName as is if CaseSensitive is true, or return it lowercase if false
            /// </summary>
            /// <param name="sectionName">The name to be set.</param>
            private string SetNameCase(string sectionName)
            {
                if ((CaseSensitive))
                {
                    return sectionName;
                }
                else
                {
                    return sectionName.ToLower();
                }
            }

            /// <summary>
            /// Returns the root element of the XML document
            /// </summary>
            private XmlElement GetRoot()
            {
                return m_XmlDoc.DocumentElement;
            }

            /// <summary>
            /// Gets the last section.
            /// </summary>
            /// <return>The function returns the last section as System.Xml.XmlElement.</return>
            private XmlElement GetLastSection()
            {
                if (m_sections.Count == 0)
                {
                    return GetRoot();
                }
                else
                {
                    return GetSection(m_sections.Last());
                }
            }

            /// <summary>
            /// Gets a section as System.Xml.XmlElement.
            /// </summary>
            /// <param name="sectionName">The name of a section.</param>
            /// <return>The function returns a section as System.Xml.XmlElement.</return>
            private XmlElement GetSection(string sectionName)
            {
                if (!string.IsNullOrWhiteSpace(sectionName))
                {
                    sectionName = SetNameCase(sectionName);
                    return (XmlElement) m_XmlDoc.SelectSingleNode("//section[@name='" + sectionName + "']");
                }
                return null;
            }

            /// <summary>
            /// Gets an item.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="keyName">The name of the key.</param>
            /// <return>The function returns a XML element.</return>
            private XmlElement GetItem(string sectionName, string keyName)
            {
                XmlElement section = default(XmlElement);
                if (!string.IsNullOrWhiteSpace(keyName))
                {
                    keyName = SetNameCase(keyName);
                    section = GetSection(sectionName);
                    if (((section != null)))
                    {
                        return (XmlElement) section.SelectSingleNode("item[@key='" + keyName + "']");
                    }
                }
                return null;
            }

            /// <summary>
            /// Legacy function name; calls SetXMLSection
            /// </summary>
            public bool SetIniSection(string oldSection, string newSection)
            {
                return SetXMLSection(oldSection, newSection);
            }

            /// <summary>
            /// Sets the ini section name.
            /// </summary>
            /// <param name="oldSection">The name of the old ini section name.</param>
            /// <param name="newSection">The new name for the ini section.</param>
            /// <return>The function returns a boolean that shows if the change was done.</return>
            public bool SetXMLSection(string oldSection, string newSection)
            {
                XmlElement section = default(XmlElement);
                if (!Initialized)
                {
                    throw new XMLFileReaderNotInitializedException();
                }
                if (!string.IsNullOrWhiteSpace(newSection))
                {
                    section = GetSection(oldSection);
                    if (((section != null)))
                    {
                        section.SetAttribute("name", SetNameCase(newSection));
                        UpdateSections();
                        return true;
                    }
                }
                return false;
            }

            /// <summary>
            /// Legacy function name; calls SetXMLValue
            /// </summary>
            public bool SetIniValue(string sectionName, string keyName, string newValue)
            {
                return SetXMLValue(sectionName, keyName, newValue);
            }

            /// <summary>
            /// Sets a new value for the "value" attribute.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="keyName">The name of the key.</param>
            /// <param name="newValue">The new value for the "value".</param>
            /// <return>The function returns a boolean that shows if the change was done.</return>
            public bool SetXMLValue(string sectionName, string keyName, string newValue)
            {
                XmlElement item = default(XmlElement);
                XmlElement section = default(XmlElement);
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                section = GetSection(sectionName);
                if (section == null)
                {
                    if (CreateSection(sectionName))
                    {
                        section = GetSection(sectionName);
                        // exit if keyName is Nothing or blank
                        if (string.IsNullOrEmpty(keyName))
                        {
                            return true;
                        }
                    }
                    else
                    {
                        // can't create section
                        return false;
                    }
                }
                if (keyName == null)
                {
                    // Delete the section
                    return DeleteSection(sectionName);
                }

                item = GetItem(sectionName, keyName);
                if ((item != null))
                {
                    if (newValue == null)
                    {
                        // delete this item
                        return DeleteItem(sectionName, keyName);
                    }
                    else
                    {
                        // add or update the value attribute
                        item.SetAttribute("value", newValue);
                        return true;
                    }
                }
                else
                {
                    // try to create the item
                    if ((!string.IsNullOrEmpty(keyName)) && ((newValue != null)))
                    {
                        // construct a new item (blank values are OK)
                        item = m_XmlDoc.CreateElement("item");
                        item.SetAttribute("key", SetNameCase(keyName));
                        item.SetAttribute("value", newValue);
                        section.AppendChild(item);
                        return true;
                    }
                }
                return false;
            }

            /// <summary>
            /// The function deletes a section in the file.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <return>The function returns a boolean that shows if the delete was completed.</return>
            private bool DeleteSection(string sectionName)
            {
                XmlElement section = GetSection(sectionName);
                if ((section != null))
                {
                    section.ParentNode.RemoveChild(section);
                    UpdateSections();
                    return true;
                }
                return false;
            }

            /// <summary>
            /// The function deletes a item in a specific section.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="keyName">The name of the key.</param>
            /// <return>The function returns a boolean that shows if the delete was completed.</return>
            private bool DeleteItem(string sectionName, string keyName)
            {
                XmlElement item = GetItem(sectionName, keyName);
                if ((item != null))
                {
                    item.ParentNode.RemoveChild(item);
                    return true;
                }
                return false;
            }

            /// <summary>
            /// Legacy function name; calls SetXmlKey
            /// </summary>
            public bool SetIniKey(string sectionName, string keyName, string newValue)
            {
                return SetXmlKey(sectionName, keyName, newValue);
            }

            /// <summary>
            /// Sets a new value for the "key" attribute.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="keyName">The name of the key.</param>
            /// <param name="newValue">The new value for the "key".</param>
            /// <return>The function returns a boolean that shows if the change was done.</return>
            public bool SetXmlKey(string sectionName, string keyName, string newValue)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                XmlElement item = GetItem(sectionName, keyName);
                if ((item != null))
                {
                    item.SetAttribute("key", SetNameCase(newValue));
                    return true;
                }
                return false;
            }

            /// <summary>
            /// Legacy function name; calls GetXMLValue
            /// </summary>
            public string GetIniValue(string sectionName, string keyName)
            {
                return GetXMLValue(sectionName, keyName);
            }

            /// <summary>
            /// Get the value for the given key in the given section
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="keyName">The name of the key.</param>
            /// <return>The function returns the name of the "value" attribute.</return>
            public string GetXMLValue(string sectionName, string keyName)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                XmlNode node = GetItem(sectionName, keyName);
                if ((node != null))
                {
                    return (node.Attributes.GetNamedItem("value").Value);
                }
                return null;
            }

            /// <summary>
            /// Legacy function name; calls GetXmlSectionComments
            /// </summary>
            public List<string> GetIniComments(string sectionName)
            {
                return GetXmlSectionComments(sectionName);
            }

            /// <summary>
            /// Gets the comments for a section name.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            ///<return>The function returns a string collection with comments</return>
            public List<string> GetXmlSectionComments(string sectionName)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                var sc = new List<string>();
                XmlNode target = default(XmlNode);

                if (sectionName == null)
                {
                    target = m_XmlDoc.DocumentElement;
                }
                else
                {
                    target = GetSection(sectionName);
                }
                if ((target != null))
                {
                    var nodes = target.SelectNodes("comment");
                    if (nodes.Count > 0)
                    {
                        foreach (XmlNode node in nodes)
                        {
                            sc.Add(node.InnerText);
                        }
                    }
                }
                return sc;
            }

            /// <summary>
            /// Legacy function name; calls SetXMLComments
            /// </summary>
            public bool SetIniComments(string sectionName, List<string> comments)
            {
                return SetXMLComments(sectionName, comments);
            }

            /// <summary>
            /// Sets a the comments for a section name.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="comments">A string collection.</param>
            ///<return>The function returns a Boolean that shows if the change was done.</return>
            public bool SetXMLComments(string sectionName, List<string> comments)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                XmlNode target = default(XmlNode);

                if (sectionName == null)
                {
                    target = m_XmlDoc.DocumentElement;
                }
                else
                {
                    target = GetSection(sectionName);
                }

                if ((target != null))
                {
                    var nodes = target.SelectNodes("comment");
                    foreach (XmlNode node in nodes)
                    {
                        target.RemoveChild(node);
                    }

                    foreach (var s in comments)
                    {
                        var node = m_XmlDoc.CreateElement("comment");
                        node.InnerText = s;
                        var nodeLastComment = (XmlElement) target.SelectSingleNode("comment[last()]");
                        if (nodeLastComment == null)
                        {
                            target.PrependChild(node);
                        }
                        else
                        {
                            target.InsertAfter(node, nodeLastComment);
                        }
                    }
                    return true;
                }
                return false;
            }

            /// <summary>
            /// The subroutine updades the sections.
            /// </summary>
            private void UpdateSections()
            {
                m_sections = new List<string>();
                foreach (XmlElement node in m_XmlDoc.SelectNodes("sections/section"))
                {
                    m_sections.Add(node.GetAttribute("name"));
                }
            }

            /// <summary>
            /// The subroutine gets the sections.
            /// </summary>
            /// <return>The subroutine returns a strin collection of sections.</return>
            public List<string> AllSections
            {
                get
                {
                    if (!Initialized)
                    {
                        return new List<string>();
                    }
                    else
                    {
                        return m_sections;
                    }
                }
            }

            /// <summary>
            /// Gets a collection of items for a section name.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="itemType">Item type.</param>
            /// <return>The function returns a string colection of items in a section.</return>
            private List<string> GetItemsInSection(string sectionName, XMLItemTypeEnum itemType)
            {
                var items = new List<string>();
                XmlNode section = GetSection(sectionName);

                if (section == null)
                {
                    return null;
                }
                else
                {
                    var nodes = section.SelectNodes("item");
                    if (nodes.Count > 0)
                    {
                        foreach (XmlNode node in nodes)
                        {
                            switch (itemType)
                            {
                                case XMLItemTypeEnum.GetKeys:
                                    items.Add(node.Attributes.GetNamedItem("key").Value);
                                    break;
                                case XMLItemTypeEnum.GetValues:
                                    items.Add(node.Attributes.GetNamedItem("value").Value);
                                    break;
                                case XMLItemTypeEnum.GetKeysAndValues:
                                    items.Add(node.Attributes.GetNamedItem("key").Value + "=" +
                                              node.Attributes.GetNamedItem("value").Value);
                                    break;
                            }
                        }
                    }
                    return items;
                }
            }

            /// <summary>The funtions gets a collection of keys in a section.</summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <return>The function returns a string colection of all the keys in a section.</return>
            public List<string> AllKeysInSection(string sectionName)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                return GetItemsInSection(sectionName, XMLItemTypeEnum.GetKeys);
            }

            /// <summary>The funtions gets a collection of values in a section.</summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <return>The function returns a string colection of all the values in a section.</return>
            public List<string> AllValuesInSection(string sectionName)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                return GetItemsInSection(sectionName, XMLItemTypeEnum.GetValues);
            }

            /// <summary>The funtions gets a collection of items in a section.</summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <return>The function returns a string colection of all the items in a section.</return>
            public List<string> AllItemsInSection(string sectionName)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                return (GetItemsInSection(sectionName, XMLItemTypeEnum.GetKeysAndValues));
            }

            /// <summary>The funtions gets a custom attribute name.</summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="keyName">The name of the key.</param>
            /// <param name="attributeName">The name of the attribute.</param>
            /// <return>The function returns a string.</return>
            public string GetCustomIniAttribute(string sectionName, string keyName, string attributeName)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                if (!string.IsNullOrEmpty(attributeName))
                {
                    var node = GetItem(sectionName, keyName);
                    if ((node != null))
                    {
                        attributeName = SetNameCase(attributeName);
                        return node.GetAttribute(attributeName);
                    }
                }
                return null;
            }

            /// <summary>The funtions sets a custom attribute name.</summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="keyName">The name of the key.</param>
            /// <param name="attributeName">The name of the attribute.</param>
            /// <param name="attributeValue">The value of the attribute.</param>
            /// <return>The function returns a Boolean.</return>
            public bool SetCustomIniAttribute(string sectionName, string keyName, string attributeName, string attributeValue)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                if (!string.IsNullOrEmpty(attributeName))
                {
                    var node = GetItem(sectionName, keyName);
                    if ((node != null))
                    {
                        try
                        {
                            if (attributeValue == null)
                            {
                                // delete the attribute
                                node.RemoveAttribute(attributeName);
                                return true;
                            }
                            else
                            {
                                attributeName = SetNameCase(attributeName);
                                node.SetAttribute(attributeName, attributeValue);
                                return true;
                            }
                        }
                        catch (Exception)
                        {
                            if (NotifyOnException)
                            {
                                throw new Exception("Failed to create item.");
                            }
                        }
                    }
                    return false;
                }
                return false;
            }

            /// <summary>The funtions creates a section name.</summary>
            /// <param name="sectionName">The name of the section to be created.</param>
            /// <return>The function returns a Boolean.</return>
            private bool CreateSection(string sectionName)
            {
                if (!string.IsNullOrEmpty(sectionName))
                {
                    sectionName = SetNameCase(sectionName);
                    try
                    {
                        var node = m_XmlDoc.CreateElement("section");
                        var attribute = m_XmlDoc.CreateAttribute("name");
                        attribute.Value = SetNameCase(sectionName);
                        node.Attributes.SetNamedItem(attribute);
                        m_XmlDoc.DocumentElement.AppendChild(node);
                        m_sections.Add(attribute.Value);
                        return true;
                    }
                    catch (Exception)
                    {
                        if (NotifyOnException)
                        {
                            throw new Exception("Failed to create item.");
                        }
                        return false;
                    }
                }
                return false;
            }

            /// <summary>
            /// Manually read a XML or .INI settings file line-by-line, extracting out any settings in the expected format
            /// </summary>
            /// <param name="strFilePath"></param>
            /// <returns></returns>
            /// <remarks></remarks>
            public bool ManualParseXmlOrIniFile(string strFilePath)
            {
                // Create a new, blank XML document
                m_XmlDoc.LoadXml(@"<?xml version=""1.0"" encoding=""UTF-8""?><sections></sections>");

                try
                {
                    var fi = new FileInfo(strFilePath);

                    if ((fi.Exists))
                    {
                        // Read strFilePath line-by-line to see if it has any .Ini style settings
                        // For example:
                        //   [SectionName]
                        //   Setting1=ValueA
                        //   Setting2=ValueB

                        // Also look for XML-style entries
                        // For example:
                        //   <section name="SectionName">
                        //     <item key="Setting1" value="ValueA" />
                        //   </section>

                        using (var srInFile = new StreamReader(new FileStream(fi.FullName, FileMode.Open, FileAccess.Read, FileShare.ReadWrite)))
                        {
                            while (srInFile.Peek() > -1)
                            {
                                var s = srInFile.ReadLine();

                                // Try to manually parse this line
                                ParseLineManual(s, m_XmlDoc);
                            }

                            m_XmlFilename = strFilePath;
                            m_initialized = true;
                        }
                    }
                    else
                    {
                        // File doesn't exist; create a new, blank .XML file
                        m_XmlFilename = strFilePath;
                        m_XmlDoc.Save(m_XmlFilename);
                        m_initialized = true;
                    }

                    return true;
                }
                catch (Exception)
                {
                    if (NotifyOnException)
                    {
                        throw new Exception("Failed to read XML file.");
                    }
                }

                return false;
            }

            /// <summary>Manually parses a line to extract the settings information
            /// Supports the traditional .Ini file format
            /// Also supports the 'key="KeyName" value="Value"' method used in XML settings files
            /// If success, then adds attributes to the doc object</summary>
            /// <param name="strLine">The name of the string to be parse.</param>
            /// <param name="doc">The name of the System.Xml.XmlDocument.</param>
            /// <remarks>Returns True for blank lines</remarks>
            private void ParseLineManual(string strLine, XmlDocument doc)
            {
                const string SECTION_NAME_TAG = "<section name=";
                const string KEY_TAG = "key=";
                const string VALUE_TAG = "value=";

                strLine = strLine.TrimStart();
                if (string.IsNullOrWhiteSpace(strLine))
                {
                    return;
                }

                switch ((strLine.Substring(0, 1)))
                {
                    case "[":
                        // this is a section
                        // trim the first and last characters
                        strLine = strLine.TrimStart('[');
                        strLine = strLine.TrimEnd(']');
                        // create a new section element
                        CreateSection(strLine);
                        return;
                    case ";":
                        // new comment
                        var node = doc.CreateElement("comment");
                        node.InnerText = strLine.Substring(1);
                        GetLastSection().AppendChild(node);
                        return;
                    default:
                        // Look for typical XML settings file elements

                        string strKey = string.Empty;
                        if (ParseLineManualCheckTag(strLine, SECTION_NAME_TAG, out strKey))
                        {
                            // This is an XML-style section

                            // Create a new section element
                            CreateSection(strKey);
                            return;
                        }
                        else
                        {
                            string strValue = string.Empty;
                            if (ParseLineManualCheckTag(strLine, KEY_TAG, out strKey))
                            {
                                // This is an XML-style key

                                ParseLineManualCheckTag(strLine, VALUE_TAG, out strValue);
                            }
                            else
                            {
                                // split the string on the "=" sign, if present
                                if ((strLine.IndexOf('=') > 0))
                                {
                                    var parts = strLine.Split('=');
                                    strKey = parts[0].Trim();
                                    strValue = parts[1].Trim();
                                }
                                else
                                {
                                    strKey = strLine;
                                    strValue = string.Empty;
                                }
                            }

                            if (string.IsNullOrEmpty(strKey))
                            {
                                strKey = string.Empty;
                            }

                            if (string.IsNullOrEmpty(strValue))
                            {
                                strValue = string.Empty;
                            }

                            if (string.IsNullOrWhiteSpace(strKey))
                            {
                                return;
                            }
                            else
                            {
                                var blnAddSetting = true;

                                switch (strKey.ToLower().Trim())
                                {
                                    case "<sections>":
                                    case "</section>":
                                    case "</sections>":
                                        // Do not add a new key
                                        if (string.IsNullOrEmpty(strValue))
                                        {
                                            blnAddSetting = false;
                                        }

                                        break;
                                }

                                if (blnAddSetting)
                                {
                                    var node1 = doc.CreateElement("item");
                                    var nodeAttribute = doc.CreateAttribute("key");
                                    nodeAttribute.Value = SetNameCase(strKey);
                                    node1.Attributes.SetNamedItem(nodeAttribute);

                                    nodeAttribute = doc.CreateAttribute("value");
                                    nodeAttribute.Value = strValue;
                                    node1.Attributes.SetNamedItem(nodeAttribute);

                                    GetLastSection().AppendChild(node1);
                                }

                                return;
                            }
                        }

                        break;
                }
            }

            private bool ParseLineManualCheckTag(
                string strLine,
                string strTagTofind,
                out string strTagValue)
            {
                int intMatchIndex = 0;
                int intNextMatchIndex = 0;

                strTagValue = string.Empty;
                intMatchIndex = strLine.ToLower().IndexOf(strTagTofind, StringComparison.Ordinal);

                if (intMatchIndex >= 0)
                {
                    strTagValue = strLine.Substring(intMatchIndex + strTagTofind.Length);

                    if (strTagValue.StartsWith("\""))
                    {
                        strTagValue = strTagValue.Substring(1);
                    }

                    intNextMatchIndex = strTagValue.IndexOf("\"");
                    if (intNextMatchIndex >= 0)
                    {
                        strTagValue = strTagValue.Substring(0, intNextMatchIndex);
                    }

                    return true;
                }
                else
                {
                    return false;
                }
            }

            /// <summary>It Sets or Gets the output file name.</summary>
            public string OutputFilename
            {
                get
                {
                    if (!Initialized)
                    {
                        return string.Empty;
                    }
                    else
                    {
                        return m_SaveFilename;
                    }
                }
                set
                {
                    FileInfo fi = default(FileInfo);
                    if (!Initialized)
                        throw new XMLFileReaderNotInitializedException();
                    fi = new FileInfo(value);
                    if (!fi.Directory.Exists)
                    {
                        if (NotifyOnException)
                        {
                            throw new Exception("Invalid path for output file.");
                        }
                    }
                    else
                    {
                        m_SaveFilename = value;
                    }
                }
            }

            /// <summary>
            /// It saves the data to the Xml output file.
            /// </summary>
            public void Save()
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                if ((OutputFilename != null) && (m_XmlDoc != null))
                {
                    var fi = new FileInfo(OutputFilename);
                    if (!fi.Directory.Exists)
                    {
                        if (NotifyOnException)
                        {
                            throw new Exception("Invalid path.");
                        }
                        return;
                    }
                    if (fi.Exists)
                    {
                        fi.Delete();
                        m_XmlDoc.Save(OutputFilename);
                    }
                    else
                    {
                        m_XmlDoc.Save(OutputFilename);
                    }
                    if (NotifyOnEvent)
                    {
                        if (InformationMessage != null)
                        {
                            InformationMessage("File save complete.");
                        }
                    }
                }
                else
                {
                    if (NotifyOnException)
                    {
                        throw new Exception("Not Output File name specified.");
                    }
                }
            }

            /// <summary>
            /// Gets the System.Xml.XmlDocument.
            /// </summary>
            public XmlDocument XmlDoc
            {
                get
                {
                    if (!Initialized)
                    {
                        return new XmlDocument();
                    }
                    else
                    {
                        return m_XmlDoc;
                    }
                }
            }

            /// <summary>
            /// Converts an XML document to a string.
            /// </summary>
            /// <return>It returns the XML document formatted as a string.</return>
            public string XML
            {
                get
                {
                    if (!Initialized)
                    {
                        return string.Empty;
                    }

                    var sb = new StringBuilder();
                    using (var sw = new StringWriter(sb))
                    {
                        using (var xw = new XmlTextWriter(sw))
                        {
                            xw.Indentation = 3;
                            xw.Formatting = Formatting.Indented;
                            m_XmlDoc.WriteContentTo(xw);
                        }
                    }

                    return sb.ToString();
                }
            }
        }

        public class XMLFileReaderNotInitializedException : ApplicationException
        {
            public override string Message
            {
                get { return "The XMLFileReader instance has not been properly initialized."; }
            }
        }
    }
}
