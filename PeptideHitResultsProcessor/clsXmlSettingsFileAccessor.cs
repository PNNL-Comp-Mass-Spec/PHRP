using System;
using System.Collections;
using System.Collections.Specialized;
using System.IO;
using System.Text;
using System.Xml;

// This class is used to read or write settings in an Xml settings file
// Based on a class from the DMS Analysis Manager software written by Dave Clark and Gary Kiebel (PNNL, Richland, WA)
// Additional features added by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) in October 2003
// Copyright 2005, Battelle Memorial Institute
//
// Updated in October 2004 to truly be case-insensitive if IsCaseSensitive = False when calling LoadSettings()
// Updated in August 2007 to remove the PRISM.Logging functionality and to include class XMLFileReader inside class XmlSettingsFileAccessor
// Updated in December 2010 to rename objects from Ini to XML

namespace PeptideHitResultsProcessor
{
    public class XmlSettingsFileAccessor
    {
        public XmlSettingsFileAccessor()
        {
            mCaseSensitive = false;
            htSectionNames = new Hashtable();

            mCachedSection = new udtRecentSectionType();
            mCachedSection.SectionName = string.Empty;
            mCachedSection.htKeys = new Hashtable();
        }

        private struct udtRecentSectionType
        {
            public string SectionName;              // Stores the section name whose keys are cached; the section name is capitalized identically to that actually present in the Xml file
            public Hashtable htKeys;
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

        // When mCaseSensitive = False, then htSectionNames stores mapping between lowercase section name and actual section name stored in file
        //   If section is present more than once in file, then only grabs the last occurence of the section
        // When mCaseSensitive = True, then the mappings in htSectionNames are effectively not used
        private readonly Hashtable htSectionNames;
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
        /// <return>The function returns a boolean that shows if the file was successfully loaded.</return>
        public bool LoadSettings(string XmlSettingsFilePath)
        {
            return LoadSettings(XmlSettingsFilePath, false);
        }

        /// <summary>
        /// Loads the settings for the defined Xml Settings File
        /// </summary>
        /// <param name="XmlSettingsFilePath">The path to the XML settings file.</param>
        /// <param name="IsCaseSensitive">Case sensitive names if True.  Non-case sensitive if false.</param>
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

        /// <summary>Checks if a section is present in the settings file.</summary>
        /// <param name="sectionName">The name of the section to look for.</param>
        /// <return>The function returns a boolean that shows if the section is present.</return>
        public bool SectionPresent(string sectionName)
        {
            StringCollection strSections = default(StringCollection);

            strSections = m_XMLFileAccessor.AllSections;

            for (var intIndex = 0; intIndex <= strSections.Count - 1; intIndex++)
            {
                if (SetNameCase(strSections[intIndex]) == SetNameCase(sectionName))
                    return true;
            }

            return false;
        }

        private bool CacheKeyNames(string sectionName)
        {
            // Looks up the Key Names for the given section, storing them in mCachedSection
            // This is done so that this class will know the correct capitalization for the key names

            StringCollection strKeys = default(StringCollection);

            string sectionNameInFile = null;
            string strKeyNameToStore = null;

            // Lookup the correct capitalization for sectionName (only truly important if mCaseSensitive = False)
            sectionNameInFile = GetCachedSectionName(sectionName);
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
            mCachedSection.htKeys.Clear();

            for (var intIndex = 0; intIndex <= strKeys.Count - 1; intIndex++)
            {
                if (mCaseSensitive)
                {
                    strKeyNameToStore = string.Copy(strKeys[intIndex]);
                }
                else
                {
                    strKeyNameToStore = string.Copy(strKeys[intIndex].ToLower());
                }

                if (!mCachedSection.htKeys.Contains(strKeyNameToStore))
                {
                    mCachedSection.htKeys.Add(strKeyNameToStore, strKeys[intIndex]);
                }
            }

            return true;
        }

        private void CacheSectionNames()
        {
            // Looks up the Section Names in the XML file
            // This is done so that this class will know the correct capitalization for the section names

            StringCollection strSections = default(StringCollection);
            string strSectionNameToStore = null;

            strSections = m_XMLFileAccessor.AllSections;

            htSectionNames.Clear();

            for (var intIndex = 0; intIndex <= strSections.Count - 1; intIndex++)
            {
                if (mCaseSensitive)
                {
                    strSectionNameToStore = string.Copy(strSections[intIndex]);
                }
                else
                {
                    strSectionNameToStore = string.Copy(strSections[intIndex].ToLower());
                }

                if (!htSectionNames.Contains(strSectionNameToStore))
                {
                    htSectionNames.Add(strSectionNameToStore, strSections[intIndex]);
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
                if (mCachedSection.htKeys.ContainsKey(keyNameToFind))
                {
                    return Convert.ToString(mCachedSection.htKeys[keyNameToFind]);
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

            string sectionNameToFind = null;

            sectionNameToFind = SetNameCase(sectionName);
            if (htSectionNames.ContainsKey(sectionNameToFind))
            {
                return Convert.ToString(htSectionNames[sectionNameToFind]);
            }
            else
            {
                return string.Empty;
            }
        }

        private string SetNameCase(string aName)
        {
            // Changes aName to lowercase if mCaseSensitive = False

            if (mCaseSensitive)
            {
                return aName;
            }
            else
            {
                return aName.ToLower();
            }
        }

        /// <summary>
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a String.</return>
        public string GetParam(string sectionName, string keyName, string valueIfMissing)
        {
            bool blah = false;
            return GetParam(sectionName, keyName, valueIfMissing, ref blah);
        }

        /// <summary>
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.  Returned ByRef.</param>
        /// <return>The function returns the name of the "value" attribute as a String.</return>
        public string GetParam(string sectionName, string keyName, string valueIfMissing, ref bool valueNotPresent)
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
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns boolean True if the "value" attribute is "true".  Otherwise, returns boolean False.</return>
        public bool GetParam(string sectionName, string keyName, bool valueIfMissing)
        {
            bool blah = false;
            return GetParam(sectionName, keyName, valueIfMissing, ref blah);
        }

        /// <summary>
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.  Returned ByRef.</param>
        /// <return>The function returns boolean True if the "value" attribute is "true".  Otherwise, returns boolean False.</return>
        public bool GetParam(string sectionName, string keyName, bool valueIfMissing, ref bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), ref blnNotFound);
            if (strResult == null || blnNotFound)
            {
                valueNotPresent = true;
                return valueIfMissing;
            }
            else
            {
                valueNotPresent = false;
                if (strResult.ToLower() == "true")
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
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a Short.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public short GetParam(string sectionName, string keyName, short valueIfMissing)
        {
            bool blah = false;
            return GetParam(sectionName, keyName, valueIfMissing, ref blah);
        }

        /// <summary>
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.  Returned ByRef.</param>
        /// <return>The function returns the name of the "value" attribute as a Short.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public short GetParam(string sectionName, string keyName, short valueIfMissing, ref bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;
            short intValue = 0;

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), ref blnNotFound);
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
                    else if (strResult.ToLower() == "true")
                    {
                        return -1;
                    }
                    else if (strResult.ToLower() == "false")
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
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as an Integer.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public int GetParam(string sectionName, string keyName, int valueIfMissing)
        {
            bool blah = false;
            return GetParam(sectionName, keyName, valueIfMissing, ref blah);
        }

        /// <summary>
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.  Returned ByRef.</param>
        /// <return>The function returns the name of the "value" attribute as an Integer.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public int GetParam(string sectionName, string keyName, int valueIfMissing, ref bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;
            int intValue = 0;

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), ref blnNotFound);
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
                    else if (strResult.ToLower() == "true")
                    {
                        return -1;
                    }
                    else if (strResult.ToLower() == "false")
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
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a Long.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public Int64 GetParam(string sectionName, string keyName, long valueIfMissing)
        {
            bool blah = false;
            return GetParam(sectionName, keyName, valueIfMissing, ref blah);
        }

        /// <summary>
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.  Returned ByRef.</param>
        /// <return>The function returns the name of the "value" attribute as a Long.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public Int64 GetParam(string sectionName, string keyName, long valueIfMissing, ref bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;
            Int64 intValue = default(Int64);

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), ref blnNotFound);
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
                    else if (strResult.ToLower() == "true")
                    {
                        return -1;
                    }
                    else if (strResult.ToLower() == "false")
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
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a Single.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public float GetParam(string sectionName, string keyName, float valueIfMissing)
        {
            bool blah = false;
            return GetParam(sectionName, keyName, valueIfMissing, ref blah);
        }

        /// <summary>
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.  Returned ByRef.</param>
        /// <return>The function returns the name of the "value" attribute as a Single.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public float GetParam(string sectionName, string keyName, float valueIfMissing, ref bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;
            float sngValue = 0;

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), ref blnNotFound);
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
                    else if (strResult.ToLower() == "true")
                    {
                        return -1;
                    }
                    else if (strResult.ToLower() == "false")
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
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <return>The function returns the name of the "value" attribute as a Double.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public double GetParam(string sectionName, string keyName, double valueIfMissing)
        {
            bool blah = false;
            return GetParam(sectionName, keyName, valueIfMissing, ref blah);
        }

        /// <summary>
        /// The function gets the name of the "value" attribute in section "sectionName".
        /// </summary>
        /// <param name="sectionName">The name of the section.</param>
        /// <param name="keyName">The name of the key.</param>
        /// <param name="valueIfMissing">Value to return if "sectionName" or "keyName" is missing.</param>
        /// <param name="valueNotPresent">Set to True if "sectionName" or "keyName" is missing.  Returned ByRef.</param>
        /// <return>The function returns the name of the "value" attribute as a Double.  If "value" is "true" returns -1.  If "value" is "false" returns 0.</return>
        public double GetParam(string sectionName, string keyName, double valueIfMissing, ref bool valueNotPresent)
        {
            string strResult = null;
            var blnNotFound = false;
            double dblValue = 0;

            strResult = this.GetParam(sectionName, keyName, valueIfMissing.ToString(), ref blnNotFound);
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
                    else if (strResult.ToLower() == "true")
                    {
                        return -1;
                    }
                    else if (strResult.ToLower() == "false")
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
        public void SetIniFilePath(string XmlSettingsFilePath)
        {
            SetXMLFilePath(XmlSettingsFilePath);
        }

        /// <summary>
        /// The function sets the path to the Xml Settings File.
        /// </summary>
        /// <param name="XmlSettingsFilePath">The path to the XML settings file.</param>
        public void SetXMLFilePath(string XmlSettingsFilePath)
        {
            m_XMLFilePath = XmlSettingsFilePath;
        }

        /// <summary>
        /// The function sets a new String value for the "value" attribute.
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
        /// The function sets a new Boolean value for the "value" attribute.
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
        /// The function sets a new Short value for the "value" attribute.
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
        /// The function sets a new Integer value for the "value" attribute.
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
        /// The function sets a new Long value for the "value" attribute.
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
        /// The function sets a new Single value for the "value" attribute.
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
        /// The function sets a new Double value for the "value" attribute.
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
        /// The function renames a section.
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

            private StringCollection sections = new StringCollection();
            private bool m_CaseSensitive = false;
            private string m_SaveFilename;
            private bool m_initialized = false;

            public bool NotifyOnEvent;
            public bool NotifyOnException;

            public event InformationMessageEventHandler InformationMessage;
            public delegate void InformationMessageEventHandler(string msg);

            /// <summary>Initializes a new instance of the XMLFileReader (non case-sensitive)</summary>
            /// <param name="XmlFilename">The name of the XML file.</param>
            public XMLFileReader(string XmlFilename)
            {
                NotifyOnException = false;
                InitXMLFileReader(XmlFilename, false);
            }

            /// <summary>Initializes a new instance of the XMLFileReader.</summary>
            /// <param name="XmlFilename">The name of the XML file.</param>
            /// <param name="IsCaseSensitive">Case sensitive as boolean.</param>
            public XMLFileReader(string XmlFilename, bool IsCaseSensitive)
            {
                NotifyOnException = true;
                InitXMLFileReader(XmlFilename, IsCaseSensitive);
            }

            /// <summary>
            /// This routine is called by each of the constructors to make the actual assignments.
            /// </summary>
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
            /// This routine returns a boolean showing if the file was initialized or not.
            /// </summary>
            /// <return>The function returns a Boolean.</return>
            public bool Initialized
            {
                get { return m_initialized; }
            }

            /// <summary>
            /// This routine returns a boolean showing if the name is case sensitive or not.
            /// </summary>
            /// <return>The function returns a Boolean.</return>
            public bool CaseSensitive
            {
                get { return m_CaseSensitive; }
            }

            /// <summary>
            /// This routine sets a name.
            /// </summary>
            /// <param name="aName">The name to be set.</param>
            /// <return>The function returns a string.</return>
            private string SetNameCase(string aName)
            {
                if ((CaseSensitive))
                {
                    return aName;
                }
                else
                {
                    return aName.ToLower();
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
            /// The function gets the last section.
            /// </summary>
            /// <return>The function returns the last section as System.Xml.XmlElement.</return>
            private XmlElement GetLastSection()
            {
                if (sections.Count == 0)
                {
                    return GetRoot();
                }
                else
                {
                    return GetSection(sections[sections.Count - 1]);
                }
            }

            /// <summary>
            /// The function gets a section as System.Xml.XmlElement.
            /// </summary>
            /// <param name="sectionName">The name of a section.</param>
            /// <return>The function returns a section as System.Xml.XmlElement.</return>
            private XmlElement GetSection(string sectionName)
            {
                if ((!(sectionName == null)) && (!string.IsNullOrEmpty(sectionName)))
                {
                    sectionName = SetNameCase(sectionName);
                    return (XmlElement)m_XmlDoc.SelectSingleNode("//section[@name='" + sectionName + "']");
                }
                return null;
            }

            /// <summary>
            /// The function gets an item.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="keyName">The name of the key.</param>
            /// <return>The function returns a XML element.</return>
            private XmlElement GetItem(string sectionName, string keyName)
            {
                XmlElement section = default(XmlElement);
                if (((keyName != null)) && (!string.IsNullOrEmpty(keyName)))
                {
                    keyName = SetNameCase(keyName);
                    section = GetSection(sectionName);
                    if (((section != null)))
                    {
                        return (XmlElement)section.SelectSingleNode("item[@key='" + keyName + "']");
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
            /// The function sets the ini section name.
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
                if (((newSection != null)) && (!string.IsNullOrEmpty(newSection)))
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
            /// The function sets a new value for the "value" attribute.
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
                        if ((keyName == null) || (string.IsNullOrEmpty(keyName)))
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
                    // delete the section
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
            /// The function sets a new value for the "key" attribute.
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
            /// The function gets the name of the "value" attribute.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="keyName">The name of the key.</param>
            ///<return>The function returns the name of the "value" attribute.</return>
            public string GetXMLValue(string sectionName, string keyName)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                XmlNode N = GetItem(sectionName, keyName);
                if ((N != null))
                {
                    return (N.Attributes.GetNamedItem("value").Value);
                }
                return null;
            }

            /// <summary>
            /// Legacy function name; calls GetXmlSectionComments
            /// </summary>
            public StringCollection GetIniComments(string sectionName)
            {
                return GetXmlSectionComments(sectionName);
            }

            /// <summary>
            /// The function gets the comments for a section name.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            ///<return>The function returns a string collection with comments</return>
            public StringCollection GetXmlSectionComments(string sectionName)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                var sc = new StringCollection();
                XmlNode target = default(XmlNode);
                XmlNodeList nodes = default(XmlNodeList);
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
                    nodes = target.SelectNodes("comment");
                    if (nodes.Count > 0)
                    {
                        foreach (XmlNode N in nodes)
                        {
                            sc.Add(N.InnerText);
                        }
                    }
                }
                return sc;
            }

            /// <summary>
            /// Legacy function name; calls SetXMLComments
            /// </summary>
            public bool SetIniComments(string sectionName, StringCollection comments)
            {
                return SetXMLComments(sectionName, comments);
            }

            /// <summary>
            /// The function sets a the comments for a section name.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="comments">A string collection.</param>
            ///<return>The function returns a Boolean that shows if the change was done.</return>
            public bool SetXMLComments(string sectionName, StringCollection comments)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                XmlNode target = default(XmlNode);
                XmlNodeList nodes = default(XmlNodeList);
                XmlElement NLastComment = default(XmlElement);
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
                    nodes = target.SelectNodes("comment");
                    foreach (XmlNode N in nodes)
                    {
                        target.RemoveChild(N);
                    }
                    foreach (string s in comments)
                    {
                        XmlNode N = m_XmlDoc.CreateElement("comment");
                        N.InnerText = s;
                        NLastComment = (XmlElement)target.SelectSingleNode("comment[last()]");
                        if (NLastComment == null)
                        {
                            target.PrependChild(N);
                        }
                        else
                        {
                            target.InsertAfter(N, NLastComment);
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
                sections = new StringCollection();
                foreach (XmlElement N in m_XmlDoc.SelectNodes("sections/section"))
                {
                    sections.Add(N.GetAttribute("name"));
                }
            }
            /// <summary>
            /// The subroutine gets the sections.
            /// </summary>
            /// <return>The subroutine returns a strin collection of sections.</return>
            public StringCollection AllSections
            {
                get
                {
                    if (!Initialized)
                    {
                        return new StringCollection();
                    }
                    else
                    {
                        return sections;
                    }
                }
            }

            /// <summary>
            /// The function gets a collection of items for a section name.
            /// </summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="itemType">Item type.</param>
            /// <return>The function returns a string colection of items in a section.</return>
            private StringCollection GetItemsInSection(string sectionName, XMLItemTypeEnum itemType)
            {
                XmlNodeList nodes = default(XmlNodeList);
                var items = new StringCollection();
                XmlNode section = GetSection(sectionName);
                if (section == null)
                {
                    return null;
                }
                else
                {
                    nodes = section.SelectNodes("item");
                    if (nodes.Count > 0)
                    {
                        foreach (XmlNode N in nodes)
                        {
                            switch (itemType)
                            {
                                case XMLItemTypeEnum.GetKeys:
                                    items.Add(N.Attributes.GetNamedItem("key").Value);
                                    break;
                                case XMLItemTypeEnum.GetValues:
                                    items.Add(N.Attributes.GetNamedItem("value").Value);
                                    break;
                                case XMLItemTypeEnum.GetKeysAndValues:
                                    items.Add(N.Attributes.GetNamedItem("key").Value + "=" + N.Attributes.GetNamedItem("value").Value);
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
            public StringCollection AllKeysInSection(string sectionName)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                return GetItemsInSection(sectionName, XMLItemTypeEnum.GetKeys);
            }

            /// <summary>The funtions gets a collection of values in a section.</summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <return>The function returns a string colection of all the values in a section.</return>
            public StringCollection AllValuesInSection(string sectionName)
            {
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                return GetItemsInSection(sectionName, XMLItemTypeEnum.GetValues);
            }

            /// <summary>The funtions gets a collection of items in a section.</summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <return>The function returns a string colection of all the items in a section.</return>
            public StringCollection AllItemsInSection(string sectionName)
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
                XmlElement N = default(XmlElement);
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                if (((attributeName != null)) && (!string.IsNullOrEmpty(attributeName)))
                {
                    N = GetItem(sectionName, keyName);
                    if ((N != null))
                    {
                        attributeName = SetNameCase(attributeName);
                        return N.GetAttribute(attributeName);
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
                XmlElement N = default(XmlElement);
                if (!Initialized)
                    throw new XMLFileReaderNotInitializedException();
                if (!string.IsNullOrEmpty(attributeName))
                {
                    N = GetItem(sectionName, keyName);
                    if ((N != null))
                    {
                        try
                        {
                            if (attributeValue == null)
                            {
                                // delete the attribute
                                N.RemoveAttribute(attributeName);
                                return true;
                            }
                            else
                            {
                                attributeName = SetNameCase(attributeName);
                                N.SetAttribute(attributeName, attributeValue);
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
                XmlElement N = default(XmlElement);
                XmlAttribute Natt = default(XmlAttribute);
                if (((sectionName != null)) && (!string.IsNullOrEmpty(sectionName)))
                {
                    sectionName = SetNameCase(sectionName);
                    try
                    {
                        N = m_XmlDoc.CreateElement("section");
                        Natt = m_XmlDoc.CreateAttribute("name");
                        Natt.Value = SetNameCase(sectionName);
                        N.Attributes.SetNamedItem(Natt);
                        m_XmlDoc.DocumentElement.AppendChild(N);
                        sections.Add(Natt.Value);
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

            /// <summary>The funtions creates a section name.</summary>
            /// <param name="sectionName">The name of the section.</param>
            /// <param name="keyName">The name of the key.</param>
            /// <param name="newValue">The new value to be created.</param>
            /// <return>The function returns a Boolean.</return>
            private bool CreateItem(string sectionName, string keyName, string newValue)
            {
                XmlElement item = default(XmlElement);
                XmlElement section = default(XmlElement);
                try
                {
                    section = GetSection(sectionName);
                    if ((section != null))
                    {
                        item = m_XmlDoc.CreateElement("item");
                        item.SetAttribute("key", keyName);
                        item.SetAttribute("newValue", newValue);
                        section.AppendChild(item);
                        return true;
                    }
                    return false;
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

            /// <summary>
            /// Manually read a XML or .INI settings file line-by-line, extracting out any settings in the expected format
            /// </summary>
            /// <param name="strFilePath"></param>
            /// <returns></returns>
            /// <remarks></remarks>
            public bool ManualParseXmlOrIniFile(string strFilePath)
            {
                // Create a new, blank XML document
                m_XmlDoc.LoadXml("<?xml version=\"1.0\" encoding=\"UTF-8\"?><sections></sections>");

                try
                {
                    FileInfo fi = default(FileInfo);
                    string s = null;

                    fi = new FileInfo(strFilePath);
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
                                s = srInFile.ReadLine();

                                // Try to manually parse this line
                                ParseLineManual(s, ref m_XmlDoc);
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
            /// <returns>True if success, false if not a recognized line format</returns>
            /// <remarks>Returns True for blank lines</remarks>
            private bool ParseLineManual(string strLine, ref XmlDocument doc)
            {
                const string SECTION_NAME_TAG = "<section name=";
                const string KEY_TAG = "key=";
                const string VALUE_TAG = "value=";

                string strKey = string.Empty;
                string strValue = string.Empty;
                bool blnAddSetting = false;

                XmlElement N = default(XmlElement);
                XmlAttribute Natt = default(XmlAttribute);
                string[] parts = null;

                strLine = strLine.TrimStart();
                if (string.IsNullOrWhiteSpace(strLine))
                {
                    return true;
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
                        return true;
                    case ";":
                        // new comment
                        N = doc.CreateElement("comment");
                        N.InnerText = strLine.Substring(1);
                        GetLastSection().AppendChild(N);
                        return true;
                    default:
                        // Look for typical XML settings file elements

                        if (ParseLineManualCheckTag(strLine, SECTION_NAME_TAG, ref strKey))
                        {
                            // This is an XML-style section

                            // Create a new section element
                            CreateSection(strKey);
                            return true;
                        }
                        else
                        {
                            if (ParseLineManualCheckTag(strLine, KEY_TAG, ref strKey))
                            {
                                // This is an XML-style key

                                ParseLineManualCheckTag(strLine, VALUE_TAG, ref strValue);
                            }
                            else
                            {
                                // split the string on the "=" sign, if present
                                if ((strLine.IndexOf('=') > 0))
                                {
                                    parts = strLine.Split('=');
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
                                return false;
                            }
                            else
                            {
                                blnAddSetting = true;

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
                                    N = doc.CreateElement("item");
                                    Natt = doc.CreateAttribute("key");
                                    Natt.Value = SetNameCase(strKey);
                                    N.Attributes.SetNamedItem(Natt);

                                    Natt = doc.CreateAttribute("value");
                                    Natt.Value = strValue;
                                    N.Attributes.SetNamedItem(Natt);

                                    GetLastSection().AppendChild(N);
                                }

                                return true;
                            }
                        }
                }
            }

            private bool ParseLineManualCheckTag(
                string strLine,
                string strTagTofind,
                ref string strTagValue)
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

                    intNextMatchIndex = strTagValue.IndexOf('"');
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
            /// <summary>It saves the data to the Xml output file.</summary>
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

            /// <summary>It gets the System.Xml.XmlDocument.</summary>
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

            /// <summary>Converts an XML document to a string.</summary>
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