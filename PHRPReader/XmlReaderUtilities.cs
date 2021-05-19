using System.Xml.Linq;

namespace PHRPReader
{
    /// <summary>
    /// Utility methods for reading data from XML files using an XDocument
    /// </summary>
    public static class XmlReaderUtilities
    {
        /// <summary>
        /// Get the named attribute from the given element
        /// </summary>
        /// <param name="item"></param>
        /// <param name="attributeName"></param>
        /// <param name="attributeValue"></param>
        /// <returns>True if found, otherwise false</returns>
        public static bool TryGetAttribute(XElement item, string attributeName, out string attributeValue)
        {
            if (!item.HasAttributes)
            {
                attributeValue = string.Empty;
                return false;
            }

            var attribute = item.Attribute(attributeName);

            if (attribute == null)
            {
                attributeValue = string.Empty;
                return false;
            }

            attributeValue = attribute.Value;
            return true;
        }

        /// <summary>
        /// Look for the given element below the given parent
        /// If found, update elementValue with its value and return true
        /// Otherwise, return false
        /// </summary>
        /// <param name="parentItem"></param>
        /// <param name="elementName"></param>
        /// <param name="elementValue"></param>
        /// <returns>True if found, otherwise false</returns>
        /// <remarks>See also <see cref="GetElementValueOrDefault"/></remarks>
        public static bool TryGetElementValue(XElement parentItem, string elementName, out string elementValue)
        {
            var node = parentItem.Element(elementName);
            if (node == null)
            {
                elementValue = string.Empty;
                return false;
            }

            elementValue = node.Value;
            return true;
        }
    }
}
