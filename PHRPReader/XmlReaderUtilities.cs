using System.Xml.Linq;

namespace PHRPReader
{
    /// <summary>
    /// Utility methods for reading data from XML files using an XDocument
    /// </summary>
    public static class XmlReaderUtilities
    {
        /// <summary>
        /// Traverse an XML node hierarchy for the given element names
        /// If the final element is found, return its value, otherwise return an empty string
        /// </summary>
        /// <remarks>See also <see cref="TryGetElementValue"/></remarks>
        /// <param name="parentNode">Parent node</param>
        /// <param name="elementNames">Element names</param>
        /// <returns>Value if found, or an empty string</returns>
        public static string GetElementValueOrDefault(XElement parentNode, params string[] elementNames)
        {
            if (elementNames.Length == 0)
                return parentNode.Value;

            if (elementNames.Length == 1)
            {
                if (TryGetElementValue(parentNode, elementNames[0], out var value))
                    return value;

                return string.Empty;
            }

            for (var i = 0; i < elementNames.Length; i++)
            {
                var childNode = parentNode.Element(elementNames[i]);

                if (childNode == null)
                    return string.Empty;

                if (i == elementNames.Length - 1)
                {
                    return childNode.Value;
                }

                parentNode = childNode;
            }

            return string.Empty;
        }

        /// <summary>
        /// Get the named attribute from the given element
        /// </summary>
        /// <param name="item">XML element</param>
        /// <param name="attributeName">Attribute name</param>
        /// <param name="attributeValue">Attribute value</param>
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
        /// <remarks>See also <see cref="GetElementValueOrDefault"/></remarks>
        /// <param name="parentItem">Parent XML element</param>
        /// <param name="elementName">Element name</param>
        /// <param name="elementValue">Element value</param>
        /// <returns>True if found, otherwise false</returns>
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
