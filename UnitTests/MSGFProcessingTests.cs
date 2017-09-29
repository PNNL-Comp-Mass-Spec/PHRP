using System;
using System.Collections.Generic;
using NUnit.Framework;
using PeptideHitResultsProcessor;
using PHRPReader;

namespace PHRP_UnitTests
{
    [TestFixture]
    public class MSGFProcessingTests
    {
        [Test]
        public void TestReplaceMSGFModTextWithSymbol()
        {
            var modInfo = new List<clsMSGFPlusParamFileModExtractor.udtModInfoType>();
            modInfo.Add(new clsMSGFPlusParamFileModExtractor.udtModInfoType()
            {
                ModMass = "C2H3N1O1",
                ModMassVal = 57.0214619,
                ModName = "Carbamidomethyl",
                ModSymbol = '-',
                ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.StaticMod,
                Residues = "C",
            });
            modInfo.Add(new clsMSGFPlusParamFileModExtractor.udtModInfoType()
            {
                ModMass = "229.1629",
                ModMassVal = 229.1629,
                ModName = "TMT6plex",
                ModSymbol = '-',
                ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.StaticMod,
                Residues = "<",
            });
            modInfo.Add(new clsMSGFPlusParamFileModExtractor.udtModInfoType()
            {
                ModMass = "229.1629",
                ModMassVal = 229.1629,
                ModName = "TMT6plex",
                ModSymbol = '-',
                ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.StaticMod,
                Residues = "K",
            });
            modInfo.Add(new clsMSGFPlusParamFileModExtractor.udtModInfoType()
            {
                ModMass = "O1",
                ModMassVal = 15.9949141,
                ModName = "Oxidation",
                ModSymbol = '*',
                ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.DynamicMod,
                Residues = "M",
            });
            modInfo.Add(new clsMSGFPlusParamFileModExtractor.udtModInfoType()
            {
                ModMass = "-187.152366",
                ModMassVal = -187.152366,
                ModName = "AcNoTMT",
                ModSymbol = '#',
                ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.DynamicMod,
                Residues = "K",
            });
            modInfo.Add(new clsMSGFPlusParamFileModExtractor.udtModInfoType()
            {
                ModMass = "-187.152366",
                ModMassVal = -187.152366,
                ModName = "AcNoTMT",
                ModSymbol = '#',
                ModType = clsMSGFPlusParamFileModExtractor.eMSGFDBModType.DynNTermPeptide,
                Residues = "<",
            });

            var peptide = "-187.152+229.163ATK-187.152+229.163QIFDC+57.021K+229.163";
            double totalModMass;
            var processor = new clsMSGFDBResultsProcessor();

            var replaced = processor.ReplaceMSGFModTextWithSymbol(peptide, modInfo, true, out totalModMass);
            Console.WriteLine("Original: {0}", peptide);
            Console.WriteLine("Replaced: {0}", replaced);
            Assert.AreEqual("A#TK#QIFDCK", replaced);

            var peptide1 = "+229.163ITVVGVGAVGM+15.995AC+57.021AISILMK+229.163";
            var replaced1 = processor.ReplaceMSGFModTextWithSymbol(peptide1, modInfo, true, out totalModMass);
            Console.WriteLine("Original: {0}", peptide1);
            Console.WriteLine("Replaced: {0}", replaced1);
            Assert.AreEqual("ITVVGVGAVGM*ACAISILMK", replaced1);
        }
    }
}
