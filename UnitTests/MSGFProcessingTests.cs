﻿using System;
using System.Collections.Generic;
using NUnit.Framework;
using PeptideHitResultsProcessor;
using PeptideHitResultsProcessor.Processor;
using PHRPReader;

namespace PHRP_UnitTests
{
    [TestFixture]
    public class MSGFProcessingTests
    {
        // Ignore Spelling: Carbamidomethyl, PHRP

        [Test]
        public void TestReplaceMSGFModTextWithSymbol()
        {
            var modDef = new List<MSGFPlusParamFileModExtractor.ModInfo>
            {
                new()
                {
                    ModMass = "C2H3N1O1",
                    ModMassVal = 57.0214619,
                    ModName = "Carbamidomethyl",
                    ModSymbol = '-',
                    ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod,
                    Residues = "C",
                },
                new()
                {
                    ModMass = "229.1629",
                    ModMassVal = 229.1629,
                    ModName = "TMT6plex",
                    ModSymbol = '-',
                    ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod,
                    Residues = "<",
                },
                new()
                {
                    ModMass = "229.1629",
                    ModMassVal = 229.1629,
                    ModName = "TMT6plex",
                    ModSymbol = '-',
                    ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.StaticMod,
                    Residues = "K",
                },
                new()
                {
                    ModMass = "O1",
                    ModMassVal = 15.9949141,
                    ModName = "Oxidation",
                    ModSymbol = '*',
                    ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod,
                    Residues = "M",
                },
                new()
                {
                    ModMass = "-187.152366",
                    ModMassVal = -187.152366,
                    ModName = "AcNoTMT",
                    ModSymbol = '#',
                    ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.DynamicMod,
                    Residues = "K",
                },
                new()
                {
                    ModMass = "-187.152366",
                    ModMassVal = -187.152366,
                    ModName = "AcNoTMT",
                    ModSymbol = '#',
                    ModType = MSGFPlusParamFileModExtractor.MSGFPlusModType.DynNTermPeptide,
                    Residues = "<",
                }
            };

            // ReSharper disable StringLiteralTypo

            const string peptide = "-187.152+229.163ATK-187.152+229.163QIFDC+57.021K+229.163";
            var processor = new MSGFPlusResultsProcessor(new PHRPOptions());

            var replaced = processor.ReplaceMSGFModTextWithSymbol(peptide, modDef, true, out _);
            Console.WriteLine("Original: {0}", peptide);
            Console.WriteLine("Replaced: {0}", replaced);
            Assert.AreEqual("A#TK#QIFDCK", replaced);

            const string peptide1 = "+229.163ITVVGVGAVGM+15.995AC+57.021AISILMK+229.163";
            var replaced1 = processor.ReplaceMSGFModTextWithSymbol(peptide1, modDef, true, out _);
            Console.WriteLine("Original: {0}", peptide1);
            Console.WriteLine("Replaced: {0}", replaced1);
            Assert.AreEqual("ITVVGVGAVGM*ACAISILMK", replaced1);

            // ReSharper restore StringLiteralTypo
        }
    }
}
