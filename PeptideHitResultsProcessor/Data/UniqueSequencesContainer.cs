// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 11, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://github.com/PNNL-Comp-Mass-Spec/ or https://panomics.pnnl.gov/ or https://www.pnnl.gov/integrative-omics
// -------------------------------------------------------------------------------
//
// Licensed under the 2-Clause BSD License; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// https://opensource.org/licenses/BSD-2-Clause
//
// Copyright 2018 Battelle Memorial Institute

using System;
using System.Collections.Generic;

namespace PeptideHitResultsProcessor.Data
{
    /// <summary>
    /// This class is used to track peptide sequences and their modification descriptors
    /// It assigns a unique integer ID to each combination of sequence and modification description
    /// </summary>
    public class UniqueSequencesContainer
    {
        private const int DEFAULT_INITIAL_SEQ_ID = 1;
        private const char SEQUENCE_MOD_DESC_SEP = '_';

        private readonly Dictionary<string, int> mMasterSequences;

        private int mNextUniqueSeqID;

        /// <summary>
        /// Unique sequence count
        /// </summary>
        // ReSharper disable once UnusedMember.Global
        public int UniqueSequenceCount => mMasterSequences.Count;

        /// <summary>
        /// Constructor
        /// </summary>
        public UniqueSequencesContainer()
        {
            mMasterSequences = new Dictionary<string, int>();
            Clear();
        }

        /// <summary>
        /// Clear stored sequences, resetting the initial sequence ID to 1
        /// </summary>
        public void Clear()
        {
            Clear(DEFAULT_INITIAL_SEQ_ID);
        }

        /// <summary>
        /// Clear stored sequences, resetting the initial sequence ID to the given value
        /// </summary>
        /// <param name="initialSeqID"></param>
        public void Clear(int initialSeqID)
        {
            // Clears mMasterSequences and resets mNextUniqueSeqID to initialSeqID
            mMasterSequences.Clear();
            mNextUniqueSeqID = initialSeqID;
        }

        /// <summary>
        /// Look for the given sequence in the master sequence list, returning the sequence ID if found, or adding it if missing
        /// </summary>
        /// <param name="sequence"></param>
        /// <param name="modDescription"></param>
        /// <param name="existingSequenceFound"></param>
        public int GetNextUniqueSequenceID(string sequence, string modDescription, out bool existingSequenceFound)
        {
            int uniqueSeqID;

            existingSequenceFound = false;

            try
            {
                sequence ??= string.Empty;

                modDescription ??= string.Empty;

                var key = sequence + SEQUENCE_MOD_DESC_SEP + modDescription;

                if (mMasterSequences.ContainsKey(key))
                {
                    uniqueSeqID = mMasterSequences[key];
                    existingSequenceFound = true;
                }
                else
                {
                    mMasterSequences.Add(key, mNextUniqueSeqID);
                    uniqueSeqID = mNextUniqueSeqID;
                    mNextUniqueSeqID++;
                }
            }
            catch (Exception)
            {
                uniqueSeqID = 0;
            }

            return uniqueSeqID;
        }
    }
}
