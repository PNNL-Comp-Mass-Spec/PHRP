// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 11, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or https://www.pnnl.gov/sysbio/ or https://panomics.pnnl.gov/
// -------------------------------------------------------------------------------
//
// Licensed under the 2-Clause BSD License; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// https://opensource.org/licenses/BSD-2-Clause
//
// Copyright 2018 Battelle Memorial Institute

using System;
using System.Collections.Generic;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class is used to track peptide sequences and their modification descriptors
    /// It assigns a unique integer ID to each combination of sequence and modification description
    /// </summary>
    public class clsUniqueSequencesContainer
    {
        #region "Constants and Enums"
        private const int DEFAULT_INITIAL_SEQ_ID = 1;
        private const char SEQUENCE_MOD_DESC_SEP = '_';
        #endregion

        #region "Class wide Variables"

        private readonly Dictionary<string, int> mMasterSequences;

        private int mNextUniqueSeqID;

        #endregion

        #region "Properties"

        // ReSharper disable once UnusedMember.Global
        public int UniqueSequenceCount => mMasterSequences.Count;

        #endregion

        /// <summary>
        /// Constructor
        /// </summary>
        public clsUniqueSequencesContainer()
        {
            mMasterSequences = new Dictionary<string, int>();
            Clear();
        }

        public void Clear()
        {
            Clear(DEFAULT_INITIAL_SEQ_ID);
        }

        public void Clear(int initialSeqID)
        {
            // Clears mMasterSequences and resets mNextUniqueSeqID to initialSeqID
            mMasterSequences.Clear();
            mNextUniqueSeqID = initialSeqID;
        }

        public int GetNextUniqueSequenceID(string sequence, string modDescription, out bool existingSequenceFound)
        {
            int uniqueSeqID;

            existingSequenceFound = false;

            try
            {
                if (sequence == null)
                    sequence = string.Empty;
                if (modDescription == null)
                    modDescription = string.Empty;

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
