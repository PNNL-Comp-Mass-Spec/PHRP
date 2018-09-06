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
using System.Collections;

namespace PeptideHitResultsProcessor
{
    /// <summary>
    /// This class is used to track peptide sequences and their modification descriptors
    /// It assigns a unique integer ID to each combination of sequence and modification description
    /// </summary>
    public class clsUniqueSequencesContainer
    {
        #region "Constants and Enums"
        protected const int DEFAULT_INITIAL_SEQ_ID = 1;
        protected const char SEQUENCE_MOD_DESC_SEP = '_';
        #endregion

        #region "Classwide Variables"
        protected Hashtable htMasterSequences;
        protected int mNextUniqueSeqID;
        #endregion

        #region "Properties"
        public int UniqueSequenceCount => htMasterSequences.Count;

        #endregion

        public clsUniqueSequencesContainer()
        {
            InitializeLocalVariables();
        }

        public void Clear()
        {
            this.Clear(DEFAULT_INITIAL_SEQ_ID);
        }

        public void Clear(int initialSeqID)
        {
            // Clears htMasterSequences and resets mNextUniqueSeqID to initialSeqID
            mNextUniqueSeqID = initialSeqID;
            if (htMasterSequences == null)
            {
                htMasterSequences = new Hashtable();
            }
            else
            {
                htMasterSequences.Clear();
            }
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

                if (htMasterSequences.ContainsKey(key))
                {
                    uniqueSeqID = Convert.ToInt32(htMasterSequences[key]);
                    existingSequenceFound = true;
                }
                else
                {
                    htMasterSequences.Add(key, mNextUniqueSeqID);
                    uniqueSeqID = mNextUniqueSeqID;
                    mNextUniqueSeqID += 1;
                }
            }
            catch (Exception)
            {
                uniqueSeqID = 0;
            }

            return uniqueSeqID;
        }

        protected void InitializeLocalVariables()
        {
            this.Clear();
        }
    }
}
