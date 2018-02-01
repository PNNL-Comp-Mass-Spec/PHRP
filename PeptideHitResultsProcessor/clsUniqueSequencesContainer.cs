// This class is used to track peptide sequences and their modification descriptors
// It assigns a unique integer ID to each combination of sequence and modification description
//
// -------------------------------------------------------------------------------
// Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
// Program started January 11, 2006
//
// E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov
// Website: https://omics.pnl.gov/ or http://www.sysbio.org/resources/staff/
// -------------------------------------------------------------------------------
//
// Licensed under the Apache License, Version 2.0; you may not use this file except
// in compliance with the License.  You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Notice: This computer software was prepared by Battelle Memorial Institute,
// hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the
// Department of Energy (DOE).  All rights in the computer software are reserved
// by DOE on behalf of the United States Government and the Contractor as
// provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY
// WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
// SOFTWARE.  This notice including this sentence must appear on any copies of
// this computer software.

using System;
using System.Collections;

namespace PeptideHitResultsProcessor
{
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
