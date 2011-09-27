The Peptide Hit Results Processor can be used to convert an XTandem results 
file (XML format), Sequest Synopsis/First Hits file, Inspect search result 
file, or MSGF-DB search result file to a series of tab-delimited text files 
summarizing the results. It will insert modification symbols into the 
peptide sequences for modified peptides.  Parallel files are created containing 
sequence information, modification details, and protein information.  The 
user can optionally provide a modification definition file that specifies 
the symbol to use for each modification mass.

Example usage:
PeptideHitResultsProcRunner.exe /i:ExampleXTandemData_xt.xml /m:Inspect_NoMods_ModDefs.txt


Program syntax:
PeptideHitResultsProcRunner.exe InputFilePath [/O:OutputFolderPath]
 [/P:ParameterFilePath] [/M:ModificationDefinitionFilePath]
 [/T:MassCorrectionTagsFilePath] [/N:SearchToolParameterFilePath] [/SynPvalue:0.2]
 [/InsFHT:True|False] [/InsSyn:True|False]
 [/S:[MaxLevel]] [/A:AlternateOutputFolderPath] [/R] [/L:[LogFilePath]] [/Q]

The input file should be an XTandem Results file (_xt.xml), a Sequest Synopsis 
File (_syn.txt), a Sequest First Hits file (_fht.txt), an Inspect results file 
(_inspect.txt), or a MSGF-DB results file (_msgfdb.txt)  

The output folder switch is optional.  If omitted, the output file will be 
created in the same folder as the input file.

The parameter file path is optional.  If included, it should point to a 
valid XML parameter file.

Use /M to specify the file containing the modification definitions.  This file 
should be tab delimited, with the first column containing the modification 
symbol, the second column containing the modification mass, plus optionally a 
third column listing the residues that can be modified with the given mass 
(1 letter residue symbols, no need to separated with commas or spaces).

Use /T to specify the file containing the mass correction tag info.  This file 
should be tab delimited, with the first column containing the mass correction 
tag name and the second column containing the mass (the name cannot contain 
commas or colons and can be, at most, 8 characters long).

Use /N to specify the parameter file provided to the search tool.  This is 
only used when processing Inspect or MSGF-DB files.

When processing an Inspect results file, use /SynPvalue to customize the 
PValue threshold used to determine which peptides are written to the the 
synopsis file.  The default is /SynPvalue:0.2  Note that peptides with 
a TotalPRMScore >= 50 or an FScore >= 0 will also be included in the 
synopsis file.

Use /InsFHT:True or /InsFHT:False to toggle the creation of a first-hits 
file (_fht.txt) when processing Inspect or MSGF-DB results (default is /InsFHT:True)

Use /InsSyn:True or /InsSyn:False to toggle the creation of a synopsis 
file (_syn.txt) when processing Inspect or MSGF-DB results (default is /InsSyn:True)

Use /S to process all valid files in the input folder and subfolders. 
Include a number after /S (like /S:2) to limit the level of subfolders to examine.
When using /S, you can redirect the output of the results using /A.
When using /S, you can use /R to re-create the input folder hierarchy 
in the alternate output folder (if defined).

Use /L to specify that a log file should be created.  
Use /L:LogFilePath to specify the name (or full path) for the log file.

Use the optional /Q switch will suppress all error messages.

-------------------------------------------------------------------------------
Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA)
Copyright 2006, Battelle Memorial Institute.  All Rights Reserved.

E-mail: matthew.monroe@pnl.gov or matt@alchemistmatt.com
Website: http://ncrr.pnl.gov/ or http://www.sysbio.org/resources/staff/
-------------------------------------------------------------------------------

Licensed under the Apache License, Version 2.0; you may not use this file except 
in compliance with the License.  You may obtain a copy of the License at 
http://www.apache.org/licenses/LICENSE-2.0

All publications that result from the use of this software should include 
the following acknowledgment statement:
 Portions of this research were supported by the W.R. Wiley Environmental 
 Molecular Science Laboratory, a national scientific user facility sponsored 
 by the U.S. Department of Energy's Office of Biological and Environmental 
 Research and located at PNNL.  PNNL is operated by Battelle Memorial Institute 
 for the U.S. Department of Energy under contract DE-AC05-76RL0 1830.

Notice: This computer software was prepared by Battelle Memorial Institute, 
hereinafter the Contractor, under Contract No. DE-AC05-76RL0 1830 with the 
Department of Energy (DOE).  All rights in the computer software are reserved 
by DOE on behalf of the United States Government and the Contractor as 
provided in the Contract.  NEITHER THE GOVERNMENT NOR THE CONTRACTOR MAKES ANY 
WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS 
SOFTWARE.  This notice including this sentence must appear on any copies of 
this computer software.
