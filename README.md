# Peptide Hit Results Processor

The Peptide Hit Results Processor (PHRP) can be used to convert an X!Tandem results 
file (XML format), MSGF+ .tsv search result file, Sequest Synopsis/First Hits file, 
or Inspect search result into a series of tab-delimited text files 
summarizing the results. It will insert modification symbols into the 
peptide sequences for modified peptides.  Parallel files are created containing 
sequence information, modification details, and protein information.  The 
user can optionally provide a modification definition file that specifies 
the symbol to use for each modification mass.

## Example Data

Example input and output files are in the Data folder:
* MSGFPlus_Example has MSGF+ results
* XTandem_Example has X!Tandem results

For [MSGF+](https://github.com/sangtaekim/msgfplus) results, prior to running PHRP, use the 
[Mzid-To-Tsv-Converter](https://github.com/PNNL-Comp-Mass-Spec/Mzid-To-Tsv-Converter)
to convert the .mzid file to a tab-delimited .tsv file.

## Example Command line 

PeptideHitResultsProcRunner.exe /I:Dataset_msgfplus.tsv /M:MSGFDB_PartTryp_MetOx_20ppmParTol_ModDefs.txt /N:MSGFDB_PartTryp_MetOx_20ppmParTol.txt /T:Mass_Correction_Tags.txt /L /ProteinMods /F:Shewanella_oneidensis_MR1_2010-04-22_Tryp_Pig_Bov.revCat.fasta

## Console Switches

PHRP is a console application, and must be run from the Windows command prompt.

```
PeptideHitResultsProcRunner.exe InputFilePath [/O:OutputFolderPath]
 [/P:ParameterFilePath] [/M:ModificationDefinitionFilePath]
 [/ProteinMods] [/F:FastaFilePath] 
 [/ProteinModsViaPHRP] [/IgnorePepToProtMapErrors]
 [/ProteinModsIncludeReversed] [/UseExistingPepToProteinMapFile]
 [/T:MassCorrectionTagsFilePath] [/N:SearchToolParameterFilePath] 
 [/SynPvalue:0.2] [/InsFHT:True|False] [/InsSyn:True|False]
 [/S:[MaxLevel]] [/A:AlternateOutputFolderPath] [/R] [/L:[LogFilePath]]
```

The input file should be an XTandem Results file (_xt.xml), a Sequest Synopsis File 
(_syn.txt), a Sequest First Hits file (_fht.txt), an Inspect results file (_inspect.txt), an 
MSGF-DB results file (_msgfdb.txt), an MSGF+ results file (_msgfdb.tsv or _msgfplus.tsv), or 
an MSAlign results files (_MSAlign_ResultTable.txt)

The output folder switch is optional.  If omitted, the output file will be created in the same 
folder as the input file.

The parameter file path is optional.  If included, it should point to a valid XML parameter 
file.

Use /M to specify the file containing the modification definitions.  This file should be tab 
delimited, with the first column containing the modification symbol, the second column 
containing the modification mass, plus optionally a third column listing the residues that can 
be modified with the given mass (1 letter residue symbols, no need to separated with commas or 
spaces).

Use /ProteinMods to indicate that the _ProteinMods.txt file should be created.  This requires 
that either an existing _PepToProtMapMTS.txt file exist, or that the Fasta file be defined 
using /F

Use /ProteinModsViaPHRP to indicate that InputFilePath specifies a valid PHRP data file and 
thus the PHRP data files should not be re-created; only the _ProteinMods.txt file should be 
created.  This requires that either an existing _PepToProtMapMTS.txt file exist, or that the 
Fasta file be defined using /F

Use /IgnorePepToProtMapErrors to ignore peptide to protein mapping errors that occur when 
creating a missing _PepToProtMapMTS.txt file

Use /ProteinModsIncludeReversed to include Reversed proteins in the _ProteinMods.txt file

Use /UseExistingPepToProteinMapFile to use an existing _PepToProtMapMTS.txt file if it exists

Use /T to specify the file containing the mass correction tag info.  This file should be tab 
delimited, with the first column containing the mass correction tag name and the second column 
containing the mass (the name cannot contain commas or colons and can be, at most, 8 
characters long).

Use /N to specify the parameter file provided to the search tool.  This is only used when 
processing Inspect or MSGF-DB files.

When processing an Inspect results file, use /SynPvalue to customize the PValue threshold used 
to determine which peptides are written to the the synopsis file.  The default is 
/SynPvalue:0.2  Note that peptides with a TotalPRMScore >= 50 or an FScore >= 0 will also be 
included in the synopsis file.

Use /InsFHT:True or /InsFHT:False to toggle the creation of a first-hits file (_fht.txt) when 
processing Inspect or MSGF-DB results (default is /InsFHT:True)

Use /InsSyn:True or /InsSyn:False to toggle the creation of a synopsis file (_syn.txt) when 
processing Inspect or MSGF-DB results (default is /InsSyn:True)

Use /S to process all valid files in the input folder and subfolders. Include a number after 
/S (like /S:2) to limit the level of subfolders to examine.

When using /S, you can redirect the output of the results using /A.

When using /S, you can use /R to re-create the input folder hierarchy in the alternate output 
folder (if defined).

Use /L to specify that a log file should be created.  Use /L:LogFilePath to specify the name 
(or full path) for the log file.

## Contacts

Written by Matthew Monroe for the Department of Energy (PNNL, Richland, WA) \
E-mail: matthew.monroe@pnnl.gov or proteomics@pnnl.gov\
Website: https://omics.pnl.gov/ or https://panomics.pnnl.gov/

## License

Peptide Hit Results Processor is licensed under the Apache License, Version 2.0; you may not use this 
file except in compliance with the License.  You may obtain a copy of the 
License at https://opensource.org/licenses/Apache-2.0