// This file is used by Code Analysis to maintain SuppressMessage
// attributes that are applied to this project.
// Project-level suppressions either have no target or are given
// a specific target and scoped to a namespace, type, member, etc.

using System.Diagnostics.CodeAnalysis;

[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.clsPeptideCleavageStateCalculator.InitializeRegExObjects")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.clsPeptideMassCalculator.ComputeSequenceMass(System.String)~System.Double")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.clsPeptideMassCalculator.ConvertAminoAcidSequenceToEmpiricalFormula(System.String)~PHRPReader.clsEmpiricalFormula")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.clsPHRPParser.IsNumber(System.String)~System.Boolean")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.clsPHRPReader.IsNumber(System.String)~System.Boolean")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.clsPeptideModificationContainer.LookupMassCorrectionTagByMass(System.Double,System.Byte,System.Boolean,System.Byte)~System.String")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.clsPeptideModificationContainer.SetDefaultMassCorrectionTags")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.clsPeptideModificationContainer.StoreMassCorrectionTag(System.String,System.Double)")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.clsPeptideModificationContainer.UpdateDefaultModificationSymbols(System.String)")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.clsPeptideModificationContainer.ValidateModificationsVsDefaultModificationSymbols")]
[assembly: SuppressMessage("Readability", "RCS1123:Add parentheses when necessary.", Justification = "Parentheses not needed", Scope = "member", Target = "~M:PHRPReader.clsPeptideMassCalculator.ConvoluteMass(System.Double,System.Int32,System.Int32,System.Double)~System.Double")]
[assembly: SuppressMessage("Readability", "RCS1234:Duplicate enum value.", Justification = "<Pending>", Scope = "type", Target = "~T:PHRPReader.clsPHRPReader.ePeptideHitResultType")]
[assembly: SuppressMessage("Readability", "RCS1234:Duplicate enum value.", Justification = "<Pending>", Scope = "type", Target = "~T:PHRPReader.clsPHRPReader.PeptideHitResultTypes")]
