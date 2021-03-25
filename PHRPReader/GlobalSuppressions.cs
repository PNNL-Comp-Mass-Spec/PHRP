// This file is used by Code Analysis to maintain SuppressMessage
// attributes that are applied to this project.
// Project-level suppressions either have no target or are given
// a specific target and scoped to a namespace, type, member, etc.

using System.Diagnostics.CodeAnalysis;

[assembly: SuppressMessage("CodeQuality", "IDE0052:Remove unread private members", Justification = "Required for IDisposable", Scope = "member", Target = "~F:PHRPReader.clsPHRPReader.mDisposedValue")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.Data.PeptideModificationContainer.LookupMassCorrectionTagByMass(System.Double,System.Byte,System.Boolean,System.Byte)~System.String")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.Data.PeptideModificationContainer.SetDefaultMassCorrectionTags")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.Data.PeptideModificationContainer.StoreMassCorrectionTag(System.String,System.Double)")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.Data.PeptideModificationContainer.UpdateDefaultModificationSymbols(System.String)")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.Data.PeptideModificationContainer.ValidateModificationsVsDefaultModificationSymbols")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.PeptideCleavageStateCalculator.InitializeRegExObjects")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.PeptideMassCalculator.ComputeSequenceMass(System.String)~System.Double")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.PeptideMassCalculator.ConvertAminoAcidSequenceToEmpiricalFormula(System.String)~PHRPReader.Data.EmpiricalFormula")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.Reader.SynFileReaderBaseClass.IsNumber(System.String)~System.Boolean")]
[assembly: SuppressMessage("Design", "RCS1075:Avoid empty catch clause that catches System.Exception.", Justification = "Allowed", Scope = "member", Target = "~M:PHRPReader.ReaderFactory.IsNumber(System.String)~System.Boolean")]
