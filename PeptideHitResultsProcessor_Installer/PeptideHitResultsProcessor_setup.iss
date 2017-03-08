; This is an Inno Setup configuration file
; http://www.jrsoftware.org/isinfo.php

#define ApplicationVersion GetFileVersion('..\bin\PeptideHitResultsProcRunner.exe')

[CustomMessages]
AppName=Peptide Hit Results Processor
[Messages]
WelcomeLabel2=This will install [name/ver] on your computer.
; Example with multiple lines:
; WelcomeLabel2=Welcome message%n%nAdditional sentence
[Files]
Source: ..\bin\PeptideHitResultsProcRunner.exe         ; DestDir: {app}
Source: ..\bin\PeptideHitResultsProcRunner.exe.config  ; DestDir: {app}
Source: ..\bin\PeptideHitResultsProcessor.dll          ; DestDir: {app}
Source: ..\bin\PeptideToProteinMapEngine.dll           ; DestDir: {app}
Source: ..\bin\PHRPReader.dll                          ; DestDir: {app}
Source: ..\bin\PeptideHitResultsProcessor.dll          ; DestDir: {app}
Source: ..\bin\ProteinFileReader.dll                   ; DestDir: {app}
Source: ..\bin\System.Data.SQLite.dll                  ; DestDir: {app}
Source: ..\bin\x64\SQLite.Interop.dll                  ; DestDir: {app}\x64
Source: ..\bin\x86\SQLite.Interop.dll                  ; DestDir: {app}\x86

Source: ..\bin\Convert_mzid_to_tsv.bat                        ; DestDir: {app}

Source: ..\bin\Example_ModDefs.txt                            ; DestDir: {app}
Source: ..\bin\ExampleData_syn.txt                            ; DestDir: {app}
Source: ..\bin\Mass_Correction_Tags.txt                       ; DestDir: {app}
Source: ..\bin\MSGFPlus_Example_Data.zip                      ; DestDir: {app}
Source: ..\bin\PeptideHitResultsProcessorParameters.xml       ; DestDir: {app}
Source: ..\bin\RunProgram.bat                                 ; DestDir: {app}
Source: Images\delete_16x.ico                ; DestDir: {app}
Source: ..\Readme.txt                        ; DestDir: {app}
Source: ..\RevisionHistory.txt               ; DestDir: {app}

[Dirs]
Name: {commonappdata}\PeptideHitResultsProcessor; Flags: uninsalwaysuninstall

[Tasks]
Name: desktopicon; Description: {cm:CreateDesktopIcon}; GroupDescription: {cm:AdditionalIcons}; Flags: unchecked
; Name: quicklaunchicon; Description: {cm:CreateQuickLaunchIcon}; GroupDescription: {cm:AdditionalIcons}; Flags: unchecked

[Icons]
Name: {commondesktop}\Peptide Hit Results Processor; Filename: {app}\PeptideHitResultsProcRunner.exe; Tasks: desktopicon; Comment: Peptide Hit Results ProcessorGUI
Name: {group}\Peptide Hit Results Processor; Filename: {app}\PeptideHitResultsProcRunner.exe; Comment: Peptide Hit Results Processor

[Setup]
AppName=Peptide Hit Results Processor
AppVersion={#ApplicationVersion}
;AppVerName=PeptideHitResultsProcessor
AppID=PeptideHitResultsProcessorId
AppPublisher=Pacific Northwest National Laboratory
AppPublisherURL=http://omics.pnl.gov/software
AppSupportURL=http://omics.pnl.gov/software
AppUpdatesURL=http://omics.pnl.gov/software
DefaultDirName={pf}\PeptideHitResultsProcessor
DefaultGroupName=PAST Toolkit
AppCopyright=© PNNL
;LicenseFile=.\License.rtf
PrivilegesRequired=poweruser
OutputBaseFilename=PeptideHitResultsProcessor_Installer
;VersionInfoVersion=1.57
VersionInfoVersion={#ApplicationVersion}
VersionInfoCompany=PNNL
VersionInfoDescription=Peptide Hit Results Processor
VersionInfoCopyright=PNNL
DisableFinishedPage=true
ShowLanguageDialog=no
ChangesAssociations=false
EnableDirDoesntExistWarning=false
AlwaysShowDirOnReadyPage=true
UninstallDisplayIcon={app}\delete_16x.ico
ShowTasksTreeLines=true
OutputDir=.\Output
[Registry]
;Root: HKCR; Subkey: MyAppFile; ValueType: string; ValueName: ; ValueDataMyApp File; Flags: uninsdeletekey
;Root: HKCR; Subkey: MyAppSetting\DefaultIcon; ValueType: string; ValueData: {app}\wand.ico,0; Flags: uninsdeletevalue
[UninstallDelete]
Name: {app}; Type: filesandordirs
