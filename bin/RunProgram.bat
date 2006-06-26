@echo off

rem SELECT AJ.AJ_jobID, AJ.AJ_start, DS.Dataset_Num, 
rem     AJR.[Results Folder Path]
rem FROM dbo.T_Param_Files PF INNER JOIN
rem     dbo.T_Analysis_Job AJ ON 
rem     PF.Param_File_Name = AJ.AJ_parmFileName INNER JOIN
rem     dbo.T_Dataset DS ON 
rem     AJ.AJ_datasetID = DS.Dataset_ID INNER JOIN
rem     dbo.V_Analysis_Job_ReportEx AJR ON 
rem     AJ.AJ_jobID = AJR.JobNum
rem WHERE (PF.Param_File_Description LIKE '%C[_]Term[_]P%') AND 
rem     (AJ.AJ_start >= '5/28/2006') OR
rem     (PF.Param_File_Description LIKE '%C[_]Term[_]P%') AND 
rem     (AJ.AJ_StateID < 4)
rem ORDER BY AJ.AJ_jobID

rem PeptideHitResultsProcRunner.exe /i:IMAC_HSP90_H304Q_061406_syn.txt /m:sequest_N14_NE_STY_Phos_Stat_Met_Iodo_Tryp_ModDefs.txt /t:Mass_Correction_Tags.txt
rem PeptideHitResultsProcRunner.exe /i:SR_Control_screening_110505_syn.txt /m:sequest_N14_Tryp_STY_Phos_Stat_Met_ModDefs.txt /t:Mass_Correction_Tags.txt