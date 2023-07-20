#!/usr/bin/env python

###Input
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus_MoreSamples/export_A4 Thy083_depleted_stained_CD45+.csv"
OutputFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/ThymusMore/A4_ptThy083_depleted_keymarkers.csv"

def CreateLabelFile(InFile, OutFile):
	f = open(InFile, "r")
	outf = open(OutFile, "w")
	FirstLine = True
	Index = 0

	for line in f.readlines():
		line = line.rstrip()
		SL = line.split(",") # SplitLine

		if FirstLine:
			outf.write("Index"+","+line+","+"Celltype"+"\n")
			FirstLine = False
		else:
			Index = Index + 1
			if ((float(SL[14]) >= 2300.0)):
				outf.write(str(Index)+","+line+","+"B-cell"+"\n")
			elif ((float(SL[15]) >= 992.0)):
				outf.write(str(Index)+","+line+","+"Monocyte"+"\n")
			elif ((float(SL[10]) >= 1100.0)):
				outf.write(str(Index)+","+line+","+"NK-cell"+"\n")
			elif ((float(SL[22]) >= 1200.0)):
				outf.write(str(Index)+","+line+","+"cDC"+"\n")
			elif ((float(SL[19]) >= 1700.0)):
				outf.write(str(Index)+","+line+","+"pDC"+"\n")
			elif ((float(SL[11]) >= 726.0) & (float(SL[16]) >= 212.0)):
				outf.write(str(Index)+","+line+","+"DP"+"\n")
			elif ((float(SL[11]) >= 726.0) & (float(SL[21]) >= 626.0) & (float(SL[16]) < 212.0)):
				outf.write(str(Index)+","+line+","+"CD4"+"\n")
			elif ((float(SL[16]) >= 212.0)):
				outf.write(str(Index)+","+line+","+"CD8"+"\n")
			elif ((float(SL[11]) >= 726.0) & (float(SL[21]) < 626.0)):
				outf.write(str(Index)+","+line+","+"ISP"+"\n")
			elif ((float(SL[17]) >= 111.0) & (float(SL[27]) >= 456.0)):
				outf.write(str(Index)+","+line+","+"gdTcell"+"\n")
			elif ((float(SL[17]) >= 111.0)):
				outf.write(str(Index)+","+line+","+"DN"+"\n")
			else:
				outf.write(str(Index)+","+line+","+"Unknown"+"\n")
	# Close files
	outf.close()
	f.close()
	
CreateLabelFile(FlowFile, OutputFile) 
