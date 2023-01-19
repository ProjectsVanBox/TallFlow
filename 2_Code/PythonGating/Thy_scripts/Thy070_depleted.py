#!/usr/bin/env python

###Input
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus/export_A8 Thy070_depleted_stained_CD45+.csv"
OutputFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/Thymus/Thy070_depleted_keymarkers.csv"

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
			if ((float(SL[14]) >= 1200.0)):
				outf.write(str(Index)+","+line+","+"B-cell"+"\n")
			elif ((float(SL[15]) >= 237.0)):
				outf.write(str(Index)+","+line+","+"Monocyte"+"\n")
			elif ((float(SL[10]) >= 429.0)):
				outf.write(str(Index)+","+line+","+"NK-cell"+"\n")
			elif ((float(SL[22]) >= 321.0)):
				outf.write(str(Index)+","+line+","+"cDC"+"\n")
			elif ((float(SL[19]) >= 1100.0)):
				outf.write(str(Index)+","+line+","+"pDC"+"\n")
			elif ((float(SL[11]) >= 476.0) & (float(SL[16]) >= 373.0)):
				outf.write(str(Index)+","+line+","+"DP"+"\n")
			elif ((float(SL[11]) >= 476.0) & (float(SL[21]) >= 414.0) & (float(SL[16]) < 373.0)):
				outf.write(str(Index)+","+line+","+"CD4"+"\n")
			elif ((float(SL[16]) >= 373.0)):
				outf.write(str(Index)+","+line+","+"CD8"+"\n")
			elif ((float(SL[11]) >= 476.0) & (float(SL[21]) < 414.0)):
				outf.write(str(Index)+","+line+","+"ISP"+"\n")
			elif ((float(SL[17]) >= 184.0) & (float(SL[27]) >= 510.0)):
				outf.write(str(Index)+","+line+","+"gdTcell"+"\n")
			elif ((float(SL[17]) >= 184.0)):
				outf.write(str(Index)+","+line+","+"DN"+"\n")
			else:
				outf.write(str(Index)+","+line+","+"Unknown"+"\n")
	# Close files
	outf.close()
	f.close()
	
CreateLabelFile(FlowFile, OutputFile) 
