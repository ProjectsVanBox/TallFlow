#!/usr/bin/env python

###Input
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus/export_export_A2 Thy073_undepleted_stained_CD45+.csv"
OutputFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/Thymus/Thy073_undepleted_keymarkers.csv"

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
			elif ((float(SL[15]) >= 658.0)):
				outf.write(str(Index)+","+line+","+"Monocyte"+"\n")
			elif ((float(SL[10]) >= 379.0)):
				outf.write(str(Index)+","+line+","+"NK-cell"+"\n")
			elif ((float(SL[22]) >= 1500.0)):
				outf.write(str(Index)+","+line+","+"cDC"+"\n")
			elif ((float(SL[19]) >= 1000.0)):
				outf.write(str(Index)+","+line+","+"pDC"+"\n")
			elif ((float(SL[11]) >= 633.0) & (float(SL[16]) >= 1000.0)):
				outf.write(str(Index)+","+line+","+"DP"+"\n")
			elif ((float(SL[11]) >= 633.0) & (float(SL[21]) >= 711.0) & (float(SL[16]) < 1000.0)):
				outf.write(str(Index)+","+line+","+"CD4"+"\n")
			elif ((float(SL[16]) >= 1000.0)):
				outf.write(str(Index)+","+line+","+"CD8"+"\n")
			elif ((float(SL[11]) >= 633.0) & (float(SL[21]) < 711.0)):
				outf.write(str(Index)+","+line+","+"ISP"+"\n")
			elif ((float(SL[17]) >= 195.0) & (float(SL[27]) >= 386.0)):
				outf.write(str(Index)+","+line+","+"gdTcell"+"\n")
			elif ((float(SL[17]) >= 195.0)):
				outf.write(str(Index)+","+line+","+"DN"+"\n")
			else:
				outf.write(str(Index)+","+line+","+"Unknown"+"\n")
	# Close files
	outf.close()
	f.close()
	
CreateLabelFile(FlowFile, OutputFile) 
