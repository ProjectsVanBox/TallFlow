#!/usr/bin/env python

###Input
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/ThyCITE/export_A2 thy03_undepleted_stained_CD45+.csv"
OutputFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/ThyCITE/A2_ptthy03_undepleted_keymarkers.csv"

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
			if ((float(SL[22]) >= 1100.0)):
				outf.write(str(Index)+","+line+","+"B-cell-like"+"\n")
			elif ((float(SL[14]) >= 836.0)):
				outf.write(str(Index)+","+line+","+"Monocyte-like"+"\n")
			elif ((float(SL[9]) >= 1800.0)):
				outf.write(str(Index)+","+line+","+"NK-cell-like"+"\n")
			elif ((float(SL[10]) >= 1200.0)):
				outf.write(str(Index)+","+line+","+"cDC-like"+"\n")
			elif ((float(SL[18]) >= 260.0)):
				outf.write(str(Index)+","+line+","+"pDC-like"+"\n")
			elif ((float(SL[21]) >= 1700.0) & (float(SL[15]) >= 470.0)):
				outf.write(str(Index)+","+line+","+"DP-like"+"\n")
			elif ((float(SL[21]) >= 1700.0) & (float(SL[23]) >= 1100.0) & (float(SL[12]) < 342.0)):
				outf.write(str(Index)+","+line+","+"ISP-like"+"\n")
			elif ((float(SL[21]) >= 1700.0) & (float(SL[15]) < 470.0)):
				outf.write(str(Index)+","+line+","+"CD4-like"+"\n")
			elif ((float(SL[15]) >= 470.0)):
				outf.write(str(Index)+","+line+","+"CD8-like"+"\n")
			elif ((float(SL[16]) >= 127.0) & (float(SL[26]) >= 526.0)):
				outf.write(str(Index)+","+line+","+"gdTcell-like"+"\n")
			elif ((float(SL[23]) >= 1100.0) & (float(SL[16]) >= 127.0)):
				outf.write(str(Index)+","+line+","+"DN3-like"+"\n")
			elif ((float(SL[16]) >= 127.0)):
				outf.write(str(Index)+","+line+","+"DN-like"+"\n")
			else:
				outf.write(str(Index)+","+line+","+"Unknown"+"\n")
	# Close files
	outf.close()
	f.close()
	
CreateLabelFile(FlowFile, OutputFile) 
