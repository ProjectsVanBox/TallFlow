#!/usr/bin/env python

###Input
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_manGate_Thymi/export_A10 Thy075_thawed_stained_CD45+.csv.csv"
OutputFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/ThymusManGated/A10_ptThy075_thawed_keymarkers.csv"

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
			if ((float(SL[15]) >= 2300.0)):
				outf.write(str(Index)+","+line+","+"B-cell"+"\n")
			elif ((float(SL[16]) >= 992.0)):
				outf.write(str(Index)+","+line+","+"Monocyte"+"\n")
			elif ((float(SL[11]) >= 1100.0)):
				outf.write(str(Index)+","+line+","+"NK-cell"+"\n")
			elif ((float(SL[23]) >= 1200.0)):
				outf.write(str(Index)+","+line+","+"cDC"+"\n")
			elif ((float(SL[20]) >= 1700.0)):
				outf.write(str(Index)+","+line+","+"pDC"+"\n")
			elif ((float(SL[12]) >= 726.0) & (float(SL[17]) >= 212.0)):
				outf.write(str(Index)+","+line+","+"DP"+"\n")
			elif ((float(SL[12]) >= 726.0) & (float(SL[14]) >= 336.0) & (float(SL[22]) < 626.0)):
				outf.write(str(Index)+","+line+","+"ISP"+"\n")
			elif ((float(SL[12]) >= 726.0) & (float(SL[17]) < 212.0)):
				outf.write(str(Index)+","+line+","+"CD4"+"\n")
			elif ((float(SL[17]) >= 212.0)):
				outf.write(str(Index)+","+line+","+"CD8"+"\n")
			elif ((float(SL[18]) >= 111.0) & (float(SL[28]) >= 456.0)):
				outf.write(str(Index)+","+line+","+"gdTcell"+"\n")
			elif ((float(SL[18]) >= 111.0)):
				outf.write(str(Index)+","+line+","+"DN"+"\n")
			else:
				outf.write(str(Index)+","+line+","+"Unknown"+"\n")
	# Close files
	outf.close()
	f.close()
	
CreateLabelFile(FlowFile, OutputFile) 
