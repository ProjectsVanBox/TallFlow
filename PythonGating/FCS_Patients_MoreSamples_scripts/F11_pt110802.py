#!/usr/bin/env python

###Input
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Patients_MoreSamples/export_F11 110802_stained_CD45+.csv"
OutputFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/FCS_Patients_MoreSamples/F11_pt110802_keymarkers.csv"

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
			if ((float(SL[23]) >= 1200.0)):
				outf.write(str(Index)+","+line+","+"B-cell-like"+"\n")
			elif ((float(SL[15]) >= 2700.0)):
				outf.write(str(Index)+","+line+","+"Monocyte-like"+"\n")
			elif ((float(SL[10]) >= 3000.0)):
				outf.write(str(Index)+","+line+","+"NK-cell-like"+"\n")
			elif ((float(SL[11]) >= 2900.0)):
				outf.write(str(Index)+","+line+","+"cDC-like"+"\n")
			elif ((float(SL[19]) >= 1600.0)):
				outf.write(str(Index)+","+line+","+"pDC-like"+"\n")
			elif ((float(SL[22]) >= 393.0) & (float(SL[16]) >= 468.0)):
				outf.write(str(Index)+","+line+","+"DP-like"+"\n")
			elif ((float(SL[22]) >= 393.0) & (float(SL[24]) >= 11000.0) & (float(SL[13]) < 400.0)):
				outf.write(str(Index)+","+line+","+"ISP-like"+"\n")
			elif ((float(SL[22]) >= 393.0) & (float(SL[16]) < 468.0)):
				outf.write(str(Index)+","+line+","+"CD4-like"+"\n")
			elif ((float(SL[16]) >= 468.0)):
				outf.write(str(Index)+","+line+","+"CD8-like"+"\n")
			elif ((float(SL[17]) >= 483.0) & (float(SL[27]) >= 2800.0)):
				outf.write(str(Index)+","+line+","+"gdTcell-like"+"\n")
			elif ((float(SL[24]) >= 11000.0) & (float(SL[17]) >= 483.0)):
				outf.write(str(Index)+","+line+","+"DN3-like"+"\n")
			elif ((float(SL[17]) >= 483.0)):
				outf.write(str(Index)+","+line+","+"DN-like"+"\n")
			else:
				outf.write(str(Index)+","+line+","+"Unknown"+"\n")
	# Close files
	outf.close()
	f.close()
	
CreateLabelFile(FlowFile, OutputFile) 
