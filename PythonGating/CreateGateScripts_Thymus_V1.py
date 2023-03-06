#!/usr/bin/env python


### Read threshold file and return corresponding threshold per marker 
# Note that the names should be the exact same! 
def RetrieveThresholds(InFile):
	f = open(InFile, "r")
	ThreshDict = {}

	for line in f.readlines():
		line = line.rstrip()
		if line.startswith("#"):
			pass
		else:
			SL = line.split("\t") # SplitLine
			ThreshDict[SL[0]] = float(SL[1])

	# Close and return 
	f.close()
	return ThreshDict


### Read the fcs file and get the on which spot (splitLine) is which marker
def GetMarkerSpots(InFile):
	MarkerSpot = {}

	with open(InFile) as f:
		line = f.readline()
		line = line.rstrip()
		SL = line.split(",") # SplitLine

		#if line.startswith("Index"):
		index = 0 
		for marker in SL:
			if ("__" in marker):
				markerSplit = marker.split("__")[1]
				markerSplit = markerSplit.replace('-A"','')
				#print(markerSplit)
				if (markerSplit in ThreshDict):
					MarkerSpot[markerSplit] = index
			index = index + 1 
		#print MarkerSpot
		f.close()
	return MarkerSpot


### Function to create a script that automatically can label your cells to the specified files 
# Infile = flow file containing the information of fluorescence of each cell, should be in csv format 
# Threshdict = dictionairy with all thresholds for each marker 
# OutFile = provided output file that will be labelled 
# OutputScript = location of where the script is going to be written to
# LabelFile = File containing the label (to be added), the positive markers and the negative markers for that cell type 
def CreateScript(InFile, ThreshDict, MarkerSpot, OutFile, OutputScript, LabelFile):
	### Create the input variables 
	OutputString = "#!/usr/bin/env python\n\n###Input\n"
	OutputString += "FlowFile = "+'"'+str(InFile)+'"\n'
	OutputString += "OutputFile = "+'"'+str(OutFile)+'"\n'
	OutputString += "\n"

	### Create the code 
	OutputString += '''def CreateLabelFile(InFile, OutFile):
	f = open(InFile, "r")
	outf = open(OutFile, "w")
	FirstLine = True
	Index = 0

	for line in f.readlines():
		line = line.rstrip()
		SL = line.split(",") # SplitLine

		if FirstLine:
			outf.write("Index"+","+line+","+"Celltype"+"\\n")
			FirstLine = False
		else:
			Index = Index + 1'''
	OutputString += "\n"

	### Create statements from the colour file 
	FirstLine = True
	f = open(LabelFile, "r")
	for line in f.readlines():
		line = line.rstrip("\n")
		if line.startswith("#"):
			pass
		else:
			SL = line.split("\t")
			CellType = SL[0]

			# Start the statement 
			if (FirstLine):
				Statement = "\t\t\tif ("
				FirstLine = False
			else:
				Statement = "\t\t\telif ("

			# Fill statement with postive markers 
			PosCount = 0
			if (SL[1] != ""):
				PosSplit = SL[1].split(",")
				for PosMarker in PosSplit:
					PosCount += 1
					if (PosCount == len(PosSplit)):
						Statement += "(float(SL["+str(MarkerSpot[PosMarker])+"]) >= "+str(ThreshDict[PosMarker])+")"
					else: 
						Statement += "(float(SL["+str(MarkerSpot[PosMarker])+"]) >= "+str(ThreshDict[PosMarker])+") & "

			# Fill statement with negative markers
			NegCount = 0
			if (SL[2] != ""):
				NegSplit = SL[2].split(",")
				if (PosCount != 0):
					Statement += " & "
				for NegMarker in NegSplit:
					NegCount += 1 
					if (NegCount == len(NegSplit)):
						Statement += "(float(SL["+str(MarkerSpot[NegMarker])+"]) < "+str(ThreshDict[NegMarker])+")"
					else: 
						Statement += "(float(SL["+str(MarkerSpot[NegMarker])+"]) < "+str(ThreshDict[NegMarker])+") & "

			# End the statement 
			Statement += "):\n"
			Statement += '\t\t\t\toutf.write(str(Index)+","+line+","+"'+CellType+'"+"\\n")\n'
			OutputString += Statement
	f.close()
	OutputString += '\t\t\telse:\n\t\t\t\toutf.write(str(Index)+","+line+","+"Unknown"+"\\n")'

	### Close all variables 
	OutputString += '''
	# Close files
	outf.close()
	f.close()
	'''
	OutputString += "\n"

	### Create Call statement
	OutputString += "CreateLabelFile(FlowFile, OutputFile) "
	OutputString += "\n"

	### Create the script
	OutScript = open(OutputScript, "w")
	OutScript.write(OutputString)
	OutScript.close()



### Input files 
ThresholdFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/ThresholdFiles/MarkerThresholds_transposed_070.txt"
LabelCodeFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/20220801_GatingStrategy.txt"
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus/export_A6 Thy070_undepleted_stained_CD45+.csv"
FlowOutfile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/Thymus/Thy070_undepleted_keymarkers.csv"
ScriptFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/2_Code/PythonGating/Thy_scripts/Thy070_undepleted.py"
### Call all functions 
ThreshDict = RetrieveThresholds(ThresholdFile)
MarkerSpot = GetMarkerSpots(FlowFile)
CreateScript(FlowFile, ThreshDict, MarkerSpot, FlowOutfile, ScriptFile, LabelCodeFile)


### Input files 
ThresholdFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/ThresholdFiles/MarkerThresholds_transposed_070.txt"
LabelCodeFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/20220801_GatingStrategy.txt"
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus/export_A8 Thy070_depleted_stained_CD45+.csv"
FlowOutfile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/Thymus/Thy070_depleted_keymarkers.csv"
ScriptFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/2_Code/PythonGating/Thy_scripts/Thy070_depleted.py"
### Call all functions 
ThreshDict = RetrieveThresholds(ThresholdFile)
MarkerSpot = GetMarkerSpots(FlowFile)
CreateScript(FlowFile, ThreshDict, MarkerSpot, FlowOutfile, ScriptFile, LabelCodeFile)


### Input files 
ThresholdFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/ThresholdFiles//MarkerThresholds_transposed_072.txt"
LabelCodeFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/20220801_GatingStrategy.txt"
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus/export_A2 Thy072_undepleted_stained_CD45+.csv"
FlowOutfile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/Thymus/Thy072_undepleted_keymarkers.csv"
ScriptFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/2_Code/PythonGating/Thy_scripts/Thy072_undepleted.py"
ThreshDict = RetrieveThresholds(ThresholdFile)
MarkerSpot = GetMarkerSpots(FlowFile)
CreateScript(FlowFile, ThreshDict, MarkerSpot, FlowOutfile, ScriptFile, LabelCodeFile)


### Input files 
ThresholdFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/ThresholdFiles//MarkerThresholds_transposed_072.txt"
LabelCodeFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/20220801_GatingStrategy.txt"
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus/export_A6 Thy072_70Mdepleted_stained_CD45+.csv"
FlowOutfile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/Thymus/Thy072_depleted_keymarkers.csv"
ScriptFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/2_Code/PythonGating/Thy_scripts/Thy072_depleted.py"
ThreshDict = RetrieveThresholds(ThresholdFile)
MarkerSpot = GetMarkerSpots(FlowFile)
CreateScript(FlowFile, ThreshDict, MarkerSpot, FlowOutfile, ScriptFile, LabelCodeFile)


### Input files 
ThresholdFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/ThresholdFiles/MarkerThresholds_transposed_073.txt"
LabelCodeFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/20220801_GatingStrategy.txt"
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus/export_A6 Thy073_10Mand70M_depleted_stained_CD45+.csv"
FlowOutfile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/Thymus/Thy073_depleted_keymarkers.csv"
ScriptFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/2_Code/PythonGating/Thy_scripts/Thy073_depleted.py"
### Call all functions 
ThreshDict = RetrieveThresholds(ThresholdFile)
MarkerSpot = GetMarkerSpots(FlowFile)
CreateScript(FlowFile, ThreshDict, MarkerSpot, FlowOutfile, ScriptFile, LabelCodeFile)


### Input files 
ThresholdFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/ThresholdFiles/MarkerThresholds_transposed_073.txt"
LabelCodeFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/20220801_GatingStrategy.txt"
FlowFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/1_Input/FCS_Thymus/export_export_A2 Thy073_undepleted_stained_CD45+.csv"
FlowOutfile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/3_Output/PythonGating/Thymus/Thy073_undepleted_keymarkers.csv"
ScriptFile = "/Users/ricohagelaar/Documents/Thymus_Project/FlowCytometry/TallFlow/2_Code/PythonGating/Thy_scripts/Thy073_undepleted.py"
### Call all functions 
ThreshDict = RetrieveThresholds(ThresholdFile)
MarkerSpot = GetMarkerSpots(FlowFile)
CreateScript(FlowFile, ThreshDict, MarkerSpot, FlowOutfile, ScriptFile, LabelCodeFile)


#for i in $( ls *.py ); do python $i; done







