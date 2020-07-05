#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

// IMOD models converted using model2point are the input here.
// Note that naming of files is important.
// GFP_Tub_X_Y and GFP_X_Y will cause problems.
// Must be "UniqueAlnum"_X_Y for the code to run in current state.

////////////////////////////////////////////////////////////////////////
// Menu items
////////////////////////////////////////////////////////////////////////
Menu "Macros"
	Submenu	"Cell Shape"
		"Load IMOD Models...",  IMODModelAnalysis()
		"Start Over", CleanSlate()
	End
End

////////////////////////////////////////////////////////////////////////
// Master functions and wrappers
////////////////////////////////////////////////////////////////////////
Function IMODModelAnalysis()
	PreLoader()
End
// we stop beween these two. User input needed.
Function TheLoader()
	LoadIMODModels()
	ProcessAllModels()
	CollectAllMeasurements()
	ProcessAllConditions()
	ImageQuilt()
	MakeTheLayouts("p",5,3, rev = 1, saveIt = 0)
	MakeTheLayouts("quilt",1,2, rev = 1, saveIt = 0, orient = 1) // landscape
	TidyAndSave("p")
	TidyAndSave("quilt")
End

////////////////////////////////////////////////////////////////////////
// Main functions
////////////////////////////////////////////////////////////////////////
Function PreLoader()
	NewPath/O/Q/M="Please find disk folder" expDiskFolder
	if (V_flag != 0)
		DoAlert 0, "Disk folder error"
		Return -1
	endif
	
	PathInfo/S expDiskFolder
	String pathString = S_path
	Make/O/N=1/T pathWave = {pathString}	
	String FileList = IndexedFile(expDiskFolder,-1,".txt")
	TheDiviner(ReplaceString(".txt",FileList,""))
End

Function LoadIMODModels()
	WAVE/T/Z pathWave
	if(!WaveExists(pathWave))
		DoAlert 0, "Error"
		return -1
	endif
	String pathString = pathWave[0]
	String expDiskFolderName, expDataFolderName
	String FileList, ThisFile
	Variable FileLoop
	
	NewPath/O/Q ExpDiskFolder, pathString
	// we will get the list again rather than use the stored text wave
	FileList = IndexedFile(expDiskFolder,-1,".txt")
	Variable nFiles = ItemsInList(FileList)
	
	NewDataFolder/O/S root:data
	
	for (FileLoop = 0; FileLoop < nFiles; FileLoop += 1)
		ThisFile = StringFromList(FileLoop, FileList)
		expDataFolderName = ReplaceString(".txt",ThisFile,"")
		NewDataFolder/O/S $expDataFolderName
		LoadWave/A/J/D/O/K=1/V={" "," $",0,0}/L={0,0,0,1,0}/P=expDiskFolder ThisFile
		MakeObjectContourWaves()
		SetDataFolder root:data:
	endfor
End

Function MakeObjectContourWaves()
	Concatenate/O/KILL wavelist("wave*",";",""), matA
	WaveStats/Q/RMD=[][0] matA
	// pixel size is hard-coded here
	Variable pxSize = 0.227 // micron per pixel
	// Scale the coordinates to real values
	matA[][2,4] *= pxSize
	Variable nObjects = V_max + 1 // objects in IMOD is 1-based
	Variable nContours, contourVar
	String wName
	
	Variable i,j
	
	for (i = 0; i < nObjects; i += 1)
		MatrixOP/O filtObj = col(matA,0)
		filtObj[] = (filtObj[p] == i) ? matA[p][1] : NaN
		WaveTransform zapnans filtObj
		FindDuplicates/RN=uniqueContours filtObj
		nContours = numpnts(uniqueContours)
		// zero-indexed list of contours in this object
		for (j = 0; j < nContours; j += 1)
			contourVar = uniqueContours[j]
			// find the rows that correspond to each contour
			Duplicate/O/FREE matA,matB
			matB[][] = (matB[p][0] == i && matB[p][1] == contourVar) ? matB[p][q] : NaN
			MatrixOp/O/FREE xW = col(matB,2)
			MatrixOp/O/FREE yW = col(matB,3)
			// no need to take the z column here
			WaveTransform zapnans xW
			WaveTransform zapnans yW
			// Now make ObjectContour waves
			wName = "cell_" + num2str(i) + "_" + num2str(contourVar)
			Concatenate/O/NP=1 {xW,yW}, $wName
			// close contour if it's not already closed 
			Wave cW = $wName // contour wave
			if(cW[DimSize(cW,0)-1][0] != cW[0][0] || cW[DimSize(cW,0)-1][1] != cW[0][1])
				InsertPoints DimSize(cW,0), 1, cW
				cW[DimSize(cW,0)-1][] = cW[0][q]
			endif
			// delete if it's just 1-3 pixels?
			if(DimSize(cW,0) < 4)
				KillWaves cW
			endif
		endfor
	endfor
	KillWaves/Z filtObj,UniqueContours,MatA
End

// this function goes into each datafolder and run some code on the contours in there
Function ProcessAllModels()
	SetDataFolder root:data:	// relies on earlier load
	DFREF dfr = GetDataFolderDFR()
	String folderName
	Variable numDataFolders = CountObjectsDFR(dfr, 4)
	Variable nWaves
	
	Variable i
		
	for(i = 0; i < numDataFolders; i += 1)
		folderName = GetIndexedObjNameDFR(dfr, 4, i)
		SetDataFolder ":'" + folderName + "':"
		// Look at all outlines. There are centring/rotation options. Plots are called q_*
		CentreAndPlot(0)
		TakeMeasurements()
		SetDataFolder root:data:
	endfor
	SetDataFolder root:
//	MakeTheLayouts("q_",6,8)
End

///	@param	plotOpt Variable to determine if we make a plot or not
STATIC Function CentreAndPlot(plotOpt)
	Variable plotOpt
	String wList = WaveList("cell*",";","")
	Variable nWaves = ItemsInList(wList)
	String wName, w0Name
	String folderName = GetDataFolder(0)
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		w0Name = "root:data:" + folderName + ":" + wName
		Wave/Z w0 = $w0Name
		// the wave generated here is centred and/or rotated. It's called c_*
		Wave w1 = Centralise2DWave(w0,midpoint = 1) // centring/rotation option here
		if(plotOpt == 1)
			MakeCellggPlot(w0)
		endif
	endfor
End

STATIC Function MakeCellggPlot(w0, [plusMinus])
	Wave w0
	Variable plusMinus
	String plotName = "q_" + GetWavesDataFolder(w0,0) + "_" + NameOfWave(w0)
	KillWindow/Z $plotname
	Display/N=$plotName/HIDE=1 w0[][1]/TN=Outline0 vs w0[][0]
	ModifyGraph/W=$plotName rgb(Outline0)=(68*257,170*257,153*257)
	ModifyGraph/W=$plotName height={Aspect,1}
	// now work out how to display the image
	if(ParamIsDefault(plusMinus) == 1)
		Variable last2 = dimSize(w0,0)
		Variable xMinVal = WaveMin(w0,0,last2 - 1)
		Variable xMaxVal = WaveMax(w0,0,last2 - 1)
		xMaxVal = Max(xMaxVal,abs(xMinVal))
		Variable yMinVal = WaveMin(w0,last2,last2*2 - 1)
		Variable yMaxVal = WaveMax(w0,last2,last2*2 - 1)
		yMaxVal = Max(yMaxVal,abs(yMinVal))
		Variable theMaxIs = Max(xMaxVal,yMaxVal)
		SetAxis/W=$plotName left -theMaxIs,theMaxIs
		SetAxis/W=$plotName bottom theMaxIs,-theMaxIs
	else
		SetAxis/W=$plotName left -plusMinus,plusMinus
		SetAxis/W=$plotName bottom plusMinus,-plusMinus
	endif
	ModifyGraph/W=$plotName grid=1,gridRGB=(65535,65535,65535)
	ModifyGraph/W=$plotName axRGB=(65535,65535,65535),tlblRGB=(65535,65535,65535),alblRGB=(65535,65535,65535)
	ModifyGraph/W=$plotName noLabel=2,axThick=1,standoff=0
	ModifyGraph/W=$plotName margin=2
	ModifyGraph/W=$plotName gbRGB=(61166,61166,61166)
	TextBox/W=$plotName/C/N=text0/F=0/A=RB/X=0.00/Y=0.00 GetWavesDataFolder(w0,0)
End

// Centralise 2D Wave using centre of mass (default) or midpoint
STATIC Function/WAVE Centralise2DWave(w,[midpoint])
	Wave w
	Variable midpoint
	String newName = GetWavesDataFolder(w,1) + "c_" + NameOfWave(w)
	Duplicate/O w,$newName
	Wave wc = $newName
	if(ParamIsDefault(midpoint) == 1)
		// centre of mass
		MatrixOp/O/FREE meanMat = averageCols(w)
		wc[][0] -= meanMat[0][0]
		wc[][1] -= meanMat[0][1]
	elseif(midpoint == 0)
		// find midpoints in x and y
		Variable nRows = DimSize(w,0)
		Variable xMid = (WaveMax(w,0,nRows - 1) + WaveMin(w,0,nRows - 1) ) / 2
		Variable yMid = (WaveMax(w,nRows,nRows*2 - 1) + WaveMin(w,nRows,nRows*2 - 1) ) / 2
		wc[][0] -= xMid
		wc[][1] -= yMid
	elseif(midpoint == 1)
		// centre and rotate
		FindEV(wc)
	endif
	return wc
End

///	@param	m1	2D wave of xy coords
STATIC Function FindEV(m1)
	Wave m1
	MatrixOp/O/FREE xCoord = col(m1,0)
	MatrixOp/O/FREE yCoord = col(m1,1)
	
	// translate to origin
	Variable offX = mean(xCoord)
	Variable offY = mean(yCoord)
	xCoord[] -= offX
	yCoord[] -= offY
	// do PCA. Rotated points are in M_R
	PCA/ALL/SEVC/SRMT/SCMT xCoord,yCoord
	WAVE M_R
	m1[][] = M_R[p][q]
	KillWaves/Z M_R
End

Function TakeMeasurements()
	// full list of measurements:
	//  MinAxis, MajAxis, Perimeter, Area
	//  ConvexArea, Solidity, Extent
	//  AreaROI, Symmetry
	//  MaxFromCentre, MinFromCentre
	//  AspectRatio (aka Eccentricity), Circularity, Circularity2
	// Not computing irregularity measures
	String wList = WaveList("c_cell_*",";","")
	Variable nWaves = ItemsInList(wList)
	Make/O/N=(nWaves) Img_MinAxis, Img_MajAxis, Img_Perimeter, Img_Area
	Make/O/N=(nWaves) Img_ConvexArea, Img_Solidity, Img_Extent
	Make/O/N=(nWaves) Img_AreaROI, Img_Symmetry
	Make/O/N=(nWaves) Img_MaxFromCentre, Img_MinFromCentre
	Make/O/N=(nWaves,2)/T Img_Label
	String currentDF = GetDataFolder(0)
	Img_label[][0] = currentDF
	String wName
	Variable BoundingBoxArea
	
	Variable i,j
	
	for (i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Img_Label[i][1] = wName
		Wave w0 = $wName
		// MinAxis, MajAxis, Perimeter
		Img_MinAxis[i] = VesicleAxisLength(w0,1)
		Img_MajAxis[i] = VesicleAxisLength(w0,0)
		Img_Perimeter[i] = FindLengthOfXYCoords(w0)
		// Area
		MatrixOp/O/FREE w0c0 = col(w0,0)
		MatrixOp/O/FREE w0c1 = col(w0,1)
		Img_Area[i] = PolygonArea(w0c0,w0c1)
		// ConvexArea
		ConvexHull/C w0c0,w0c1
		WAVE/Z W_Xhull,W_Yhull
		Img_ConvexArea[i] = PolygonArea(W_Xhull,W_Yhull)
		KillWaves/Z W_Xhull, W_Yhull
		// Solidity
		Img_Solidity[i] = Img_Area[i] / Img_ConvexArea[i] // could do this by MatrixOp below
		// Extent
		BoundingBoxArea = (WaveMax(w0c0) - WaveMin(w0c0)) * (WaveMax(w0c1) - WaveMin(w0c1))
		Img_Extent[i] = Img_Area[i] / BoundingBoxArea
		// To calculate AreaROI and Symmetry we'll send the co-ords to a separate function
		Wave resultW = SymmetryCalculator(w0c0,w0c1)
		Img_AreaROI[i] = resultW[0]
		Img_Symmetry[i] = resultW[1] / resultW[0]
		// Max or Min from Centre
		MatrixOp/O/FREE distanceW = sqrt(sumrows(w0 * w0))
		Img_MaxFromCentre[i] = WaveMax(distanceW)
		Img_MinFromCentre[i] = WaveMin(distanceW)
	endfor
	if(numpnts(Img_Area) > 0)
		MatrixOp/O/NTHR=0 Img_AspectRatio = Img_MinAxis / Img_MajAxis
		MatrixOp/O/NTHR=0 Img_Circularity = (4 * pi * Img_Area) / (Img_Perimeter * Img_Perimeter)
		MatrixOp/O/NTHR=0 Img_Circularity2 = Img_Perimeter / (2 * sqrt(pi * Img_Area))
	else
		Print "No cells in", currentDF
	endif
End

///	@param	m1	2D wave of xy coords
///	@param	colNo	column number to use for search
STATIC Function VesicleAxisLength(m1,colNo)
	Wave m1
	Variable colNo
	MatrixOp/O/FREE m1c0 = col(m1,0)
	MatrixOp/O/FREE m1c1 = col(m1,1)
	Variable V_Value,len
	if(colNo == 1)
		FindLevel/Q/EDGE=1/P m1c0, 0
		len = abs(m1c1(V_LevelX))
		FindLevel/Q/EDGE=2/P m1c0, 0
		len += abs(m1c1(V_LevelX))
	else
		FindLevel/Q/EDGE=1/P m1c1, 0
		len = abs(m1c0(V_LevelX))
		FindLevel/Q/EDGE=2/P m1c1, 0
		len += abs(m1c0(V_LevelX))
	endif
	return len
End

///	@param	m1	2D wave of xy coords
STATIC Function FindLengthOfXYCoords(m1)
	Wave m1
	// make new 2D wave of xy coords
	Duplicate/O/FREE m1,tempDist
	// offset to zero
	tempDist[][0] -= m1[0][0]
	tempDist[][1] -= m1[0][1]
	// Differentiate, backward difference
	Differentiate/METH=2 tempDist
	// find norm, cumlative distance
	MatrixOp/O/FREE/NTHR=0 tempNorm = sqrt(sumRows(tempDist * tempDist))
	tempNorm[0] = 0 // first point is garbage
	// return the sum of distances
	return sum(tempNorm)
End

STATIC Function/WAVE SymmetryCalculator(xw,yw)
	Wave xw,yw // 1d column of x and y coords of centered/rotated cell shape
	Make/O/N=(2)/FREE resultW = 0
	// convert back to pixels (hard-coded)
	Variable pxSize = 0.227
	xw /= pxSize
	yw /= pxSize
	// these coords are centred at the origin, we need the furthest point
	Duplicate/O/FREE xw, xw1
	Duplicate/O/FREE yw, yw1
	xw1[] = abs(xw[p])
	yw1[] = abs(yw[p])
	// this will find the most extreme point and add 2 to it
	// we need to do this so that we can set the background seed and not encounter the ROI
	Variable xOff = WaveMax(xw1) + 2
	Variable yOff = WaveMax(yw1) + 2
	// offset a copy of original coordinate set
	xw1[] = xOff + xw[p]
	yw1[] = yOff + yw[p]

	// make a mask with value of 1
	Wave maskMat = GenerateTheMask(2 * xOff, 2 * yOff, xw1, yw1)

	// store integer representation of ROI
	resultW[0] = sum(maskMat)
	// The resulting image is ceil(requestedSize) - 1
	
	// now make the two mirror image images
	Duplicate/O/FREE maskMat, upMask, downMask
	// Do the mirroring. If even, the image has duplicated middle rows. If odd, the new image will not
	Variable nCol = DimSize(maskMat,1) // the height of the image
	if(mod(nCol,2) == 0)
		upMask[][nCol / 2,nCol - 1] = maskMat[p][nCol - q - 1]
		downMask[][0,nCol / 2 - 1] = maskMat[p][nCol - 1 - q]
	else
		upMask[][(nCol + 1) / 2 , nCol - 1] = maskMat[p][nCol - q - 1]
		downMask[][0,(nCol - 1) / 2 - 1] = maskMat[p][nCol - 1 - q]
	endif
	// Intersection will be 2
	MatrixOp/O MirrorResult = upMask + downMask
	MirrorResult[][] = (MirrorResult[p][q] == 2) ? 1 : 0
	resultW[1] = sum(MirrorResult) // this is the area not the ratio
	// convert result back to sq pixels
	resultW[] *= pxSize ^ 2
	KillWaves/Z maskMat
	return resultW
End

STATIC Function/WAVE GenerateTheMask(ww,hh,xw,yw)
	Variable ww, hh
	Wave xw, yw
	// make an image of boundary (background is 0, boundary is 255)
	ImageBoundaryToMask width=ww, height=hh, xwave=xw, ywave=yw
	WAVE/Z M_ROIMask
	Duplicate/O/FREE M_ROIMask, boundMat
	// make an image with background and boundary filled with 1, interior of ROI is 0
	ImageBoundaryToMask width=ww, height=hh, xwave=xw, ywave=yw, seedX=1, seedY=1
	Duplicate/O/FREE M_ROIMask, outFill, combineMat
	// convert background to grey
	combineMat[][] = (outFill[p][q] == 1 && boundMat[p][q] == 0) ? 127 : 0
	// make the ROI and boundary 1
	combineMat[][] = ((outFill[p][q] == 0 && boundMat[p][q] == 255) || outFill[p][q] == 0 && boundMat[p][q] == 0) ? 1 : combineMat[p][q]
	// now eliminate the background
	combineMat[][] = (combineMat[p][q] == 127) ? 0 : combineMat[p][q]
	
	return combineMat
End

STATIC Function CollectAllMeasurements()
	SetDataFolder root:data:	// relies on earlier load
	DFREF dfr = GetDataFolderDFR()
	String folderName
	Variable numDataFolders = CountObjectsDFR(dfr, 4)
	String wList = ""
	
	Variable i,j
	// assemble a string of semi-colon separated targets in the data folder
	for(i = 0; i < numDataFolders; i += 1)
		folderName = GetIndexedObjNameDFR(dfr, 4, i)
		wList += "root:data:" + folderName + ":thisWave;"
	endfor
	
	// we need to concatenate these waves into root (p_all_*)
	String targetWaveList = "Img_MinAxis;Img_MajAxis;Img_Perimeter;Img_Area;Img_Label;Img_AspectRatio;Img_Circularity;"
	targetWaveList += "Img_ConvexArea;Img_Solidity;Img_Extent;Img_AreaROI;Img_Symmetry;Img_MaxFromCentre;Img_MinFromCentre;"
	Variable nTargets = ItemsInList(targetWaveList)
	String targetName, tList, conName
	
	SetDataFolder root:
	String fullName,modtList
	
	for(i = 0; i < nTargets; i += 1)
		targetName = StringFromList(i,targetWaveList)
		tList = ReplaceString("thisWave",wList,targetName)
		modtList = tList
		// because some waves might not exist
		for(j = 0; j < numDataFolders; j +=1)
			fullName = StringFromList(j, tList)
			Wave testW = $fullName
			if(!WaveExists(testW))
				modtList = RemoveFromList(fullName,modtList)
			endif
		endfor
		// if there were no waves of that type in any of the data folders
		if(ItemsInList(modtList) > 0)
			conName = "all_" + targetName
			Concatenate/O/NP=0 modtList, $conName
		endif
	endfor
End

Function ProcessAllConditions()
	WAVE/Z/T condWave
	WAVE/Z/T all_Img_Label
	
	Variable cond = numpnts(condWave)
	Variable nWaves = dimsize(all_Img_Label,0)
	// all_Img_index will hold index of condition [0]
	// number of contours [1] possibly redundant
	// index of all_Img_Label
	Make/O/N=(nWaves,3) all_Img_Index
	all_Img_Index[][2] = p
	Variable counter
	String condition, dataFolderName
	
	Variable i,j
	
	for(i = 0; i < cond; i += 1)
		counter = 0
		condition = condWave[i]
		for(j = 0; j < nWaves; j += 1)
			dataFolderName = All_Img_Label[j][0]
			if(stringmatch(dataFolderName,condition + "_*") == 1)
				all_Img_Index[j][0] = i
				all_Img_Index[j][1] = counter
				counter += 1
			endif
		endfor
	endfor
	
	// now that all_Img_Index is built, for each condition
	// make a wave to lookup (using index of all_Img_Label) the wave
	String w0Name
	
	for(i = 0; i < cond; i += 1)
		w0Name = "cond" + num2str(i) + "_Img_Index"
		Make/O/N=(nWaves) $w0Name
		Wave w0 = $w0Name
		w0[] = (all_Img_Index[p][0] == i) ? all_Img_Index[p][2] : NaN
		WaveTransform zapnans w0
	endfor
	
	// now we have these waves let's make the versions with the results
	String targetWaveList = "Img_MinAxis;Img_MajAxis;Img_Perimeter;Img_Area;Img_AspectRatio;Img_Circularity;"
	targetWaveList += "Img_ConvexArea;Img_Solidity;Img_Extent;Img_AreaROI;Img_Symmetry;Img_MaxFromCentre;Img_MinFromCentre;"
	Variable nTargets = ItemsInList(targetWaveList)
	String w1Name, tName,w2Name
	Variable nRows
	// which condition has the most cells? We need to know this to make the plots
	WaveStats/Q/RMD=[][1] all_Img_Index
	Variable mostCells = V_Max + 1 // 0-based
	
	for(i = 0; i < cond; i += 1)
		w0Name =  "cond" + num2str(i) + "_Img_Index"
		Wave w0 = $w0Name
		nRows = numpnts(w0)
		for(j = 0; j < nTargets; j += 1)
			tName = "all_" + StringFromList(j,targetWaveList)
			Wave targetW = $tName
			w1Name = "cond" + num2str(i) + "_" + StringFromList(j,targetWaveList)
			Make/O/N=(nRows) $w1Name
			Wave w1 = $w1Name
			w1[] = targetW[w0[p]]
			WaveTransform zapnans w1
			// now store in a 2d matrix
			w2Name = ReplaceString("cond",w1Name,"vb")
			Make/O/N=(mostCells,cond) $w2Name = NaN
			Wave w2 = $w2Name
			w2[0,nRows - 1][i] = w1[p]
		endfor
	endfor
	
	String plotName
	
	for(i = 0; i < nTargets; i += 1)
		plotName =  "p_" + StringFromList(i,targetWaveList)
		KillWindow/Z $plotName
		Display/N=$plotName
		for(j = 0; j < cond; j += 1)
			w0Name = "vb" + num2str(j) + "_" + StringFromList(i,targetWaveList)
			Wave w0 = $w0Name
			BuildBoxOrViolinPlot(w0,plotName,j)
		endfor
		SetAxis/A/N=1/E=1/W=$plotName left
		ModifyGraph/W=$plotName toMode=-1
		ModifyGraph/W=$plotName margin(left)=40
	endfor
	
	// Label y-axes
	Label/W=p_Img_Area left "Area (\u03BCm\S2\M)"
	Label/W=p_Img_Perimeter left "Perimeter (\u03BCm)"
	Label/W=p_Img_minAxis left "Minor axis (\u03BCm)"
	Label/W=p_Img_majAxis left "Major axis (\u03BCm)"
	Label/W=p_Img_AspectRatio left "Aspect Ratio"
	Label/W=p_Img_Circularity left "Circularity"
	Label/W=p_Img_ConvexArea left "Convex area (\u03BCm\S2\M)"
	Label/W=p_Img_Solidity left "Solidity"
	Label/W=p_Img_Extent left "Extent"
	Label/W=p_Img_AreaROI left "Mask area (\u03BCm\S2\M)"
	Label/W=p_Img_Symmetry left "Symmetry"
	Label/W=p_Img_MaxFromCentre left "Max from centre (\u03BCm)"
	Label/W=p_Img_MinFromCentre left "Min from centre (\u03BCm)"
	
	// Look at maj/minor axes on a plot
	plotName = "p_Img_Axes"
	KillWindow/Z $plotName
	Display/N=$plotName
	Variable alphaLevel = DecideOpacity(mostCells)
	WAVE/Z All_Img_MajAxis = root:All_Img_MajAxis
	Variable maxVal = RoundFunction(WaveMax(all_Img_MajAxis),20)
	WAVE/Z colorWave = root:colorWave
	WAVE/Z all_Img_Area = root:all_Img_Area
	Variable biggestCell = WaveMax(all_Img_Area)
	for(i = 0; i < cond; i += 1)
		w0Name = "cond" + num2str(i) + "_Img_MajAxis"
		w1Name = "cond" + num2str(i) + "_Img_MinAxis"
		AppendToGraph/W=$plotName $w1Name vs $w0Name
		ModifyGraph/W=$plotName rgb($w1Name)=(colorWave[i][0],colorWave[i][1],colorWave[i][2],alphaLevel)
		ModifyGraph/W=$plotName rgb($w1Name)=(colorWave[i][0],colorWave[i][1],colorWave[i][2],alphaLevel)
		w2Name = ReplaceString("MinAxis",w1Name,"Area")
		ModifyGraph/W=$plotName zmrkSize($w1Name)={$w2Name,0,biggestCell,1,8}
		ModifyGraph/W=$plotName mrkThick($w1Name)=0
	endfor
	SetAxis/W=$plotName left 0,maxVal
	SetAxis/W=$plotName bottom 0,maxVal
	ModifyGraph/W=$plotName mode=3,marker=19
	ModifyGraph/W=$plotName height={Aspect,1}
	Label/W=$plotName left "Minor axis (\u03BCm)"
	Label/W=$plotName bottom "Major axis (\u03BCm)"
End

Function ImageQuilt()
	// First let's work out how big a quilt to make (they are square)
	Variable qSize = QuiltCalculator()
	// if qSize is 7, we will take 49 outlines and plot them in a 7 x 7 grid
	// however we will use a grid of 9 x 9 to allow a border of one, i.e. 0-8
	// In CellMigration a large matrix was build for each condition with everything offset
	// We'll do everything by offsetting using c_cell* waves
	// Each condition has a condn_IMG_Index wave which gives the index position in All_Img_Label
	WAVE/Z/T condWave = root:condWave
	Variable cond = numpnts(condWave)
	String w0Name, w1Name, w2Name
	WAVE/Z all_Img_Area = root:all_Img_area
	String plotName, cellDF, cellName, plotWave, tName
	WAVE/Z/T all_Img_Label = root:all_Img_Label
	
	Variable i,j
	
	for(i = 0; i < cond; i += 1)
		w0Name = "cond" + num2str(i) + "_Img_Index"
		Wave w0 = $w0Name
		// Sample qSize^2 from each cond index wave
		StatsSample/N=(qSize^2) w0
		WAVE/Z W_Sampled
		w1Name = ReplaceString("Img_Index",w0Name,"Smp_Index")
		Duplicate/O W_Sampled, $w1Name
		Wave w1 = $w1Name
		// Now collect the area values for this sample from all_Img_Area
		w2Name = ReplaceString("Img_Index",w0Name,"Smp_Area")
		Make/O/N=(qSize^2) $w2Name
		Wave w2 = $w2Name
		w2[] = all_Img_Area[w1[p]]
		// Now sort these outlines by their Area
		Sort w2, w2, w1
		// Make the quilt window
		plotName = "quilt_cond" + num2str(i) + "_sample"
		KillWindow/Z $plotName
		Display/N=$plotName
		Variable biggestVal = 0
		// Now add them all to the quilt
		for(j = 0; j < qSize^2; j += 1)
			// w1 holds the index in all_* for a sample of cells ranked by area
			cellDF = all_Img_Label[w1[j]][0]
			cellName = all_Img_Label[w1[j]][1]
			plotWave = "root:data:" + cellDF + ":" + cellName
			Wave w3 = $plotWave
			tName = "outL" + num2str(j)
			AppendToGraph/W=$plotName w3[][1]/TN=$tName vs w3[][0]
			// while we are here let's grab the biggest value so that we can offset
			biggestVal = max(WaveMax(w3),biggestVal)
		endfor
	endfor
	// the next step is to offset everything
	biggestVal = RoundFunction(biggestVal, 10) * 2
	Variable maxVal = (qSize + 1) * biggestVal
	WAVE/Z colorWave = root:colorWave
	for(i = 0; i < cond; i += 1)
		plotName = "quilt_cond" + num2str(i) + "_sample"
		for(j = 0; j < qSize^2; j += 1)
			tName = "outL" + num2str(j)
			ModifyGraph/W=$plotName offset($tName)={(mod(j,qSize) + 1) * biggestVal,(floor(j / qSize) + 1) * biggestVal}
		endfor
		ModifyGraph/W=$plotName rgb=(colorWave[i][0],colorWave[i][1],colorWave[i][2])
		SetAxis/W=$plotName left maxVal,0
		SetAxis/W=$plotName bottom 0, maxVal
		ModifyGraph/W=$plotName height={Aspect,1}
		ModifyGraph/W=$plotName noLabel=2,axThick=0,standoff=0
		ModifyGraph/W=$plotName margin=7
	endfor
End

STATIC Function QuiltCalculator()
	WAVE/Z/T labelWave = root:labelWave
	Variable cond = numpnts(labelWave) // use labelWave rather than condWave for printing
	Variable nCells= 0, leastCells = 10000
	String w0Name
	
	Variable i
	
	for(i = 0; i < cond; i += 1)
		w0Name = "cond" + num2str(i) + "_Img_Index"
		Wave w0 = $w0Name
		nCells = numpnts(w0)
		Print labelWave[i], "has", nCells, "outlines"
		leastCells = min(leastCells,nCells)
	endfor
	Variable qSize = floor(sqrt(leastCells))
	Print "Smallest number of outlines is", leastCells
	Print "Quilt of", qSize, "x", qSize
	return qSize
End

// This function will make a "multicolumn" boxplot or violinplot (Igor >8 only) 
///	@param	matA	matrix of points to be appended
///	@param	plotName	string to tell igor which graph window to work on
///	@param	ii	variable to indicate which condition (for coloring)
STATIC Function BuildBoxOrViolinPlot(matA,plotName,ii)
	WAVE matA
	String plotName
	Variable ii
	
	String wName = NameOfWave(matA)
	Wave/T/Z labelWave = root:labelWave
	Wave/Z colorWave = root:colorWave
	//  This works because all matrices passed to this function have the same dimensions
	Variable nCells = DimSize(matA,0)
	if(nCells < 100)
		AppendBoxPlot/W=$plotName matA vs labelWave
		ModifyBoxPlot/W=$plotName trace=$wName,markers={19,-1,19},markerSizes={2,2,2}
		ModifyBoxPlot/W=$plotName trace=$wName,whiskerMethod=4
	else
		AppendViolinPlot/W=$plotName matA vs labelWave
		ModifyViolinPlot/W=$plotName trace=$wName,ShowMean,MeanMarker=19,CloseOutline
		ModifyViolinPlot/W=$plotName trace=$wName,DataMarker=19
	endif
	Variable alphaLevel = DecideOpacity(nCells)
	ModifyGraph/W=$plotName rgb($wName)=(colorWave[ii][0],colorWave[ii][1],colorWave[ii][2],alphaLevel)
End

////////////////////////////////////////////////////////////////////////
// Utility functions
////////////////////////////////////////////////////////////////////////
Function CleanSlate()
	String fullList = WinList("*", ";","WIN:71")
	Variable allItems = ItemsInList(fullList)
	String name
	Variable i
 
	for(i = 0; i < allItems; i += 1)
		name = StringFromList(i, fullList)
		KillWindow/Z $name		
	endfor
	
	KillDataFolder/Z root:data:
		
	// Kill waves in root
	KillWaves/A/Z
	// Look for data folders and kill them
	DFREF dfr = GetDataFolderDFR()
	allItems = CountObjectsDFR(dfr, 4)
	for(i = 0; i < allItems; i += 1)
		name = GetIndexedObjNameDFR(dfr, 4, i)
		KillDataFolder $name		
	endfor
End

STATIC Function KillTheseWaves(wList)
	String wList
	Variable nWaves = ItemsInList(wList)
	String wName
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Wave w0 = $wName
		KillWaves/Z w0
	endfor
End

STATIC Function TheDiviner(fList)
	String fList
	Wave/T fileNameFWave = ListToTextWave(fList,";")
	MoveWave fileNameFWave, root:fileNameWave // save a copy in root
	
	Variable nFiles = numpnts(fileNameFWave)
	Make/O/N=(nFiles)/T/FREE shortNameWave
	// expression here. In future this could be determined by looking for separators first
	String expr="([[:alnum:]]+)\\w([[:digit:]]+)\\w([[:digit:]]+)"
	String cond, expt, cell
	
	Variable i
	
	for(i = 0; i < nFiles; i += 1)
		// put all the conditions in to shortNameWave
		SplitString/E=(expr) fileNameFWave[i], cond, expt, cell
		shortNameWave[i] = cond
	endfor
	FindDuplicates/RT=condWave shortNameWave
	WAVE/Z/T condWave
	if(numpnts(condWave) == 1)
		return 0
	elseif(numpnts(condWave) == 0)
		return -1 // error no conditions found
	endif
	// now send condWave to the dialog and if successful return 0
	if(Diviner2User(condWave) != 0)
		return 0
	endif
End

// Colours are taken from Paul Tol SRON stylesheet
// Colours updated. Brighter palette for up to 6 colours, then palette of 12 for > 6
// Define colours
StrConstant SRON_1 = "0x4477aa;"
StrConstant SRON_2 = "0x4477aa;0xee6677;"
StrConstant SRON_3 = "0x4477aa;0xccbb44;0xee6677;"
StrConstant SRON_4 = "0x4477aa;0x228833;0xccbb44;0xee6677;"
StrConstant SRON_5 = "0x4477aa;0x66ccee;0x228833;0xccbb44;0xee6677;"
StrConstant SRON_6 = "0x4477aa;0x66ccee;0x228833;0xccbb44;0xee6677;0xaa3377;"
StrConstant SRON_7 = "0x332288;0x88ccee;0x44aa99;0x117733;0xddcc77;0xcc6677;0xaa4499;"
StrConstant SRON_8 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0xcc6677;0xaa4499;"
StrConstant SRON_9 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_10 = "0x332288;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_11 = "0x332288;0x6699cc;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0x882255;0xaa4499;"
StrConstant SRON_12 = "0x332288;0x6699cc;0x88ccee;0x44aa99;0x117733;0x999933;0xddcc77;0x661100;0xcc6677;0xaa4466;0x882255;0xaa4499;"

/// @param hex		variable in hexadecimal
Function hexcolor_red(hex)
	Variable hex
	return byte_value(hex, 2) * 2^8
End

/// @param hex		variable in hexadecimal
Function hexcolor_green(hex)
	Variable hex
	return byte_value(hex, 1) * 2^8
End

/// @param hex		variable in hexadecimal
Function hexcolor_blue(hex)
	Variable hex
	return byte_value(hex, 0) * 2^8
End

/// @param data	variable in hexadecimal
/// @param byte	variable to determine R, G or B value
STATIC Function byte_value(data, byte)
	Variable data
	Variable byte
	return (data & (0xFF * (2^(8*byte)))) / (2^(8*byte))
End

/// @param	cond	variable for number of conditions
Function MakeColorWave(cond)
	Variable cond
	
	// Pick colours from SRON palettes
	String pal
	if(cond == 1)
		pal = SRON_1
	elseif(cond == 2)
		pal = SRON_2
	elseif(cond == 3)
		pal = SRON_3
	elseif(cond == 4)
		pal = SRON_4
	elseif(cond == 5)
		pal = SRON_5
	elseif(cond == 6)
		pal = SRON_6
	elseif(cond == 7)
		pal = SRON_7
	elseif(cond == 8)
		pal = SRON_8
	elseif(cond == 9)
		pal = SRON_9
	elseif(cond == 10)
		pal = SRON_10
	elseif(cond == 11)
		pal = SRON_11
	else
		pal = SRON_12
	endif
	
	Variable color
	Make/O/N=(cond,3) root:colorWave
	WAVE colorWave = root:colorWave
	Variable i
	
	for(i = 0; i < cond; i += 1)
		// specify colours
		color = str2num(StringFromList(mod(i, 12),pal))
		colorwave[i][0] = hexcolor_red(color)
		colorwave[i][1] = hexcolor_green(color)
		colorwave[i][2] = hexcolor_blue(color)
	endfor
End

STATIC Function MakeTheLayouts(prefix,nRow,nCol,[iter, filtVar, rev, saveIt, orient])
	String prefix
	Variable nRow, nCol
	Variable iter	// this is if we are doing multiple iterations of the same layout
	Variable filtVar // this is the object we want to filter for
	Variable rev // optional - reverse plot order
	Variable saveIt
	Variable orient //optional 1 = landscape, 0 or default is portrait
	if(ParamIsDefault(filtVar) == 0)
		String filtStr = prefix + "_*_" + num2str(filtVar) + "_*"	// this is if we want to filter for this string from the prefix
	endif
	
	String layoutName = "all"+prefix+"Layout"
	DoWindow/K $layoutName
	NewLayout/N=$layoutName
	String allList = WinList(prefix+"*",";","WIN:1") // edited this line from previous version
	String modList = allList
	Variable nWindows = ItemsInList(allList)
	String plotName
	
	Variable i
	
	if(ParamIsDefault(filtVar) == 0)
		modList = "" // reinitialise
		for(i = 0; i < nWindows; i += 1)
			plotName = StringFromList(i,allList)
			if(stringmatch(plotName,filtStr) == 1)
				modList += plotName + ";"
			endif
		endfor
	endif
	nWindows = ItemsInList(modList)
	Variable PlotsPerPage = nRow * nCol
	String exString = "Tile/A=(" + num2str(ceil(PlotsPerPage/nCol)) + ","+num2str(nCol)+")"
	
	Variable pgNum=1
	
	for(i = 0; i < nWindows; i += 1)
		if(ParamIsDefault(rev) == 0)
			if(rev == 1)
				plotName = StringFromList(nWindows - 1 - i,modList)
			else
				plotName = StringFromList(i,modList)
			endif
		else
			plotName = StringFromList(i,modList)
		endif
		AppendLayoutObject/W=$layoutName/PAGE=(pgnum) graph $plotName
		if(mod((i + 1),PlotsPerPage) == 0 || i == (nWindows -1)) // if page is full or it's the last plot
			if(ParamIsDefault(orient) == 0)
				if(orient == 1)
					LayoutPageAction size(-1)=(842,595), margins(-1)=(18, 18, 18, 18)
				endif
			else
				// default is for portrait
				LayoutPageAction/W=$layoutName size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
			endif
			ModifyLayout/W=$layoutName units=0
			ModifyLayout/W=$layoutName frame=0,trans=1
			Execute /Q exString
			if (i != nWindows -1)
				LayoutPageAction/W=$layoutName appendpage
				pgNum += 1
				LayoutPageAction/W=$layoutName page=(pgNum)
			endif
		endif
	endfor
	
	String fileName
	// if anthing is passed here we save an iteration, otherwise usual name
	if(!ParamIsDefault(iter))
		fileName = layoutName + num2str(iter) + ".pdf"
	else
		fileName = layoutName + ".pdf"
	endif
	// if anthing is passed here we save the filtered version
	if(ParamIsDefault(filtVar) == 0)
		fileName = ReplaceString(".pdf",fileName, "_" + num2str(filtVar) + ".pdf")
	endif
	if(ParamIsDefault(saveIt) == 0)
		if(saveIt == 1)
			SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
		endif
	else
		// default is to save
		SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
	endif
End

STATIC Function TidyAndSave(prefix)
	String prefix
	String layoutName = "all"+prefix+"Layout"
	// go to first page
	LayoutPageAction/W=$layoutName page=(1)
	// build the key
	WAVE/Z/T labelWave = root:labelWave
	WAVE/Z colorWave = root:colorWave
	Variable cond = numpnts(labelWave)
	String boxString = ""
	
	Variable i
	
	for(i = 0; i < cond; i += 1)
		// add text colour for condition
		boxString += "\\K(" + num2str(colorWave[i][0]) + "," + num2str(colorWave[i][1]) + "," + num2str(colorWave[i][2])
		boxString += ")" + labelWave[i]
		if (i < cond - 1)
			boxString += "\r"
		endif
	endfor
	TextBox/W=$layoutName/C/N=text0/F=0/A=RB/X=5.00/Y=5.00 boxString
	String fileName = layoutName + ".pdf"
	SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
End

STATIC Function DecideOpacity(nTrace)
	Variable nTrace
	Variable alpha
	if(nTrace < 10)
		alpha = 1
	elseif(nTrace < 50)
		alpha = 0.5
	elseif(nTrace < 100)
		alpha = 0.3
	else
		alpha = 0.2
	endif
	alpha = round(65535 * alpha)
	return alpha
End

// for axis scaling
///	@param	value	this is the input value that requires rounding up
///	@param	roundto	round to the nearest...
STATIC Function RoundFunction(value,roundTo)
	Variable value, roundTo
	
	value /= roundTo
	Variable newVal = ceil(value)
	newVal *= roundTo
	return newVal
End

////////////////////////////////////////////////////////////////////////
// Panel functions
////////////////////////////////////////////////////////////////////////
STATIC Function Diviner2User(condWave)
	WAVE/T condWave
	DoWindow/F SelectPanel
	if (V_Flag != 0)
		return 0
	endif
	Make/O/N=(numpnts(condWave)) root:testSel
	NewPanel/N=SelectPanel/K=1/W=(350,125,650,325) as "Put conditions in order"
	ListBox list0,pos={1,2},size={298,157},proc=ListBoxProc_DragNDropLB
	ListBox list0,listWave=root:condWave,selWave=root:testSel,mode= 1,selRow= 1
	Button DoIt,pos={100,177},size={100,20},proc=ButtonProc,title="Do It"
End
 
Function ListBoxProc_DragNDropLB(lba) : ListBoxControl
	STRUCT WMListboxAction &lba
 
	Variable row = lba.row
	Variable col = lba.col
	WAVE/T/Z listWave = lba.listWave
	WAVE/Z selWave = lba.selWave
 
	switch( lba.eventCode )
		case -1: // control being killed
			break
		case 1: // mouse down
			Variable/G V_MouseDownRow = row
			break
		case 2: // mouse up
			//quit if V_MouseDownRow isn't a string or numeric variable
			//otherwise create a pointer to it.
			if( exists( "V_MouseDownRow" ) == 2 )
				NVAR V_MouseDownRow
			else
				break
			endif
			if(row != V_MouseDownRow)						// dragged?
				String item = listWave[V_MouseDownRow]
				DeletePoints V_MouseDownRow, 1, listWave	// do swap
				InsertPoints row, 1, listWave
				listWave[row] = item
			endif
			KillVariables V_MouseDownRow	// cleanup variable
			break
		case 3: // double click
			break
		case 4: // cell selection
		case 5: // cell selection plus shift key
			break
		case 6: // begin edit
			break
		case 7: // finish edit
			break
	endswitch
 	//print "code 1: ", lba.eventcode
	return 1
End

// Define button
Function ButtonProc(ctrlName) : ButtonControl
	String ctrlName
		strswitch(ctrlName) 
			case "DoIt" :
				// run the rest of the program
				KillWindow/Z SelectPanel
				OptionalAliases()
				//break
		endswitch
End

// This is the Alias check to make a label wave
Function OptionalAliases()
	Wave/T/Z condWave = root:condWave
	Variable cond = numpnts(condWave)
	MakeColorWave(cond)
	Wave/Z colorWave = root:colorWave
	// duplicate the condWave to make a tentative labelWave
	Duplicate/O condWave, root:labelWave
	Wave/T/Z labelWave = root:labelWave

	DoWindow/K AliasCheck
	NewPanel/N=AliasCheck/K=1/W=(40, 40, 460, 150+30*cond)
	// labelling of columns
	DrawText/W=AliasCheck 10,30,"Name"
	DrawText/W=AliasCheck 160,30,"Alias (a nice name for the plot labels)"
	DrawText/W=AliasCheck 10,100+30*cond,"Cell Shape Analysis"
	// do it button
	Button DoIt,pos={280,70+30*cond},size={100,20},proc=DoItButtonProc,title="Do It"
	// insert rows
	String buttonName1a,buttonName1b,buttonName2a,buttonName2b,boxName0,boxName1,boxName2
	Variable i
	
	for(i = 0; i < cond; i += 1)
		boxName0 = "box0_" + num2str(i)
		boxName1 = "box1_" + num2str(i)
		// row label
		DrawText/W=AliasCheck 10,68+i*30,num2str(i+1)
		// condition label
		SetVariable $boxName0,pos={30,53+i*30},size={100,14},value= condWave[i], title=" "
		// file or dir box
		SetVariable $boxName1,pos={160,53+i*30},size={220,14},value= labelWave[i], title=" "
		SetDrawEnv fillfgc=(colorWave[i][0],colorWave[i][1],colorWave[i][2])
		DrawOval/W=AliasCheck 130,50+i*30,148,68+i*30
	endfor
End

// define buttons
Function DoItButtonProc(ctrlName) : ButtonControl
	String ctrlName
 	
 	WAVE/T/Z labelWave = root:labelWave
	Variable okvar = 0
	
	strswitch(ctrlName)	
		case "DoIt" :
			// check CondWave
			okvar = WaveChecker(labelWave)
			if (okvar == -1)
				Print "Error: Not all conditions have a name."
				break
			endif
			okvar = NameChecker(labelWave)
			if (okvar == -1)
				Print "Error: Two conditions have the same name."
				break
			else
				KillWindow/Z AliasCheck
				// now execute main program
				KillWaves/Z testSel
				TheLoader()
				return 0
			endif
	endswitch	
End

STATIC function WaveChecker(TextWaveToCheck)
	Wave/T TextWaveToCheck
	Variable nRows = numpnts(TextWaveToCheck)
	Variable len
	
	Variable i
	
	for(i = 0; i < nRows; i += 1)
		len = strlen(TextWaveToCheck[i])
		if(len == 0)
			return -1
		elseif(numtype(len) == 2)
			return -1
		endif
	endfor
	return 1
End

STATIC function NameChecker(TextWaveToCheck)
	Wave/T TextWaveToCheck
	Variable nRows = numpnts(TextWaveToCheck)
	Variable len
	
	Variable i,j
	
	for(i = 0; i < nRows; i += 1)
		for(j = 0; j < nRows; j += 1)
			if(j > i)
				if(cmpstr(TextWaveToCheck[i], TextWaveToCheck[j], 0) == 0)
					return -1
				endif
			endif
		endfor
	endfor
	return 1
End