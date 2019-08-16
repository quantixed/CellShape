#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.

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

Function TheLoader()
	LoadIMODModels()
	ProcessAllModels()
//	MakeTheLayouts("q_",6,8)
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
		CentreAndPlot()
		TakeMeasurements()
		SetDataFolder root:data:
	endfor
	SetDataFolder root:
End

STATIC Function CentreAndPlot()
	String wList = WaveList("cell*",";","")
	Variable nWaves = ItemsInList(wList)
	String wName, w0Name
	String folderName = GetDataFolder(0)
	
	Variable i
	
	for(i = 0; i < nWaves; i += 1)
		wName = StringFromList(i,wList)
		w0Name = "root:data:" + folderName + ":" + wName
		Wave/Z w0 = $w0Name
		MakeCellggPlot(w0)
	endfor
End

STATIC Function MakeCellggPlot(w0, [plusMinus])
	Wave w0
	Variable plusMinus
	// the wave generated here is centred and/or rotated. It's called c_*
	Wave w2 = Centralise2DWave(w0,midpoint = 1) // centring/rotation option here
	String plotName = "q_" + GetWavesDataFolder(w0,0) + "_" + NameOfWave(w0)
	KillWindow/Z $plotname
	Display/N=$plotName/HIDE=1 w2[][1]/TN=Outline0 vs w2[][0]
	ModifyGraph/W=$plotName rgb(Outline0)=(68*257,170*257,153*257)
	ModifyGraph/W=$plotName width={Aspect,1}
	// now work out how to display the image
	if(ParamIsDefault(plusMinus) == 1)
		Variable last2 = dimSize(w2,0)
		Variable xMinVal = WaveMin(w2,0,last2 - 1)
		Variable xMaxVal = WaveMax(w2,0,last2 - 1)
		xMaxVal = Max(xMaxVal,abs(xMinVal))
		Variable yMinVal = WaveMin(w2,last2,last2*2 - 1)
		Variable yMaxVal = WaveMax(w2,last2,last2*2 - 1)
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
	String wList = WaveList("c_cell_*",";","")
	Variable nWaves = ItemsInList(wList)
	Make/O/N=(nWaves) Img_MinAxis, Img_MajAxis, Img_Perimeter, Img_Area
	String currentDF = GetDataFolder(0)
	String wName,tName
	
	Variable i,j
	
	for (i = 0; i < nWaves; i += 1)
		wName = StringFromList(i, wList)
		Wave w0 = $wName
		Img_MinAxis[i] = VesicleAxisLength(w0,1)
		Img_MajAxis[i] = VesicleAxisLength(w0,0)
		Img_Perimeter[i] = FindLengthOfXYCoords(w0)
		MatrixOp/O/FREE w0c0 = col(w0,0)
		MatrixOp/O/FREE w0c1 = col(w0,1)
		Img_Area[i] = PolygonArea(w0c0,w0c1)
	endfor
	if(numpnts(Img_Area) > 0)
		MatrixOp/O/NTHR=0 Img_AspectRatio = Img_MinAxis / Img_MajAxis
		MatrixOp/O/NTHR=0 Img_Circularity = (4 * pi * Img_Area) / (Img_Perimeter * Img_Perimeter)
//	elseif(numpnts(Img_Area) == 1)
//		MatrixOp/O/NTHR=0 Img_AspectRatio = Img_MinAxis / Img_MajAxis
//		MatrixOp/O/FREE/NTHR=0 tempMat = Img_Perimeter * Img_Perimeter
//		MatrixOp/O/NTHR=0 Img_Circularity = 4 * pi * tempMat
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
	Make/O/N=(cond,3) root:colorwave
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

STATIC Function MakeTheLayouts(prefix,nRow,nCol,[iter, filtVar])
	String prefix
	Variable nRow, nCol
	Variable iter	// this is if we are doing multiple iterations of the same layout
	Variable filtVar // this is the object we want to filter for
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
		plotName = StringFromList(i,modList)
		AppendLayoutObject/W=$layoutName/PAGE=(pgnum) graph $plotName
		if(mod((i + 1),PlotsPerPage) == 0 || i == (nWindows -1)) // if page is full or it's the last plot
			LayoutPageAction/W=$layoutName size(-1)=(595, 842), margins(-1)=(18, 18, 18, 18)
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
	if(!ParamIsDefault(iter))
		fileName = layoutName + num2str(iter) + ".pdf"
	else
		fileName = layoutName + ".pdf"
	endif
	if(ParamIsDefault(filtVar) == 0)
		fileName = ReplaceString(".pdf",fileName, "_" + num2str(filtVar) + ".pdf")
	endif
	String folderStr
	SavePICT/O/WIN=$layoutName/PGR=(1,-1)/E=-2/W=(0,0,0,0) as fileName
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
				LoadIMODModels()
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