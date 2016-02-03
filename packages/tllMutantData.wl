(* ::Package:: *)

BeginPackage["tllMutantData`", { "utilities`"}]

dataWT::usage="Janssens et al WT data"
dataMutant::usage="Janssens et al mutant data"
probeNames::usage="Janssens et al probenames"

showPred::usage="Plot Janssens et al predictions"


Begin["`Private`"] (* Begin Private Context *) 



parseJanssens[filename_]:=Module[{gene,timepoint,text,data},
	{{gene,timepoint}}=StringCases[filename,("dmwt_"|"dmtllg_"|"dmtg_")~~geneStr:Shortest[__]~~"_t"~~timepointStr:DigitCharacter~~(".100"|".txt"):>{geneStr,timepointStr},1];
	text=Import[filename,"Text"];
	data=ImportString[StringReplace[text,StartOfLine~~(""|("#"~~Shortest[___]))~~"\n"->""] ,"Table"];
	<|timepoint->gene->data[[All,3]]|>
]

datasetTimePoint[ds_Association,timepoint_String,evedelay_Integer]:=Module[{firstDs,delayedEve},
	firstDs=Dataset[(Transpose[Thread/@ds[[timepoint]]] )];
	delayedEve="eve" /.ds[[ToString[Interpreter["Integer"][timepoint]+evedelay]]];
	MapThread[Append,{Normal[firstDs[All,KeyDrop["eve"]]],Thread["eve"->delayedEve]}]//Dataset
]

normaliseJanssens[data_]:=Module[{byGene,normByGene},
	byGene=Transpose[Values@data];
	normByGene=Thread[#[[All,1]]->Rescale[#[[All,2]]]*1.1]&/@byGene;
	Association[Thread[Keys@data->Transpose[normByGene]]]
]


showPred[predfunction3_,predfunction4_,inputData_Association,timepoint_Integer,delay_Integer]:=Module[{ds,newdata},
	ds=datasetTimePoint[inputData,ToString@timepoint,delay] ;
	newdata=ds[All,Function[Evaluate[<|#,"pred3"->predfunction3,"pred4"->predfunction4|>]]];
	ListPlot[newdata[28;;92,{"eve","hb","kni","tll","pred3","pred4"}]//Normal //Values //Transpose,Joined->True,PlotRange->All,ImageSize->250,PlotStyle->{{Gray,Thin},{Blue,Thin},{Orange,Thin},{Green,Thin},{Pink,Thick},{Thick,RGBColor[0.22,0.57,0.6]}},Ticks->None,AspectRatio->1/3]
]


(* Import data *)
SetDirectory[FileNameJoin[{ParentDirectory[DirectoryName[$InputFileName]],"externalData/JanssensJaeger"}]]
filenamesWT=FileNames[("dmwt_")~~__~~".txt"];
filenamesMutant=FileNames[("dmtg_")~~__~~".100"];
dataWT=normaliseJanssens[Merge[parseJanssens /@ filenamesWT ,Join]];
dataMutant=Append[#,"tll"->ConstantArray[0,100]]& /@normaliseJanssens[Merge[parseJanssens /@ filenamesMutant ,Join]];
probeNames=Normal[datasetTimePoint[dataWT,"4",0][1,Keys]];
ResetDirectory[]


End[] (* End Private Context *)   
EndPackage[]
