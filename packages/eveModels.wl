(* ::Package:: *)

(* Mathematica Package *)

BeginPackage["eveModels`", { "BDTNP`","utilities`"}]

embryoDS::usage = "The BDTNP embryo in Dataset format"

(* AveragedEmbryo *)
AveragedEmbryo::usage = "AveragedEmbryo[\!\(\*
StyleBox[\"dataset\",\nFontSlant->\"Italic\"]\)] returns an averaged embryo."
Options[AveragedEmbryo] = {infoCohort -> 3, mRNACohorts -> {1, 2}, proteinCohorts -> {1, 2},outputGene->"eve"} 

(* MutantEmbryo *)
MutantEmbryo::usage = "MutantEmbryo[Dataset,gene,expr] returns an embryo with mutant data."

(* nonlinearFit *)
nonlinearFit::usage = "nonlinearFit[Dataset,regs] fits a nonlinear model to the Dataset."
Options[nonlinearFit]={outputGene->"eve",SQregulators -> {},scale->False}

(* logisticFit *)
logisticFit::usage = "logisticFit[Dataset,regs] fits a logistic regression model to the Dataset."
Options[logisticFit]={outputGene->"eve",SQregulators -> {},outputThreshold->0.2}  

(* FlyModel *)
FlyModel::usage=="FlyModel[stripes, regulators]"
FlyModel::modelMethod="`1` is not a defined model method."
FlyModel::trainingSelection="`1` is not a defined training data selection method."
Options[FlyModel] = Union[Options[logisticFit],Options[nonlinearFit],Options[AveragedEmbryo], {modelMethod -> "logisticRegression", trainingSelection->"standard",trainingFilter -> True, trainingExcludeStripes->{},trainingModifyData->{}}];

(* PlotEmbryo *)
PlotEmbryo::usage = "PlotEmbryo[Dataset,gene]"
Options[PlotEmbryo]=Join[Options[ListPlot],{colour->Red,maxcolor->Green,stripes->{},dimensions->2,pointSize2D->0.012,maxvalue->1,view->"side"}]

(* PlotEmbryoTll *)
PlotEmbryoTll::usage = "PlotEmbryo[Dataset,gene]"
Options[PlotEmbryoTll]=Options[PlotEmbryo]

(* PlotModel *)
PlotModel::usage = "PlotModel[Dataset,model,stripes] plots a fitted model." (* There are other usage(s) *)
PlotModel::plotType="`1` is not a valid plotType."
Options[PlotModel]=Complement[Join[Options[PlotEmbryo],{plotType->"embryo"}],{"maxvalue"}]

(* PlotMutantSeries *)
PlotMutantSeries::usage=="PlotMutantSeries[flyResults, probe]"
Options[PlotMutantSeries] = Join[Options[PlotModel],{SeriesRange->Range[1,0, -0.25]}]

(* PlotTALE *)
PlotTALE::usage=="PlotTALE[flyResults, range]"
Options[PlotTALE] = Join[Options[PlotModel],{CompositeImageSize->1500,CompositeOrientation->"horizontal"}]


Begin["`Private`"] (* Begin Private Context *) 


embryoDS=Dataset[BDTNP`embryoAssoc]
nucleiStripesDS=Dataset[BDTNP`nucleiStripesAssoc]


stripesAndNeighbouringNucIDs[stripes_List]:=Module[{stripesNuclei},
	stripesNuclei=nucleiStripesDS[Select[MemberQ[stripes, #stripe] &], "nucleus_id"] // Normal;
	VertexList@NeighborhoodGraph[nucleiGraph,stripesNuclei]
]

standardTrainingNucIDs[stripes_List] := Module[{otherStripes,exclNuclei},
	(* This excludes the other stripes and the neighbouring nuclei *)
	otherStripes=Complement[Range@7,stripes];
	exclNuclei=stripesAndNeighbouringNucIDs[otherStripes];
	Complement[Normal@nucleiStripesDS[All,"nucleus_id"],exclNuclei]
]


AveragedEmbryo[data_Dataset, OptionsPattern[]] :=
    Module[ {infocolumns, exclcolumns,allcolumns, proteincolumns,mRNAcolumns,infoData, proteinData, mRNAData},
      		allcolumns = data[1, Keys] // Normal; (* assumes keys are the same throughout *)
      		infocolumns = {"nucleus_id", "x", "y", "z","stripe"};
      		exclcolumns = {"cohort","Nx", "Ny", "Nz"};
      		proteincolumns = Cases[colname_String /;StringMatchQ[colname, (__ ~~ "P")| (OptionValue[outputGene])]][allcolumns];
      		mRNAcolumns=Complement[allcolumns,infocolumns,proteincolumns,exclcolumns];
      			
      		infoData = data[Select[#cohort==OptionValue[infoCohort] &],infocolumns];
			  proteinData = data[GroupBy[KeyTake["nucleus_id"]], Select[MemberQ[OptionValue[proteinCohorts], #cohort ] &],proteincolumns][All, Mean][Normal][All, Apply[Association]];
			  mRNAData = data[GroupBy[KeyTake["nucleus_id"]], Select[MemberQ[OptionValue[mRNACohorts], #cohort ] &],mRNAcolumns][All, Mean][Normal][All, Apply[Association]];
  			JoinAcross[JoinAcross[Normal@infoData,Normal@proteinData,"nucleus_id"],Normal@mRNAData,"nucleus_id"] //Dataset
    ]

MutantEmbryo[ds_Dataset,gene_,fun_] :=
    	ds[All,<|#,gene-> fun|> &]

MutantEmbryo[ds_Dataset,gene_,fun_:Function[0]] :=
    	ds[All,<|#,gene-> fun|> &]

ModifyEmbryo[data_Dataset,List[gene_,expr_]] := MutantEmbryo[data,gene,expr]


nonlinearFit[ds_Dataset,regs_List,opts: OptionsPattern[]] :=
    Module[ {allregs,data,parms,betaparms,vars,sqvars,fitformula},
        allregs = DeleteDuplicates[Join[regs,OptionValue[SQregulators]]];
        vars = Symbol /@allregs;
        sqvars = Symbol[#]^2& /@ OptionValue[SQregulators];
        betaparms = Map[Subscript[\[Beta],#]&,Range[0,Length[allregs]+Length[OptionValue[SQregulators]]]];
		parms=If[OptionValue[scale],Join[betaparms,{\[Alpha]}],betaparms];
        data = ds[All,Append[allregs,OptionValue[outputGene]]] // Normal // Values;
        fitformula = If[OptionValue[scale],
			\[Alpha]/(1+Exp[-(Subscript[\[Beta], 0]+Rest[betaparms] . Union[vars,sqvars])]),
        	1/(1+Exp[-(Subscript[\[Beta], 0]+Rest[betaparms] . Union[vars,sqvars])])];
        NonlinearModelFit[data,fitformula,parms,vars,MaxIterations->5500, Method-> {"NMinimize",Method->"NelderMead"}]
    ]

logisticFit[ds_Dataset,regs_List,opts: OptionsPattern[]] :=
    Module[ {allregs,data,parms,vars,sqvars},
        allregs = DeleteDuplicates[Join[regs,OptionValue[SQregulators]]];
        vars = Symbol /@allregs;
        sqvars = Symbol[#]^2& /@ OptionValue[SQregulators];
        parms = Map[Subscript[\[Beta],#]&,Range[0,Length[allregs]+Length[OptionValue[SQregulators]]]];
		data=ds[All,Append[allregs,OptionValue[outputGene]]] // Normal // Values;
        data[[All,-1]] = If[ # > OptionValue[outputThreshold],
                             1,
                             0
                         ]&/@data[[All,-1]];
        LogitModelFit[data,Join[vars,sqvars],vars]
    ]



trainDataResetRest[ds_Dataset,stripes_List] :=
    Module[ {},
		ds[All,<|#,"eve"->If[!MemberQ[stripes,#stripe],0,#eve]|>&]
    ]


FlyModel[stripes_, regulators_, opts : OptionsPattern[]] := Module[{averagedEmbryo,modifiedEmbryo,baseEmbryo, baseNucIDs, exclNucIDs,filteredEmbryo,trainData, model,flyResults},
  averagedEmbryo=AveragedEmbryo[embryoDS,FilterRules[{opts}, Options[AveragedEmbryo]]];
  modifiedEmbryo=Fold[ModifyEmbryo,averagedEmbryo,OptionValue[trainingModifyData]];
  flyResults["embryo"] = modifiedEmbryo;
  baseEmbryo=modifiedEmbryo[Select[OptionValue[trainingFilter] &]];
  trainData=Switch[OptionValue[trainingSelection],
  	(* standard and resetRest are written for eve as the output gene only.  *)
  	 "eLife",(baseNucIDs=standardTrainingNucIDs[stripes];baseEmbryo[Select[MemberQ[baseNucIDs, #["nucleus_id"]]&]]),
  	 "resetRest",(exclNucIDs = stripesAndNeighbouringNucIDs[OptionValue[trainingExcludeStripes]];
				   filteredEmbryo=baseEmbryo[Select[!MemberQ[exclNucIDs, #["nucleus_id"]] &]];
  				 trainDataResetRest[filteredEmbryo, stripes]),
	   "all",baseEmbryo,
  	_,Return[Message[FlyModel::trainingSelection,OptionValue[trainingSelection]]]];
  model = Switch[OptionValue[modelMethod],
  	"nonlinear",nonlinearFit[trainData,regulators, FilterRules[{opts},Options[nonlinearFit]]],
  	"logisticRegression",logisticFit[trainData, regulators, FilterRules[{opts},Options[logisticFit]]],
  	_,Return[Message[FlyModel::modelMethod,OptionValue[modelMethod]]]]; 
  flyResults["model"] = model;
  flyResults["trainData"] = trainData;
  flyResults
]


PlotEmbryo[ds_Dataset,value_String,opts:OptionsPattern[]] :=
    Module[ {plotData,profileData, profileDataSplit,colournucids,offset,maxpredvalue},
        plotData = If[ OptionValue[dimensions]==3,
                       ds[All,{"x","y","z",value}],
						If[OptionValue[view]=="side",
                       ds[All,{"x","z",value}],
					   ds[Select[#z<0&],{"x","y",value}]]
                   ] //Normal //Values;
        offset = If[ OptionValue[dimensions]==3,
                     Min[plotData[[All,3]]]-140, (* was -100 *)
                     Min[plotData[[All,2]]]-140
                 ];
		maxpredvalue=OptionValue[maxvalue];
        If[ Length[plotData]>6078,
            Print["More than one cohort."],
            If[ Length[OptionValue[stripes]]>0,
                (
                colournucids = nucleiStripesDS[Select[MemberQ[OptionValue[stripes], #stripe] &], "nucleus_id"] //  Normal;
				profileData = ds[Select[Abs[#z]<20&],{#x,#[value]*120+offset,#["nucleus_id"]}&];
				profileDataSplit = GatherBy[profileData,MemberQ[colournucids,#[[3]]]&];
                Show[
                  Graphics[({Blend[{{0,LightPink},{maxpredvalue,OptionValue[colour]},{2*maxpredvalue,OptionValue[maxcolor]}} ,#[[3]]],PointSize[OptionValue[pointSize2D]],Point[#[[1;;2]]]}) & /@ plotData],
                  ListPlot[profileDataSplit[[1,All,1;;2]],FilterRules[{opts},Options[ListPlot]],Axes->False,PlotStyle->Directive[Black,PointSize[Small]],PlotRange->Full],
                  ListPlot[profileDataSplit[[2,All,1;;2]],FilterRules[{opts},Options[ListPlot]],Axes->False,PlotStyle->Directive[Darker[Green],PointSize[Small]],PlotRange->Full]
                ]	
                ),
                If[ OptionValue[dimensions]==3,
                    Graphics3D[({Blend[{{0,LightPink},{maxpredvalue,OptionValue[colour]},{2*maxpredvalue,OptionValue[maxcolor]}} ,#[[4]]],PointSize[Large],Point[#[[1;;3]]]}) & /@ plotData,Boxed->False],
					Graphics[({Blend[{{0,LightPink},{maxpredvalue,OptionValue[colour]},{2*maxpredvalue,OptionValue[maxcolor]}} ,#[[3]]],PointSize[OptionValue[pointSize2D]],Point[#[[1;;2]]]}) & /@ plotData]                    
                ]
            ]
        ]
    ] 


PlotEmbryoTll[ds_Dataset,value_String,opts:OptionsPattern[]] :=
    Module[ {plotData,profileData, profileDataSplit,colournucids,offset,maxpredvalue},
        plotData = If[ OptionValue[dimensions]==3,
                       ds[All,{"x","y","z",value}],
						If[OptionValue[view]=="side",
                       ds[All,{"x","z",value}],
					   ds[Select[#z<0&],{"x","y",value}]]
                   ] //Normal //Values;
        offset = If[ OptionValue[dimensions]==3,
                     Min[plotData[[All,3]]]-140, (* was -100 *)
                     Min[plotData[[All,2]]]-140
                 ];
		maxpredvalue=OptionValue[maxvalue];
        If[ Length[plotData]>6078,
            Print["More than one cohort."],
            If[ Length[OptionValue[stripes]]>0,
                (
                colournucids = nucleiStripesDS[Select[MemberQ[OptionValue[stripes], #stripe] &], "nucleus_id"] //  Normal;
				profileData = ds[Select[Abs[#z]<20&],{#x,#[value]*120+offset,#["nucleus_id"]}&];
				profileDataSplit = GatherBy[profileData,MemberQ[colournucids,#[[3]]]&];
                Show[
				  Graphics[({Blend[{{0,White},{0.03,Lighter[Red,0.95]},{0.1,Lighter[Red,0.2]},{maxpredvalue,OptionValue[colour]},{2*maxpredvalue,OptionValue[maxcolor]}} ,#[[3]]],PointSize[OptionValue[pointSize2D]],Point[#[[1;;2]]]}) & /@ plotData],                
                  ListPlot[profileDataSplit[[1,All,1;;2]],FilterRules[{opts},Options[ListPlot]],Axes->False,PlotStyle->Directive[Black,PointSize[Small]],PlotRange->Full],
                  ListPlot[profileDataSplit[[2,All,1;;2]],FilterRules[{opts},Options[ListPlot]],Axes->False,PlotStyle->Directive[Darker[Green],PointSize[Small]],PlotRange->Full]
                ]	
                ),
                If[ OptionValue[dimensions]==3,
                    Graphics3D[({Blend[{{0,White},{0.03,Lighter[Red,0.95]},{0.1,Lighter[Red,0.2]},{maxpredvalue,OptionValue[colour]},{2*maxpredvalue,OptionValue[maxcolor]}},#[[4]]],PointSize[Large],Point[#[[1;;3]]]}) & /@ plotData,Boxed->False],
	(*(* more subtle *)   Graphics[({Blend[{{0,White},{0.03,Lighter[Red,0.95]},{0.1,Lighter[Red,0.6]},{maxpredvalue,OptionValue[colour]},{2*maxpredvalue,OptionValue[maxcolor]}} ,#[[3]]],PointSize[OptionValue[pointSize2D]],Point[#[[1;;2]]]}) & /@ plotData]                    *)
					Graphics[({Blend[{{0,White},{0.03,Lighter[Red,0.95]},{0.1,Lighter[Red,0.2]},{maxpredvalue,OptionValue[colour]},{2*maxpredvalue,OptionValue[maxcolor]}} ,#[[3]]],PointSize[OptionValue[pointSize2D]],Point[#[[1;;2]]]}) & /@ plotData]                    
                ]
            ]
        ]
    ]


PlotModel[data_Dataset,model_FittedModel,opts:OptionsPattern[]] :=
    PlotModel[data,model["BestFit"],opts];
PlotModel[data_Dataset,bestfit_,opts:OptionsPattern[]] :=
    Module[ {predfunction,regs,regSymbols,newdata,lambda,pred,initpred,reprules,maxpredvalue},
		(* the regulators are assumed to be in the first row of the Dataset *)
		predfunction=bestfit /. reg_Symbol/;MemberQ[Normal[data[1,Keys]],ToString[reg]] :>Slot[ToString[reg]]; 
		maxpredvalue=bestfit /. Times[numerator_.,1/(1+Exp[__])] :>numerator;
		newdata=data[All,Function[Evaluate[<|#,"pred"->predfunction|>]]];
        Switch[OptionValue[plotType],
        	"embryo",PlotEmbryo[newdata,"pred",FilterRules[{opts}, Options[PlotEmbryo]],maxvalue->maxpredvalue],
        	"dataset",newdata,
        	_,Message[PlotModel::plotType,OptionValue[plotType]]]
    ] 

PlotModel[flyresults_,opts:OptionsPattern[]] :=
	PlotModel[flyresults["embryo"],flyresults["model"],opts]


PlotMutantSeries[flyResults_, probe_, opts : OptionsPattern[]] := 
  GraphicsRow[Table[PlotModel[MutantEmbryo[flyResults["embryo"], probe, Slot[probe]*scale],flyResults["model"], FilterRules[{opts},Options[PlotModel]]], {scale, OptionValue[SeriesRange]}],ImageSize -> 1500]

PlotTALE[ds_Dataset,fitFunOrModel_, range_, opts : OptionsPattern[]] := Module[{plots},
	plots=Table[PlotModel[ds,ChangePropensity[fitFunOrModel, offset], FilterRules[{opts},Options[PlotModel]]], {offset,range}];
	If[OptionValue[CompositeOrientation]=="horizontal",
		GraphicsRow[plots,ImageSize->OptionValue["CompositeImageSize"]],
		GraphicsColumn[plots,ImageSize->OptionValue["CompositeImageSize"]]
	]]

		
PlotTALE[flyResults_, range_, opts : OptionsPattern[]] := PlotTALE[flyResults["embryo"],flyResults["model"],range,opts]


 
End[] (* End Private Context *)   
EndPackage[]
