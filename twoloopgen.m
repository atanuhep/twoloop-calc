(* ::Package:: *)

(* ::Input:: *)
(*Z->Z Boson*)


description = "El Ael -> Mu Amu, QED, total cross section, tree"; 
If[$FrontEnd === Null, $FeynCalcStartupMessages = False; Print[description]; ]; 
If[$Notebooks === False, $FeynCalcStartupMessages = False];
$LoadAddOns = {"FeynArts"}; 
Get["FeynCalc`"];
$FAVerbose = 0; 
FCCheckVersion[9, 3, 0];


MakeBoxes[q1,TraditionalForm]:="\*SubscriptBox[q,1]";
MakeBoxes[q2,TraditionalForm]:="\*SubscriptBox[q,2]";
MakeBoxes[k1,TraditionalForm]:="\*SubscriptBox[k,1]";
MakeBoxes[k2,TraditionalForm]:="\*SubscriptBox[k,2]";
MakeBoxes[t1,TraditionalForm]:="\*SubscriptBox[t,1]";
MakeBoxes[t2,TraditionalForm]:="\*SubscriptBox[t,2]";
MakeBoxes[t3,TraditionalForm]:="\*SubscriptBox[t,3]";
MakeBoxes[t4,TraditionalForm]:="\*SubscriptBox[t,4]";
MakeBoxes[t5,TraditionalForm]:="\*SubscriptBox[t,5]";


(*DoPaint[diags_, file_] := Paint[diags, ColumnsXRows -> {3,1}, ImageSize -> {1366,768},
	DisplayFunction -> (Export[ToFileName["/home/atanu/Dropbox\Mathematica\diagrams", file <> ".pdf"],#]&)]*)


SetOptions[InsertFields, Model -> "SM", InsertionLevel -> {Particles}];
SetOptions[Paint, PaintLevel -> {Particles}, SheetHeader -> None, Numbering -> Simple];
SetOptions[FCFAConvert, IncomingMomenta -> {q}, OutgoingMomenta -> {q}, LoopMomenta -> {k1,k2},
	UndoChiralSplittings -> True, ChangeDimension -> D, List -> False, SMP -> True,
	Contract -> True, DropSumOver -> True];
process = {V[2]}->{V[2]};



topslist =List @@ CreateTopologies[2, 1 -> 1, ExcludeTopologies -> {Tadpoles}];
tops = CreateTopologies[2, 1 -> 1, ExcludeTopologies -> {Tadpoles}];
ins = InsertFields[TopologyList[topslist[[4]], topslist[[5]]], process, ExcludeParticles ->{U,S}];
(*DoPaint[ins, "oneloop" ];  U,S,V[1],V[3],F[1],F[2,{2}],F[2,{3}],F[3],F[4]
DiagramExtract[tops,5][[1]]
*)
Paint[ins, ColumnsXRows -> {12,3}, ImageSize -> {1000, 256}];



nod = 141; (* No. of Diagrams *)
Do[
amp1[i]= FCFAConvert[CreateFeynAmp[DiagramExtract[ins,i], Truncated -> True], LorentzIndexNames -> {\[Mu], \[Nu], \[Rho], \[Sigma]}],
 {i,1,nod}];
 
Do[
amp2[ii] = DiracSubstitute67[amp1[ii]]/.DiracGamma[5]->0;
amp2[ii] = DiracSimplify[Contract[amp2[ii]*(MTD[\[Mu],\[Nu]] - (FVD[q,\[Mu]]FVD[q,\[Nu]])/SPD[q]) 1/(D-1)]]/.{SMP["m_e"] -> 0,
	SMP["m_mu"] -> 0, SMP["m_tau"] -> 0, SMP["m_u"] -> 0, SMP["m_s"] -> 0, SMP["m_b"] -> 0, SMP["m_d"] -> 0,SMP["m_c"] -> 0,
	SMP["m_t"] -> mt, SMP["m_Z"] -> mZ, SMP["m_W"] -> mW, SMP["m_H"] -> mH, SMP["sin_W"] -> sw,  SMP["cos_W"] -> cw, SMP["e"]->e}//Simplify,
	{ii,1,nod}]



(*
n = 2;
exp = amp2[n]
Paint[DiagramExtract[ins,n],ColumnsXRows -> {1,1}, ImageSize -> {512, 100}]


dlistFA = Cases2[exp, PropagatorDenominator];
dlist = Expand[Join[Cases[dlistFA, PropagatorDenominator[Momentum[mom_,_],mass_]-> mom^2 - mass^2],
		Cases[dlistFA, PropagatorDenominator[sign_*Momentum[mom_,_],mass_]-> mom^2 - mass^2]]];
mplist = Join[Cases[dlistFA, PropagatorDenominator[Momentum[mom_,_],mass_]-> mass],
		Cases[dlistFA, PropagatorDenominator[sign_*Momentum[mom_,_],mass_]-> mass]]
clist = Join[Cases[dlistFA, PropagatorDenominator[Momentum[mom_,_],mass_]-> mom^2],
		Cases[dlistFA, PropagatorDenominator[sign_*Momentum[mom_,_],mass_]-> mom^2]]
mlist = mplist/.{mt -> t, mZ -> z, mW -> w, mH -> h};
dlist2 = {};
Do[dlist2 = Append[dlist2, dlist[[i]] == ToExpression["t"<>ToString[i]]],{i,1,Length[mplist]}];
sol =Solve2[dlist2,{k1^2, k2^2, k1 q, k2 q, k1 k2}]
dlist2*)


Do[
dlistFA = Cases2[amp2[i], PropagatorDenominator];
clist = Expand[Join[Cases[dlistFA, PropagatorDenominator[Momentum[mom_,_],mass_]-> mom^2],
		Cases[dlistFA, PropagatorDenominator[sign_*Momentum[mom_,_],mass_]-> mom^2]]];
amp[i] = amp2[i];
If[ContainsAll[clist, Expand[{(k1+q)^2}]] == True, amp[i] = amp2[i]/.k1->-k1];
If[ContainsAll[clist, Expand[{(k1+k2)^2, (k1-q)^2}]]||ContainsAll[clist, Expand[{(k1-q)^2, (k2+q)^2}]]||
	ContainsExactly[clist, Expand[{k1^2,k2^2,(k1+k2)^2}]] == True, amp[i] = amp2[i]/.k2->-k2];
If[ContainsAll[clist, Expand[{(k1-q)^2, (-k1+k2+q)^2}]] == True, amp[i] = amp2[i]/.{k2 -> k2-q}//ScalarProductExpand];
If[ContainsAll[clist, Expand[{(q+k1+k2)^2}]] == True, amp[i] = amp2[i]/.{k1->k1-q, k2->-k2}//ScalarProductExpand],
{i,1,nod}]
amp[100]


(* ::Text:: *)
(*######*)


(*
n=2
exp1 = amp[n];
dlistFA = Cases2[exp1, PropagatorDenominator];
clist2 = Expand[Join[Cases[dlistFA, PropagatorDenominator[Momentum[mom_,_],mass_]-> mom^2],
		Cases[dlistFA, PropagatorDenominator[sign_*Momentum[mom_,_],mass_]-> mom^2]]]
masslist = Join[Cases[dlistFA, PropagatorDenominator[Momentum[mom_,_],mass_]-> mass],
		Cases[dlistFA, PropagatorDenominator[sign_*Momentum[mom_,_],mass_]-> mass]]
clist1 = {k1^2, k2^2, (k1-k2)^2, (k1-q)^2, (k2-q)^2}//Expand
mplist = {0,0,0,0,0};
Do[
If[ContainsAll[clist2, {clist1[[i]]}]==True, Do[If[clist2[[j]]==clist1[[i]], s=j],{j,1,Length[clist2]}];
 mplist[[i]] =  masslist[[s]]],
{i,1,Length[clist1]}];
mplist
mlist = mplist/.{mt -> t, mZ -> z, mW -> w, mH -> h}
*)


stringjoin[I_,mlist_]:= 
Module[{strexp = I},
Do[strexp = ToString[strexp]<>ToString[mlist[[i]]],{i,1,Length[mlist]}];
strexp
]
stringexp[I_,mlist_]:= 
Module[{strexp = I},
Do[strexp = ToString[strexp]<>ToString[mlist[[i]]],{i,1,Length[mlist]}];
ToExpression[strexp]
]


\quad


amplist2 = {};
clist1 = {k1^2, k2^2, (k1-k2)^2, (k1-q)^2, (k2-q)^2}//Expand;

Do[
If[amp[ii]==0, Continue[]];
dlistFA = Cases2[amp[ii], PropagatorDenominator];
clist2 = Expand[Join[Cases[dlistFA, PropagatorDenominator[Momentum[mom_,_],mass_]-> mom^2],
		Cases[dlistFA, PropagatorDenominator[sign_*Momentum[mom_,_],mass_]-> mom^2]]];
masslist = Join[Cases[dlistFA, PropagatorDenominator[Momentum[mom_,_],mass_]-> mass],
		Cases[dlistFA, PropagatorDenominator[sign_*Momentum[mom_,_],mass_]-> mass]];

mplist = {0,0,0,0,0};

Do[

If[ContainsAll[clist2, {clist1[[i]]}]==True, Do[If[clist2[[j]]==clist1[[i]], s=j],{j,1,Length[clist2]}];
mplist[[i]] =  masslist[[s]]],
{i,1,Length[clist1]}]; (* Rearrange the masses*)

If[ContainsExactly[{mplist[[1]],mplist[[4]],mplist[[5]]},{mt}]==True, mplist[[2]] = mt]; 
If[mplist[[1]]==mw;mplist[[4]]==mw;mplist[[5]]==mt, mplist[[2]] = mt];  (* taking nonzero mass for powerless propagator *)

mlist = mplist/.{mt -> t, mZ -> z, mW -> w, mH -> h};
If[ContainsAny[amplist2,{stringjoin[I,mlist]}]== False, ToExpression[ToString[stringjoin[f,mlist]]<>"="<>"0"]];

amplist2 = Append[amplist2, stringjoin[I,mlist]];
ampsimp = Simplify[ToExpression[ToString[stringjoin[f,mlist]]<>"+="<>"amp[ii]"]];
ToExpression[ToString[stringjoin[f,mlist]]<>"="<>"ampsimp"],
{ii,1,nod}]
intlist = Tally[amplist2]
{intlist1,intlist2};
(*fw00ww
mplist
masslist
clist2*)
f00000


ftt0tt


Sum[intlist[[i]][[2]],{i,Length[intlist]}]


(*
fermintlist = {};
If[ContainsAny[Delete[Characters[ToString[intlist[[1]][[1]]]],1],{"t","0"}]== True, fermintlist = Append[fermintlist, intlist[[1]][[1]]]]

step1 = ampr[1]/.{FeynAmpDenominator[_,_,_,PropagatorDenominator[Momentum[k2-q,_],_],_] -> 1/(t1*t2*t3*t4*t5),
	FeynAmpDenominator[_,_,PropagatorDenominator[Momentum[k1-q,_],_],PropagatorDenominator[Momentum[k1-q,_],_],_] -> 1/(t1*t3*t4^2*t5),
	FeynAmpDenominator[_]*FeynAmpDenominator[_,_,_]->1/(t1*t2*t4^2),FeynAmpDenominator[_,_,_,_]->1/(t1*t2*t3*t4),
	FeynAmpDenominator[PropagatorDenominator[Momentum[k1,_],_]]*FeynAmpDenominator[PropagatorDenominator[Momentum[k2,_],_]]->1/(t1*t2),
	FeynAmpDenominator[PropagatorDenominator[Momentum[k2,_],_]]*FeynAmpDenominator[_,_]->1/(t1*t2*t4), 
	FeynAmpDenominator[_,_,_]->1/(t1*t2*t3),FeynAmpDenominator[PropagatorDenominator[Momentum[_,_],_],_]*
	FeynAmpDenominator[PropagatorDenominator[_,_],_] ->1/(t1*t2*t4*t5)}
step1 = step1/.FeynAmpDenominator[PropagatorDenominator[Momentum[q_,_],qm_]]->1/(q^2-qm^2)

dlist2 = {};
Do[dlist2 = Append[dlist2, clist2[[i]] - ToExpression["m"<>ToString[i]] == ToExpression["t"<>ToString[i]]],{i,1,Length[mplist]}];
sol =Solve2[dlist2,{k1^2, k2^2, k1 q, k2 q, k1 k2}]
*)


Do[
m = Delete[Characters[ToString[intlist[[i]][[1]]]],1];
mp = m/.{"t" -> mt, "z" -> mZ, "w" -> mW, "h" -> mH, "0"->0};
step1 = stringexp[f,m]/.{FeynAmpDenominator[_,_,_,PropagatorDenominator[Momentum[k2-q,_],_],_] -> 1/(t1*t2*t3*t4*t5),
	FeynAmpDenominator[_,_,PropagatorDenominator[Momentum[k1-q,_],_],PropagatorDenominator[Momentum[k1-q,_],_],_] -> 1/(t1*t3*t4^2*t5),
	FeynAmpDenominator[_]*FeynAmpDenominator[_,_,_]->1/(t1*t2*t4^2),FeynAmpDenominator[_,_,_,_]->1/(t1*t2*t3*t4),
	FeynAmpDenominator[PropagatorDenominator[Momentum[k1,_],_]]*FeynAmpDenominator[PropagatorDenominator[Momentum[k2,_],_]]->1/(t1*t2),
	FeynAmpDenominator[PropagatorDenominator[Momentum[k2,_],_]]*FeynAmpDenominator[_,_]->1/(t1*t2*t4), 
	FeynAmpDenominator[_,_,_]->1/(t1*t2*t3),FeynAmpDenominator[PropagatorDenominator[Momentum[_,_],_],_]*
	FeynAmpDenominator[PropagatorDenominator[_,_],_] ->1/(t1*t2*t4*t5)};
step1 = step1/.{FeynAmpDenominator[PropagatorDenominator[Momentum[q_,_],qm_]]->1/(q^2-qm^2),
	FeynAmpDenominator[PropagatorDenominator[sign_ Momentum[q_,_],qm_]]->1/(q^2-qm^2)};
step2 = step1/.{Pair[Momentum[k1,D], Momentum[k1,D]] -> t1 + mp[[1]]^2, Pair[Momentum[k2,D], Momentum[k2,D]] -> t2 + mp[[2]]^2, 
	Pair[Momentum[k1,D], Momentum[k2,D]] -> (t1+t2-t3+mp[[1]]^2+mp[[2]]^2-mp[[3]]^2)/2, Pair[Momentum[k1,D],
	Momentum[q,D]] -> (t1-t4+q^2+mp[[1]]^2-mp[[4]]^2)/2, Pair[Momentum[k2,D], Momentum[q,D]] -> (t2-t5+q^2+mp[[2]]^2-mp[[5]]^2),
	Pair[Momentum[q,D],Momentum[q,D]]-> q^2};
step3 = t1^#*t2^#*t3^#*t4^#*t5^#*step2//Expand;
step4 = step3/.t1^n1_*t2^n2_*t3^n3_*t4^n4_*t5^n5_ -> stringexp[I,m][-n1,-n2,-n3,-n4,-n5]/.#->0;
int[i] = step4//Simplify;
ToExpression[ToString[stringjoin[F,m]]<>"="<>"int[i]"],
{i,1,Length[intlist]}]

Feynints[int_]:=
Module[{},
m = Delete[Characters[ToString[int]],1];
list = Cases2[stringexp[F,m], stringexp[I,m]]
]



step4


Length[intlist];
intlist
Feynints[Ittztt]
Fttztt
fttztt










