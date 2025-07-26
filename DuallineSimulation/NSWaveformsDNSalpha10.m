(* ::Package:: *)

(*fundamental physical parameters*)


bfun[Ps_,\[Gamma]_]:=2*\[Pi]/Ps*Cos[ArcTan[\[Gamma]]]
afun[\[Epsilon]_,\[Gamma]_,Ps_]:=(bfun[Ps,\[Gamma]]*\[Gamma])/(1-\[Epsilon])
\[CapitalOmega]pfun[\[Epsilon]_,\[Gamma]_,\[Kappa]_,Ps_]:=(\[Pi]*bfun[Ps,\[Gamma]])/(2*EllipticK[16 \[Kappa] \[Gamma]^2]) Sqrt[\[Epsilon]^2/((1-\[Epsilon]) (1-\[Epsilon]+16 \[Kappa] (1-\[Epsilon])))]
\[Alpha][\[Gamma]_,\[Kappa]_]:=Re[EllipticF[ArcSin[I/\[Gamma]],16 \[Kappa] \[Gamma]^2]/(2 I EllipticK[16 \[Kappa] \[Gamma]^2])]
q[\[Gamma]_,\[Kappa]_]:=E^{-Pi EllipticK[1-16 \[Kappa] \[Gamma]^2]/EllipticK[16 \[Kappa] \[Gamma]^2]}
ITheta[\[Gamma]_,\[Kappa]_]:=I EllipticThetaPrime[4,I Pi \[Alpha][\[Gamma],\[Kappa]],q[\[Gamma],\[Kappa]]]/EllipticTheta[4,I Pi \[Alpha][\[Gamma],\[Kappa]],q[\[Gamma],\[Kappa]]]
\[CapitalOmega]rfun[\[Epsilon]_,\[Gamma]_,\[Kappa]_,Ps_]:=Re[(bfun[Ps,\[Gamma]] Sqrt[1+\[Gamma]^2])/(1-\[Epsilon])+\[CapitalOmega]pfun[\[Epsilon],\[Gamma],\[Kappa],Ps]*ITheta[\[Gamma],\[Kappa]]-\[CapitalOmega]pfun[\[Epsilon],\[Gamma],\[Kappa],Ps]][[1]]


M[m1_,m2_]:=m1+m2
\[Mu][m1_,m2_]:=m1*m2/M[m1,m2]
a[m1_,m2_,porb_,G_]:=((G*(m1+m2)*porb^2)/(4*\[Pi]^2))^(1/3)
Spin[I3_,Ps_,\[Gamma]_]:=I3*bfun[Ps,\[Gamma]]*Sqrt[1+\[Gamma]^2]
Lorb[m1_,m2_,porb_,e_,G_]:=Sqrt[G*M[m1,m2]*\[Mu][m1,m2]^2*a[m1,m2,porb,G]*(1-e^2)]
\[Theta]L[I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_]:=ArcSin[Spin[I3,Ps,\[Gamma]]*Sin[\[Theta]S]/Lorb[m1,m2,porb,e,G]]
Jtot[I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_]:=Lorb[m1,m2,porb,e,G]*Cos[\[Theta]L[I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G]]+Spin[I3,Ps,\[Gamma]]*Cos[\[Theta]S]
\[CapitalOmega]precess[I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_,c_]:=(G*Jtot[I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G])/(2*c^2*a[m1,m2,porb,G]^3*(1-e^2)^(3/2))*(1+(3 M[m1,m2])/m1)


(*waveforms from spinning NS under spin-orbit coupling*)


hplus1Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=(1-8 \[Kappa])/(4 c^4 r) b^2 G I3 \[Gamma] \[Epsilon] (Cos[t (\[CapitalOmega]p+\[CapitalOmega]r)] (Sin[2 \[Theta]S] (-(3+Cos[2 \[Iota]]) Cos[2 t \[CapitalOmega]pre]+6 Sin[\[Iota]]^2)+4 Cos[2 \[Theta]S] Cos[t \[CapitalOmega]pre] Sin[2 \[Iota]])+2 (-2 Cos[\[Theta]S] Sin[2 \[Iota]] Sin[t \[CapitalOmega]pre]+(3+Cos[2 \[Iota]]) Sin[\[Theta]S] Sin[2 t \[CapitalOmega]pre]) Sin[t (\[CapitalOmega]p+\[CapitalOmega]r)])
hplus2Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=-(1-16 \[Kappa])(1/(c^4 r))16 b^2 G I3 \[Epsilon] \[Kappa] (Cos[\[Theta]S/2]^4 (3+Cos[2 \[Iota]]) Cos[2 t (\[CapitalOmega]pre+\[CapitalOmega]r)]+(3+Cos[2 \[Iota]]) Cos[2 t (\[CapitalOmega]pre-\[CapitalOmega]r)] Sin[\[Theta]S/2]^4+Cos[2 t \[CapitalOmega]r] (3 Sin[\[Theta]S]^2 Sin[\[Iota]]^2+Cos[t \[CapitalOmega]pre] Sin[2 \[Theta]S] Sin[2 \[Iota]])-2 Sin[\[Theta]S] Sin[2 \[Iota]] Sin[t \[CapitalOmega]pre] Sin[2 t \[CapitalOmega]r])
hcross1Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=(1-8 \[Kappa])/(c^4 r) b^2 G I3 \[Gamma] \[Epsilon] (Cos[t (\[CapitalOmega]p+\[CapitalOmega]r)] (2 Cos[2 \[Theta]S] Sin[\[Iota]] Sin[t \[CapitalOmega]pre]-Cos[\[Iota]] Sin[2 \[Theta]S] Sin[2 t \[CapitalOmega]pre])+2 (-Cos[\[Iota]] Cos[2 t \[CapitalOmega]pre] Sin[\[Theta]S]+Cos[\[Theta]S] Cos[t \[CapitalOmega]pre] Sin[\[Iota]]) Sin[t (\[CapitalOmega]p+\[CapitalOmega]r)])
hcross2Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=(1-16 \[Kappa])/(c^4 r) 16 b^2 G I3 \[Epsilon] \[Kappa] (Cos[2 t \[CapitalOmega]r] (-2 Sin[2 \[Theta]S] Sin[\[Iota]] Sin[t \[CapitalOmega]pre]-(3+Cos[2 \[Theta]S]) Cos[\[Iota]] Sin[2 t \[CapitalOmega]pre])-4 (Cos[\[Theta]S] Cos[\[Iota]] Cos[2 t \[CapitalOmega]pre]+Cos[t \[CapitalOmega]pre] Sin[\[Theta]S] Sin[\[Iota]]) Sin[2 t \[CapitalOmega]r])
hplus3Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=(1+(64 \[Kappa]^2)/\[Gamma]^2)/(4 c^4 r) b^2 G I3 \[Gamma]^2 \[Epsilon] (4 Cos[\[Theta]S/2]^4 (3+Cos[2 \[Iota]]) Cos[2 t (\[CapitalOmega]p+\[CapitalOmega]pre+\[CapitalOmega]r)]+4 (3+Cos[2 \[Iota]]) Cos[2 t (\[CapitalOmega]p-\[CapitalOmega]pre+\[CapitalOmega]r)] Sin[\[Theta]S/2]^4+12 Cos[2 t (\[CapitalOmega]p+\[CapitalOmega]r)] Sin[\[Theta]S]^2 Sin[\[Iota]]^2+Sin[2 \[Iota]] (-2 Sin[\[Theta]S+t \[CapitalOmega]pre-2 t (\[CapitalOmega]p+\[CapitalOmega]r)]+Sin[2 \[Theta]S+t \[CapitalOmega]pre-2 t (\[CapitalOmega]p+\[CapitalOmega]r)]-2 Sin[\[Theta]S-t \[CapitalOmega]pre+2 t (\[CapitalOmega]p+\[CapitalOmega]r)]+Sin[2 \[Theta]S-t \[CapitalOmega]pre+2 t (\[CapitalOmega]p+\[CapitalOmega]r)]+2 Sin[\[Theta]S+t \[CapitalOmega]pre+2 t (\[CapitalOmega]p+\[CapitalOmega]r)]+Sin[2 \[Theta]S+t \[CapitalOmega]pre+2 t (\[CapitalOmega]p+\[CapitalOmega]r)]+2 Sin[\[Theta]S-t (2 \[CapitalOmega]p+\[CapitalOmega]pre+2 \[CapitalOmega]r)]+Sin[2 \[Theta]S-t (2 \[CapitalOmega]p+\[CapitalOmega]pre+2 \[CapitalOmega]r)]))
hcross3Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=(1+(64 \[Kappa]^2)/\[Gamma]^2)/(c^4 r) b^2 G I3 \[Gamma]^2 \[Epsilon] (Cos[2 t (\[CapitalOmega]p+\[CapitalOmega]r)] (2 Sin[2 \[Theta]S] Sin[\[Iota]] Sin[t \[CapitalOmega]pre]+(3+Cos[2 \[Theta]S]) Cos[\[Iota]] Sin[2 t \[CapitalOmega]pre])+4 (Cos[\[Theta]S] Cos[\[Iota]] Cos[2 t \[CapitalOmega]pre]+Cos[t \[CapitalOmega]pre] Sin[\[Theta]S] Sin[\[Iota]]) Sin[2 t (\[CapitalOmega]p+\[CapitalOmega]r)])


hplus4Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=1/(2 c^4 r) 5 b^2 G I3 \[Gamma] \[Epsilon] \[Kappa] (Cos[t (\[CapitalOmega]p-\[CapitalOmega]r)] (Sin[2 \[Theta]S] (-((3+Cos[2 \[Iota]]) Cos[2 t \[CapitalOmega]pre])+6 Sin[\[Iota]]^2)+4 Cos[2 \[Theta]S] Cos[t \[CapitalOmega]pre] Sin[2 \[Iota]])+2 (2 Cos[\[Theta]S] Sin[2 \[Iota]] Sin[t \[CapitalOmega]pre]-(3+Cos[2 \[Iota]]) Sin[\[Theta]S] Sin[2 t \[CapitalOmega]pre]) Sin[t (\[CapitalOmega]p-\[CapitalOmega]r)])
hplus5Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=1/(2 c^4 r) b^2 G I3 \[Gamma] \[Epsilon] \[Kappa] (4 Cos[\[Theta]S/2] (3+Cos[2 \[Iota]]) Sin[\[Theta]S/2] (Cos[\[Theta]S/2]^2 Cos[t (3 \[CapitalOmega]p+2 \[CapitalOmega]pre+\[CapitalOmega]r)]-Cos[t (3 \[CapitalOmega]p-2 \[CapitalOmega]pre+\[CapitalOmega]r)] Sin[\[Theta]S/2]^2)-6 Cos[t (3 \[CapitalOmega]p+\[CapitalOmega]r)] Sin[2 \[Theta]S] Sin[\[Iota]]^2+4 (Cos[\[Theta]S/2]^2 (1-2 Cos[\[Theta]S]) Cos[t (3 \[CapitalOmega]p+\[CapitalOmega]pre+\[CapitalOmega]r)]+(1+2 Cos[\[Theta]S]) Cos[t (-3 \[CapitalOmega]p+\[CapitalOmega]pre-\[CapitalOmega]r)] Sin[\[Theta]S/2]^2) Sin[2 \[Iota]])
hplus6Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=-(1/(c^4 r))8 b^2 G I3 \[Epsilon] \[Kappa]^2 (-Cos[\[Theta]S-2 (\[Iota]+t (\[CapitalOmega]p+\[CapitalOmega]pre-\[CapitalOmega]r))]-2 Cos[\[Theta]S-2 \[Iota]-t (2 \[CapitalOmega]p+\[CapitalOmega]pre-2 \[CapitalOmega]r)]+Cos[2 \[Theta]S-2 \[Iota]-t (2 \[CapitalOmega]p+\[CapitalOmega]pre-2 \[CapitalOmega]r)]+2 Cos[\[Theta]S+2 \[Iota]-t (2 \[CapitalOmega]p+\[CapitalOmega]pre-2 \[CapitalOmega]r)]-Cos[2 \[Theta]S+2 \[Iota]-t (2 \[CapitalOmega]p+\[CapitalOmega]pre-2 \[CapitalOmega]r)]-2 Cos[\[Theta]S-2 \[Iota]+t (2 \[CapitalOmega]p+\[CapitalOmega]pre-2 \[CapitalOmega]r)]+Cos[2 \[Theta]S-2 \[Iota]+t (2 \[CapitalOmega]p+\[CapitalOmega]pre-2 \[CapitalOmega]r)]+2 Cos[\[Theta]S+2 \[Iota]+t (2 \[CapitalOmega]p+\[CapitalOmega]pre-2 \[CapitalOmega]r)]-Cos[2 \[Theta]S+2 \[Iota]+t (2 \[CapitalOmega]p+\[CapitalOmega]pre-2 \[CapitalOmega]r)]+8 Cos[\[Theta]S/2]^4 (3+Cos[2 \[Iota]]) Cos[2 t (\[CapitalOmega]p-\[CapitalOmega]pre-\[CapitalOmega]r)]+24 Cos[2 t (\[CapitalOmega]p+\[CapitalOmega]pre-\[CapitalOmega]r)] Sin[\[Theta]S/2]^4+12 Cos[2 t (\[CapitalOmega]p-\[CapitalOmega]r)] Sin[\[Theta]S]^2+Sin[2 \[Iota]] (4 Cos[t (2 \[CapitalOmega]p-\[CapitalOmega]pre-2 \[CapitalOmega]r)] (2 Sin[\[Theta]S]+Sin[2 \[Theta]S])+Sin[\[Theta]S-2 t (\[CapitalOmega]p+\[CapitalOmega]pre-\[CapitalOmega]r)])+Cos[2 \[Iota]] ((3-3 Cos[\[Theta]S]+Cos[2 \[Theta]S]) Cos[2 t (\[CapitalOmega]p+\[CapitalOmega]pre-\[CapitalOmega]r)]+Sin[\[Theta]S] (-12 Cos[2 t (\[CapitalOmega]p-\[CapitalOmega]r)] Sin[\[Theta]S]+Sin[2 t (\[CapitalOmega]p+\[CapitalOmega]pre-\[CapitalOmega]r)])))


hcross4Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=1/(c^4 r) 20 b^2 G I3 \[Gamma] \[Epsilon] \[Kappa] (Sin[\[Iota]] (Cos[2 \[Theta]S] Cos[t (\[CapitalOmega]p-\[CapitalOmega]r)] Sin[t \[CapitalOmega]pre]-Cos[\[Theta]S] Cos[t \[CapitalOmega]pre] Sin[t (\[CapitalOmega]p-\[CapitalOmega]r)])+1/2 Cos[\[Iota]] Sin[\[Theta]S] ((1+Cos[\[Theta]S]) Sin[t (\[CapitalOmega]p-2 \[CapitalOmega]pre-\[CapitalOmega]r)]-(-1+Cos[\[Theta]S]) Sin[t (\[CapitalOmega]p+2 \[CapitalOmega]pre-\[CapitalOmega]r)]))
hcross5Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=1/(c^4 r) 2 b^2 G I3 \[Gamma] \[Epsilon] \[Kappa] (Sin[\[Iota]] ((Cos[\[Theta]S]-Cos[2 \[Theta]S]) Sin[t (-3 \[CapitalOmega]p+\[CapitalOmega]pre-\[CapitalOmega]r)]-(Cos[\[Theta]S]+Cos[2 \[Theta]S]) Sin[t (3 \[CapitalOmega]p+\[CapitalOmega]pre+\[CapitalOmega]r)])+Cos[\[Iota]] Sin[\[Theta]S] (-((-1+Cos[\[Theta]S]) Sin[t (3 \[CapitalOmega]p-2 \[CapitalOmega]pre+\[CapitalOmega]r)])+(1+Cos[\[Theta]S]) Sin[t (3 \[CapitalOmega]p+2 \[CapitalOmega]pre+\[CapitalOmega]r)]))
hcross6Func[t_,G_,c_,r_,b_,I3_,\[Epsilon]_, \[Kappa]_,\[Gamma]_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[CapitalOmega]pre_,\[Theta]S_,\[Iota]_]:=1/(c^4 r) 16 b^2 G I3 \[Epsilon] \[Kappa]^2 (4 Sin[2 \[Theta]S] Sin[\[Iota]] (Sin[t (2 \[CapitalOmega]p-\[CapitalOmega]pre-2 \[CapitalOmega]r)]-Sin[t (2 \[CapitalOmega]p+\[CapitalOmega]pre-2 \[CapitalOmega]r)])+8 Sin[\[Theta]S] Sin[\[Iota]] (Sin[t (2 \[CapitalOmega]p-\[CapitalOmega]pre-2 \[CapitalOmega]r)]+Sin[t (2 \[CapitalOmega]p+\[CapitalOmega]pre-2 \[CapitalOmega]r)])+2 Cos[\[Iota]] Sin[\[Theta]S]^2 (-Sin[2 t (\[CapitalOmega]p-\[CapitalOmega]pre-\[CapitalOmega]r)]+Sin[2 t (\[CapitalOmega]p+\[CapitalOmega]pre-\[CapitalOmega]r)])+Csc[\[Iota]] Sin[2 \[Iota]] ((1+Cos[\[Theta]S]) (3+Cos[\[Theta]S]) Sin[2 t (\[CapitalOmega]p-\[CapitalOmega]pre-\[CapitalOmega]r)]-(-3+Cos[\[Theta]S]) (-1+Cos[\[Theta]S]) Sin[2 t (\[CapitalOmega]p+\[CapitalOmega]pre-\[CapitalOmega]r)]))


(*elliptical orbit*)


MA[t_,t0_,porb_]:=2*\[Pi]*(t-t0)/porb
ecos\[Phi]ta[t_,t0_,porb_,e_]:=-e^2+(1-e^2)*Sum[2*BesselJ[n,n e]*Cos[n*MA[t,t0,porb]],{n,1,50}]


g[e_]:=e^(12/19)/(1-e^2)*(1+(121 e^2)/304)^(870/2299)
c0[m1_,m2_,porb0_,e0_,G_]:=a[m1,m2,porb0,G]/g[e0]
Pb[m1_,m2_,porb0_,e0_,e_,G_]:=(c0[m1,m2,porb0,e0,G]*g[e])^(3/2)*((4*\[Pi]^2)/(G*(m1+m2)))^(1/2)


(*Velocity in J-frame*)


\[Theta]Lfun[I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_]:=\[Theta]L[I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G]
\[CapitalOmega]Precesstfun[t_,I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_,c_]:=\[CapitalOmega]precess[I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c]*t
Rx[\[Theta]L_]:={{1,0,0},{0,Cos[\[Theta]L],-Sin[\[Theta]L]},{0,Sin[\[Theta]L],Cos[\[Theta]L]}}
Rz[\[Phi]N_]:={{Cos[\[Phi]N],-Sin[\[Phi]N],0},{Sin[\[Phi]N],Cos[\[Phi]N],0},{0,0,1}}
RL2J[t_,I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_,c_]:=Module[{thetaL,omegaPrecesst,RzPhiN,RxThetaL},
thetaL=\[Theta]Lfun[I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G];omegaPrecesst=\[CapitalOmega]Precesstfun[t,I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c];RzPhiN=Rz[omegaPrecesst];RxThetaL=Rx[thetaL];Dot[RzPhiN,RxThetaL]]
a1[m1_,m2_,porb_,G_]:=((G*porb^2*m2^3/(m1+m2)^2)/(4* \[Pi]^2))^(1/3)
r1[t_,t0_,m1_,m2_,porb_,e_,G_]:=(a1[m1,m2,porb,G]*(1-e^2))/(1+ecos\[Phi]ta[t,t0,porb,e])
rInL[t_,t0_,m1_,m2_,porb_,e_,G_]:={r1[t,t0,m1,m2,porb,e,G]*Cos[2\[Pi]/porb*t],r1[t,t0,m1,m2,porb,e,G]*Sin[2\[Pi]/porb*t],0}
rInJ[t_,t0_,I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_,c_]:=Dot[RL2J[t,I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c],rInL[t,t0,m1,m2,porb,e,G]]
derivaVInJ[t_,t0_,I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_,c_]:=D[rInJ[t,t0,I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c],t]
VInJ[t_,t0_,I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_,c_]=derivaVInJ[t,t0,I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c];


(*Source (NS) Doppler frequency shift*)


nvector=-{0,Sin[\[Iota]],Cos[\[Iota]]};(*from SSB*)
\[CapitalDelta]fSD[t_,f_,\[Iota]_,t0_,I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_,c_]=f*(nvector . VInJ[t,t0,I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c])/c;


(*Detector Doppler frequency shift*)


n0rE=RE*(Sin[\[Lambda]] Sin[\[Delta]d]+Cos[\[Lambda]] Cos[\[Delta]d] Cos[\[Alpha]-\[Phi]r-\[CapitalOmega]rE*t]);
n0rES=RES*(Cos[\[Alpha]] Cos[\[Delta]d] Cos[\[Phi]o+\[CapitalOmega]o*t]+(Cos[\[Epsilon]a] Sin[\[Alpha]] Cos[\[Delta]d]+Sin[\[Epsilon]a] Sin[\[Delta]d])*Sin[\[Phi]o+\[CapitalOmega]o*t]);
\[CapitalDelta]fDD[t_,f_,\[Alpha]_,\[Delta]d_,\[Lambda]_,\[Epsilon]a_,\[CapitalOmega]rE_,\[Phi]r_,\[CapitalOmega]o_,\[Phi]o_,RE_,RES_,c_]=f*D[n0rE+n0rES,t]/c;


(*Doppler frequency*)


\[CapitalOmega]rDoppler[t_,\[CapitalOmega]r_,\[Alpha]_,\[Delta]d_,\[Iota]_,t0_,I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_,c_,\[Lambda]_,\[Epsilon]a_,\[CapitalOmega]rE_,\[Phi]r_,\[CapitalOmega]o_,\[Phi]o_,RE_,RES_]:=\[CapitalOmega]r+\[CapitalDelta]fSD[t,\[CapitalOmega]r,\[Iota],t0,I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c]+\[CapitalDelta]fDD[t,\[CapitalOmega]r,\[Alpha],\[Delta]d,\[Lambda],\[Epsilon]a,\[CapitalOmega]rE,\[Phi]r,\[CapitalOmega]o,\[Phi]o,RE,RES,c]
\[CapitalOmega]pDoppler[t_,\[CapitalOmega]p_,\[Alpha]_,\[Delta]d_,\[Iota]_,t0_,I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_,c_,\[Lambda]_,\[Epsilon]a_,\[CapitalOmega]rE_,\[Phi]r_,\[CapitalOmega]o_,\[Phi]o_,RE_,RES_]:=\[CapitalOmega]p+\[CapitalDelta]fSD[t,\[CapitalOmega]p,\[Iota],t0,I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c]+\[CapitalDelta]fDD[t,\[CapitalOmega]p,\[Alpha],\[Delta]d,\[Lambda],\[Epsilon]a,\[CapitalOmega]rE,\[Phi]r,\[CapitalOmega]o,\[Phi]o,RE,RES,c]
\[CapitalOmega]preDoppler[t_,\[CapitalOmega]pre_,\[Alpha]_,\[Delta]d_,\[Iota]_,t0_,I3_,Ps_,\[Gamma]_,\[Theta]S_,m1_,m2_,porb_,e_,G_,c_,\[Lambda]_,\[Epsilon]a_,\[CapitalOmega]rE_,\[Phi]r_,\[CapitalOmega]o_,\[Phi]o_,RE_,RES_]:=\[CapitalOmega]pre+\[CapitalDelta]fSD[t,\[CapitalOmega]pre,\[Iota],t0,I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c]+\[CapitalDelta]fDD[t,\[CapitalOmega]pre,\[Alpha],\[Delta]d,\[Lambda],\[Epsilon]a,\[CapitalOmega]rE,\[Phi]r,\[CapitalOmega]o,\[Phi]o,RE,RES,c]


(*waveforms with Doppler shift*)


NonAlignedBinaryDopplerFunc[funcName_,t_,G_,c_,r_,I3_,\[Epsilon]_,\[Kappa]_,\[Gamma]_,\[Alpha]_,\[Delta]d_,\[Iota]_,t0_,Ps_,\[Theta]S_,m1_,m2_,porb_,e_,\[Lambda]_,\[Epsilon]a_,\[CapitalOmega]rE_,\[Phi]r_,\[CapitalOmega]o_,\[Phi]o_,RE_,RES_,\[CapitalOmega]r_,\[CapitalOmega]p_]:=Module[{\[CapitalOmega]pre,bf,opDoppler,orDoppler,oPreDoppler},
\[CapitalOmega]pre=\[CapitalOmega]precess[I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c];
bf=bfun[Ps,\[Gamma]];
opDoppler=\[CapitalOmega]pDoppler[t,\[CapitalOmega]p,\[Alpha],\[Delta]d,\[Iota],t0,I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c,\[Lambda],\[Epsilon]a,\[CapitalOmega]rE,\[Phi]r,\[CapitalOmega]o,\[Phi]o,RE,RES];
orDoppler=\[CapitalOmega]rDoppler[t,\[CapitalOmega]r,\[Alpha],\[Delta]d,\[Iota],t0,I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c,\[Lambda],\[Epsilon]a,\[CapitalOmega]rE,\[Phi]r,\[CapitalOmega]o,\[Phi]o,RE,RES];
oPreDoppler=\[CapitalOmega]preDoppler[t,\[CapitalOmega]pre,\[Alpha],\[Delta]d,\[Iota],t0,I3,Ps,\[Gamma],\[Theta]S,m1,m2,porb,e,G,c,\[Lambda],\[Epsilon]a,\[CapitalOmega]rE,\[Phi]r,\[CapitalOmega]o,\[Phi]o,RE,RES];
ToExpression[funcName][t,G,c,r,bf,I3,\[Epsilon],\[Kappa],\[Gamma],orDoppler,opDoppler,oPreDoppler,\[Theta]S,\[Iota]]];


hplus1NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hplus1Func",##]&;
hplus2NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hplus2Func",##]&;
hplus3NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hplus3Func",##]&;
hcross1NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hcross1Func",##]&;
hcross2NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hcross2Func",##]&;
hcross3NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hcross3Func",##]&;


hplus4NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hplus4Func",##]&;
hplus5NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hplus5Func",##]&;
hplus6NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hplus6Func",##]&;
hcross4NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hcross4Func",##]&;
hcross5NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hcross5Func",##]&;
hcross6NonAlignedBinaryDoppler=NonAlignedBinaryDopplerFunc["hcross6Func",##]&;


(*beam pattern functions Fplus and Fcross*)


at=1./16*Sin[2 \[Gamma]o]*(3-Cos[2 \[Lambda]])*(3-Cos[2 \[Delta]d])*Cos[2 (\[Alpha]-\[Phi]r-\[CapitalOmega]rE*t)]-1./4*Cos[2 \[Gamma]o]*Sin[\[Lambda]]*(3-Cos[2 \[Delta]d])*Sin[2 (\[Alpha]-\[Phi]r-\[CapitalOmega]rE*t)]+1./4*Sin[2 \[Gamma]o]*Sin[2 \[Lambda]]*Sin[2 \[Delta]d]*Cos[\[Alpha]-\[Phi]r-\[CapitalOmega]rE*t]-1./2*Cos[2 \[Gamma]o]*Cos[\[Lambda]]*Sin[2 \[Delta]d]*Sin[\[Alpha]-\[Phi]r-\[CapitalOmega]rE*t]+3./4*Sin[2 \[Gamma]o]*Cos[\[Lambda]]^2*Cos[\[Delta]d]^2;
bt=Cos[2 \[Gamma]o]*Sin[\[Lambda]]*Sin[\[Delta]d]*Cos[2 (\[Alpha]-\[Phi]r-\[CapitalOmega]rE*t)]+1./4*Sin[2 \[Gamma]o]*Sin[\[Delta]d]*(3-Cos[2 \[Lambda]])*Sin[2 (\[Alpha]-\[Phi]r-\[CapitalOmega]rE*t)]+Cos[2 \[Gamma]o]*Cos[\[Lambda]]*Cos[\[Delta]d]*Cos[\[Alpha]-\[Phi]r-\[CapitalOmega]rE*t]+1./2*Sin[2 \[Gamma]o]*Sin[2 \[Lambda]]*Cos[\[Delta]d]*Sin[\[Alpha]-\[Phi]r-\[CapitalOmega]rE*t];
Fplus[t_,\[Alpha]_,\[Delta]d_,\[Psi]p_,\[Lambda]_,\[CapitalOmega]rE_,\[Phi]r_,\[Gamma]o_,\[Zeta]_]=Sin[\[Zeta]] (at Cos[2 \[Psi]p]+bt Sin[2 \[Psi]p]);
Fcross[t_,\[Alpha]_,\[Delta]d_,\[Psi]p_,\[Lambda]_,\[CapitalOmega]rE_,\[Phi]r_,\[Gamma]o_,\[Zeta]_]=Sin[\[Zeta]] (bt Cos[2 \[Psi]p]-at Sin[2 \[Psi]p]);


(*non-aligned waveforms with antenna pattern*)


(*Define a general helper function to compute hp and hc*)
computeHPHC[hpFunc_,hcFunc_,params_]:=Module[{hp,hc},hp=hpFunc@@params;hc=hcFunc@@params;{hp,hc}];

(*Generalized function to calculate hp and hc based on the binary index*)
generalNonAlignedBinaryDoppler[index_,t_,G_,c_,r_,I3_,\[Epsilon]_,\[Kappa]_,\[Gamma]_,\[Alpha]_,\[Delta]d_,\[Iota]_,t0_,Ps_,\[Theta]S_,m1_,m2_,porb_,e_,\[Psi]p_,\[Lambda]_,\[Epsilon]a_,\[CapitalOmega]rE_,\[Phi]r_,\[CapitalOmega]o_,\[Phi]o_,RE_,RES_,\[CapitalOmega]r_,\[CapitalOmega]p_,\[Gamma]o_,\[Zeta]_]:=Module[{hp,hc,Fp,Fc,params,hpFunc,hcFunc},
(*Define the parameter list*)params={t,G,c,r,I3,\[Epsilon],\[Kappa],\[Gamma],\[Alpha],\[Delta]d,\[Iota],t0,Ps,\[Theta]S,m1,m2,porb,e,\[Lambda],\[Epsilon]a,\[CapitalOmega]rE,\[Phi]r,\[CapitalOmega]o,\[Phi]o,RE,RES,\[CapitalOmega]r,\[CapitalOmega]p};
(*Select the appropriate function based on the index*)hpFunc=ToExpression["hplus"<>ToString[index]<>"NonAlignedBinaryDoppler"];hcFunc=ToExpression["hcross"<>ToString[index]<>"NonAlignedBinaryDoppler"];
(*Compute hp and hc*){hp,hc}=computeHPHC[hpFunc,hcFunc,params];
(*Calculate Fp and Fc*)Fp=Fplus[t,\[Alpha],\[Delta]d,\[Psi]p,\[Lambda],\[CapitalOmega]rE,\[Phi]r,\[Gamma]o,\[Zeta]];Fc=Fcross[t,\[Alpha],\[Delta]d,\[Psi]p,\[Lambda],\[CapitalOmega]rE,\[Phi]r,\[Gamma]o,\[Zeta]];
(*Return the final result*)hp*Fp+hc*Fc];

(*Define the specific hp and hc functions using the general function*)
h1NonAlignedBinaryDoppler=generalNonAlignedBinaryDoppler[1,##]&;
h2NonAlignedBinaryDoppler=generalNonAlignedBinaryDoppler[2,##]&;
h3NonAlignedBinaryDoppler=generalNonAlignedBinaryDoppler[3,##]&;
h4NonAlignedBinaryDoppler=generalNonAlignedBinaryDoppler[4,##]&;
h5NonAlignedBinaryDoppler=generalNonAlignedBinaryDoppler[5,##]&;
h6NonAlignedBinaryDoppler=generalNonAlignedBinaryDoppler[6,##]&;


(*Waveform expression with simplified variables*)


Gval=6.6743*10^(-11);
cval=299792458.0;
t0val=0;
\[Epsilon]aval=23.5*\[Pi]/180;
\[CapitalOmega]rEval=2 \[Pi]/24/3600.0;
\[Phi]rval=0.0;
\[CapitalOmega]oval=2 \[Pi]/31536000.0;
\[Phi]oval=0.0;
REval=6371*10^3;
RESval=149597871*10^3;
\[Gamma]oval=1.5;
\[Lambda]val=43.8*\[Pi]/180; (*for Cosmic Explorer*)
\[Zeta]val=\[Pi]/2;


defineBinaryDetectedFunction[index_]:=Module[{hpFunc,hcFunc,generalNonAlignedBinaryDoppler},
generalNonAlignedBinaryDoppler[t_,Ps_,I3_,\[Epsilon]_,\[Kappa]_,\[Gamma]_,r_,\[Alpha]_,\[Delta]d_,\[Psi]p_,\[Iota]_,\[CapitalOmega]r_,\[CapitalOmega]p_,porb_,\[Theta]S_,e_,m1_,m2_]:=Module[{hp,hc,Fp,Fc,params},
params={t,Gval,cval,r,I3,\[Epsilon],\[Kappa],\[Gamma],\[Alpha],\[Delta]d,\[Iota],t0val,Ps,\[Theta]S,m1,m2,porb,e,\[Lambda]val,\[Epsilon]aval,\[CapitalOmega]rEval,\[Phi]rval,\[CapitalOmega]oval,\[Phi]oval,REval,RESval,\[CapitalOmega]r,\[CapitalOmega]p};
hpFunc=ToExpression["hplus"<>ToString[index]<>"NonAlignedBinaryDoppler"];
hcFunc=ToExpression["hcross"<>ToString[index]<>"NonAlignedBinaryDoppler"];
{hp,hc}=computeHPHC[hpFunc,hcFunc,params];
Fp=Fplus[t,\[Alpha],\[Delta]d,\[Psi]p,\[Lambda]val,\[CapitalOmega]rEval,\[Phi]rval,\[Gamma]oval,\[Zeta]val];
Fc=Fcross[t,\[Alpha],\[Delta]d,\[Psi]p,\[Lambda]val,\[CapitalOmega]rEval,\[Phi]rval,\[Gamma]oval,\[Zeta]val];
hp*Fp+hc*Fc];
generalNonAlignedBinaryDoppler];

h1inBinaryDetected=defineBinaryDetectedFunction[1];
h2inBinaryDetected=defineBinaryDetectedFunction[2];
h3inBinaryDetected=defineBinaryDetectedFunction[3];
h4inBinaryDetected=defineBinaryDetectedFunction[4];
h5inBinaryDetected=defineBinaryDetectedFunction[5];
h6inBinaryDetected=defineBinaryDetectedFunction[6];


(*generate equatorial ellipticity from exponential distribution*)


generate\[Epsilon]qSamples[\[Epsilon]qmax_,\[Epsilon]qmean_,numSamples_]:=Module[{(*Local variables*)\[Tau]sol,tau,p,cdf,invCdf,\[Epsilon]qsamples},
(*Define the mean value function based on \[Tau]*)mean\[Epsilon]q[\[Tau]_]:=\[Tau]-(\[Epsilon]qmax/(Exp[\[Epsilon]qmax/\[Tau]]-1));
(*Solve for \[Tau] given the target mean value \[Epsilon]qmean*)\[Tau]sol=FindRoot[mean\[Epsilon]q[\[Tau]]==\[Epsilon]qmean,{\[Tau],1.0}];tau=\[Tau]/. \[Tau]sol;
(*Define the probability distribution function*)p[\[Epsilon]q_]:=Exp[-\[Epsilon]q/tau]/(tau*(1-Exp[-\[Epsilon]qmax/tau]));
(*Define the CDF of the distribution*)cdf[\[Epsilon]q_]:=Integrate[p[\[Epsilon]],{\[Epsilon],0,\[Epsilon]q}];
(*cdf[\[Epsilon]q_]:=(1-Exp[-\[Epsilon]q/tau])/(1-Exp[-\[Epsilon]qmax/tau]);*)
(*Find the inverse of the CDF*)invCdf=InverseFunction[cdf];
(*Generate samples using the inverse CDF method*)\[Epsilon]qsamples=Table[invCdf[RandomReal[]],{numSamples}];
(*Return samples*)\[Epsilon]qsamples]


sampleSize=35;
targetFolder="DuallineSimulation";
If[SameQ[$OperatingSystem,"Windows"],basePath=FileNameJoin[{$HomeDirectory,"Desktop",targetFolder}]<>"\\",basePath=FileNameJoin[{$HomeDirectory,"Desktop",targetFolder}]<>"/"];

Msun=1.988409870698051*10^30;
M1List=Import[basePath <> "m1_alpha10.txt","List"]*Msun;
M2List=Import[basePath <> "m2_alpha10.txt","List"]*Msun;
M1kde=SmoothKernelDistribution[M1List];
m1List=RandomVariate[M1kde,sampleSize];
M2kde=SmoothKernelDistribution[M2List];
m2List=RandomVariate[M2kde,sampleSize];

kpcinm=3.086*10^19;
DNSDismFreEcc=Import[basePath <> "DNS_DismFreEcc_LISA_alpha10.txt","Table"];
DNSDismFreEcckde=SmoothKernelDistribution[DNSDismFreEcc];
DNSDismFreEccsamples=RandomVariate[DNSDismFreEcckde,sampleSize];
rList=DNSDismFreEccsamples[[All,1]]*kpcinm;
PorbList=1000/DNSDismFreEccsamples[[All,2]];
eList=DNSDismFreEccsamples[[All,3]];

\[Iota]List=ArcCos[RandomReal[{-1,1},sampleSize]];
\[Psi]List=RandomReal[{0,\[Pi]},sampleSize];
RADEC=Import[basePath <> "RA_DEC_radian_LISA_alpha10.txt","Table"];
RADECkde=SmoothKernelDistribution[RADEC];
RADECsamples=RandomVariate[RADECkde,sampleSize];
\[Alpha]List=RADECsamples[[All,1]];
\[Delta]dList=RADECsamples[[All,2]];

PsList=Table[xVal=eList[[i]];yFit=19.7363+162.954*xVal;lowerNoiseLimit=10-yFit;
noiseDist=TruncatedDistribution[{lowerNoiseLimit,Infinity},NormalDistribution[0,3*11.1145]];
ySim=yFit+RandomVariate[noiseDist];
ySim,{i,sampleSize}]/1000;

(*I3List=10^RandomReal[{38,38+Log10[3]},sampleSize];
Q22List=10^RandomReal[{Log10[4.4*10^28],Log10[2.4*10^32]},sampleSize];
generateEpsilon[\[CurlyEpsilon]q_] := Module[{\[CurlyEpsilon]},While[True,\[CurlyEpsilon] = 10^RandomReal[{-9, -5}];If[\[CurlyEpsilon]>\[CurlyEpsilon]q && \[CurlyEpsilon]q<(1-\[CurlyEpsilon])(\[CurlyEpsilon]-\[CurlyEpsilon]q), Return[\[CurlyEpsilon]]]]]
\[Gamma]List=10^RandomReal[{-3,Log10[0.05]},sampleSize];

I3List=RandomReal[{10^38,3*10^38},sampleSize];
Q22List=RandomReal[{4.4*10^28,2.4*10^32},sampleSize];
generateEpsilon[\[CurlyEpsilon]q_] := Module[{\[CurlyEpsilon]},While[True,\[CurlyEpsilon] = RandomReal[{1*10^-9, 10^-5}];If[\[CurlyEpsilon]>\[CurlyEpsilon]q && \[CurlyEpsilon]q<(1-\[CurlyEpsilon])(\[CurlyEpsilon]-\[CurlyEpsilon]q), Return[\[CurlyEpsilon]]]]]
\[Gamma]List=RandomReal[{10^-3,0.05},sampleSize];*)


GenerateSamplesNS[SampleSize_,SamplingType_]:=Module[{I3List,Q22List,\[Gamma]List,generateEpsilon},If[!MemberQ[{"Uniform","LogUniform"},SamplingType],Message[GenerateSamples::invtype,"\:65e0\:6548\:7684\:91c7\:6837\:7c7b\:578b\:ff0c\:5fc5\:987b\:4e3a\"Uniform\"\:6216\"LogUniform\""];
Return[$Failed]];
Which[SamplingType=="Uniform",(I3List=RandomReal[{10^38,3*10^38},SampleSize];
Q22List=RandomReal[{4.4*10^28,2.4*10^32},SampleSize];
\[Gamma]List=RandomReal[{10^-3,0.05},SampleSize];
generateEpsilon[\[CurlyEpsilon]q_]:=Module[{\[CurlyEpsilon]},While[True,\[CurlyEpsilon]=RandomReal[{1*10^-9,10^-5}];
If[\[CurlyEpsilon]>\[CurlyEpsilon]q&&\[CurlyEpsilon]q<(1-\[CurlyEpsilon])(\[CurlyEpsilon]-\[CurlyEpsilon]q),Return[\[CurlyEpsilon]]]]]),SamplingType=="LogUniform",(I3List=10^RandomReal[{38,38+Log10[3]},SampleSize];Q22List=10^RandomReal[{Log10[4.4*10^28],Log10[2.4*10^32]},SampleSize];\[Gamma]List=10^RandomReal[{-3,Log10[0.05]},SampleSize];generateEpsilon[\[CurlyEpsilon]q_]:=Module[{\[CurlyEpsilon]},While[True,\[CurlyEpsilon]=10^RandomReal[{-9,-5}];If[\[CurlyEpsilon]>\[CurlyEpsilon]q&&\[CurlyEpsilon]q<(1-\[CurlyEpsilon])(\[CurlyEpsilon]-\[CurlyEpsilon]q),Return[\[CurlyEpsilon]]]]])];
{I3List,Q22List,\[Gamma]List,generateEpsilon}]


{I3List,Q22List,\[Gamma]List,generateEpsilon}=GenerateSamplesNS[sampleSize,SamplingType];


\[Epsilon]qList=Sqrt[8\[Pi]/15]Q22List/I3List;
\[Epsilon]List = Table[generateEpsilon[\[Epsilon]qList[[i]]], {i, 1, sampleSize}];
\[Kappa]List=\[Epsilon]qList/(16 (1-\[Epsilon]List) (\[Epsilon]List-\[Epsilon]qList));

\[Theta]SList=ArcCos[RandomReal[{-1,1},sampleSize]]; 
bList=bfun[PsList,\[Gamma]List];
\[CapitalOmega]rList=MapThread[\[CapitalOmega]rfun,{\[Epsilon]List,\[Gamma]List,\[Kappa]List,PsList}];
\[CapitalOmega]pList=MapThread[\[CapitalOmega]pfun,{\[Epsilon]List,\[Gamma]List,\[Kappa]List,PsList}];


snrLISAList=Import[basePath <> "snr_LISA_alpha10.txt","List"];
snrLISAkde=SmoothKernelDistribution[snrLISAList];
(*SNRLISAList=RandomVariate[snrLISAkde,sampleSize];*)


ASDofCEdata=Import[basePath <> "cosmic_explorer_40km_lf.txt","Table"];
fsList=1/PsList;
fs2List=2/PsList;
(*\:5b9a\:4e49\:4e00\:4e2a\:51fd\:6570\:ff0c\:7528\:4e8e\:627e\:51fa\:6700\:63a5\:8fd1\:7684 ASD \:503c*)
findClosestASD[f_]:=Module[{differences,index},differences=Abs[ASDofCEdata[[All,1]]-f];
index=Position[differences,Min[differences]][[1,1]];
ASDofCEdata[[index,2]]];
ASDofCEat1fs[Ps_]:=findClosestASD/@(1/Ps)/Sqrt[2]
ASDofCEat2fs[Ps_]:=findClosestASD/@(2/Ps)/Sqrt[2]


h1inBinaryDetectedList[t_]:=MapThread[h1inBinaryDetected[t,#1,#2,#3,#4,#5,#6,#7,#8,#9,#10,#11,#12,#13,#14,#15,#16,#17]&,{PsList,I3List,\[Epsilon]List,\[Kappa]List,\[Gamma]List,rList,\[Alpha]List,\[Delta]dList,\[Psi]List,\[Iota]List,\[CapitalOmega]rList,\[CapitalOmega]pList,PorbList,\[Theta]SList,eList,m1List,m2List}]/ASDofCEat1fs[PsList]
h2inBinaryDetectedList[t_]:=MapThread[h2inBinaryDetected[t,#1,#2,#3,#4,#5,#6,#7,#8,#9,#10,#11,#12,#13,#14,#15,#16,#17]&,{PsList,I3List,\[Epsilon]List,\[Kappa]List,\[Gamma]List,rList,\[Alpha]List,\[Delta]dList,\[Psi]List,\[Iota]List,\[CapitalOmega]rList,\[CapitalOmega]pList,PorbList,\[Theta]SList,eList,m1List,m2List}]/ASDofCEat2fs[PsList]
h3inBinaryDetectedList[t_]:=MapThread[h3inBinaryDetected[t,#1,#2,#3,#4,#5,#6,#7,#8,#9,#10,#11,#12,#13,#14,#15,#16,#17]&,{PsList,I3List,\[Epsilon]List,\[Kappa]List,\[Gamma]List,rList,\[Alpha]List,\[Delta]dList,\[Psi]List,\[Iota]List,\[CapitalOmega]rList,\[CapitalOmega]pList,PorbList,\[Theta]SList,eList,m1List,m2List}]/ASDofCEat2fs[PsList]


h4inBinaryDetectedList[t_]:=MapThread[h4inBinaryDetected[t,#1,#2,#3,#4,#5,#6,#7,#8,#9,#10,#11,#12,#13,#14,#15,#16,#17]&,{PsList,I3List,\[Epsilon]List,\[Kappa]List,\[Gamma]List,rList,\[Alpha]List,\[Delta]dList,\[Psi]List,\[Iota]List,\[CapitalOmega]rList,\[CapitalOmega]pList,PorbList,\[Theta]SList,eList,m1List,m2List}]/ASDofCEat1fs[PsList]
h5inBinaryDetectedList[t_]:=MapThread[h5inBinaryDetected[t,#1,#2,#3,#4,#5,#6,#7,#8,#9,#10,#11,#12,#13,#14,#15,#16,#17]&,{PsList,I3List,\[Epsilon]List,\[Kappa]List,\[Gamma]List,rList,\[Alpha]List,\[Delta]dList,\[Psi]List,\[Iota]List,\[CapitalOmega]rList,\[CapitalOmega]pList,PorbList,\[Theta]SList,eList,m1List,m2List}]/ASDofCEat1fs[PsList]
h6inBinaryDetectedList[t_]:=MapThread[h6inBinaryDetected[t,#1,#2,#3,#4,#5,#6,#7,#8,#9,#10,#11,#12,#13,#14,#15,#16,#17]&,{PsList,I3List,\[Epsilon]List,\[Kappa]List,\[Gamma]List,rList,\[Alpha]List,\[Delta]dList,\[Psi]List,\[Iota]List,\[CapitalOmega]rList,\[CapitalOmega]pList,PorbList,\[Theta]SList,eList,m1List,m2List}]/ASDofCEat2fs[PsList]


(*signal-to-noise ratio function for parallel computation*)


parallelCalculateSNR[type_, td_] := ParallelMap[
  (Quiet[
      NIntegrate[#*Conjugate[#], {t, 0, td},
        WorkingPrecision -> 20,
        Method -> "DoubleExponential"
      ],
      {NIntegrate::precw, NIntegrate::ncvi, General::stop}
    ]
  )^0.5 &,
  ToExpression["h" <> ToString[type] <> "inBinaryDetectedList"][t]
]


YearInS=31536000.0;
Tobs=4.0*YearInS;
(*Off[General::stop,NIntegrate::precw,NIntegrate::ncvi,Inverse::luc]*)
snrH1aList=parallelCalculateSNR[1,Tobs];
(*Print["snrH1aList: ",snrH1aList]*)
(*snrH1bList=parallelCalculateSNR[4,Tobs];*)
(*Print["snrH1bList: ",snrH1bList]*)
(*snrH1cList=parallelCalculateSNR[5,Tobs];*)
(*Print["snrH1cList: ",snrH1cList]*)
snrH2aList=parallelCalculateSNR[2,Tobs];
(*Print["snrH2aList: ",snrH2aList]*)
snrH2bList=parallelCalculateSNR[3,Tobs];
(*Print["snrH2bList: ",snrH2bList]*)
snrH2cList=parallelCalculateSNR[6,Tobs];
(*Print["snrH2cList: ",snrH2cList]*)
(*Print["This single simulation has been done!"]*)


indices=Position[Transpose[{snrH2aList,snrH2bList,snrH2cList}],{x_,y_,z_}/;x>7&&y>7&&z>7];
Selectsnrh2a=Extract[snrH2aList,indices];
Selectsnrh2c=Extract[snrH2cList,indices];
Select\[Kappa]=Extract[\[Kappa]List,indices];
SNRLISAList=RandomVariate[snrLISAkde,Length[indices]];
\[CapitalDelta]I3toI3=Sqrt[1/Selectsnrh2a^2+(1-32 Select\[Kappa])^2 (1/Selectsnrh2a^2+1/Selectsnrh2c^2)+1/SNRLISAList^2];
