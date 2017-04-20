#### Analytical model of catchment driven lake carbon model
#### 04-18-2017
#### S.E. Jones, University of Notre Dame

rm(list=ls())

#### arguments
	# A: lake area [m2]
	# zbar:  lake mean depth [m]
	# WA2LA: watershed area to lake area ratio (drainage ratio)
	# MAP:  mean annual precipitation [m yr-1]
	# ET: mean daily evapotranspiration rate of catchment [m day-1]
	# evap:  mean daily evaporative loss from lake [m day-1]
	# Cq: total DOC concentration in inlet [g m-3]
	# Cp: total DOC concentration in precip [g m-3]
	# d1:  recalcitrant DOC processing rate [day-1]
	# d2:  labile DOC process rate [day-1]
	# fLabile:  fraction of total DOC that is labile in fluvial load
			
rateVfate<-function(A=10000,zbar=10,WA2LA=5,MAP=0.8,ET=0.0014,evap=0.0025,Cq=5,Cp=1,d1=0.005,d2=0.02,fLabile=0){
	V=zbar*A		# lake volume [m3]
	WA=A*WA2LA	# watershed area [m2]

	# hydrology
	precip=MAP/365			# daily precipitation [m day-1]
	Qin=(precip-ET)*WA		# inflow to lake [m3 day-1]
	if(Qin<0){Qin=0}		
	directP=precip*A		#total direct precip [m3 day-1]
	Etot=evap*A				#total evap [m3 day-1]
	Qout=Qin+directP-Etot	#outflow from lake [m3 day-1]

	if(Qout<0){print("Warning: A negative Qout was used")}				# doesn't allow model to run for situation where evaporation would drain the lake

	fracEvap=Etot/(Etot+Qout)	# fraction of hydrologic losses due to evaporation
	RT=V/(Qin+directP)		# residence time [days]
	
	C1q=Cq*(1-fLabile)	# concentration of recalcitrant DOC in fluvial inputs
	C2q=Cq*fLabile		# concentration of labile DOC in fluvial inputs
	C1p=Cp				# precipitation carries only recalcitrant DOC
	C2p=0

	# carbon
	C1star=(Qin*C1q+directP*C1p)/(Qout/V+d1)	# mass of recalcitrant DOC at equilibrium [g C]
	C2star=(Qin*C2q+directP*C2p)/(Qout/V+d2)	# mass of labile DOC at equilibrium [g C]
	Cstar=C1star+C2star						# mass of total DOC at equilibrium [g C]
	C1conc=C1star/V					# recalcitrant DOC concentration at equilibrium [g C m-3]
	C2conc=C2star/V					# labile DOC concentration at equilibrium [g C m-3]
	Cconc=C1conc+C2conc				# total DOC concentration at equilibrium [g C m-3]
	Cload=Qin*Cq+directP*Cp			# daily total DOC load [g day-1]
	C1load=Qin*C1q+directP*C1p			# daily recalcitrant DOC load [g day-1]
	C2load=Qin*C2q+directP*C2p			# daily labile DOC load [g day-1]
	fracCresp=1-Qout*Cconc/Cload		# fraction of DOC load that is decomposed (not exported)
	resp=C1star*d1+C2star*d2			# total daily respiration at equilibrium [g C day-1]
	respVol=resp/V				# volumetric daily respiration at equilibrium [g C m-3 day-1]
	
	return(as.list(c(A=A,zbar=zbar,WA2LA=WA2LA,MAP=MAP,ET=ET,evap=evap,Cq=Cq,Cp=Cp,d1=d1,d2=d2,V=V,WA=WA,Qin=Qin,directP=directP,Etot=Etot,Qout=Qout,fracEvap=fracEvap,RT=RT,Cstar=Cstar,Cconc=Cconc,Cload=Cload,fracCresp=fracCresp,resp=resp,respVol=respVol)))
}








#### (1) simulate over ranges of WA2LAs and zbars
WA2LAs=seq(0.5,1000,0.5)	# drainage ratios considered
zbars=1:15					# mean depths considered [m]

# matrices for storing equilibrium values
RT=matrix(NA,length(WA2LAs),length(zbars))
fracEvap=RT
fracCresp=RT
respiration=RT
respirationVol=RT
Cload=RT
fracCresp=RT
Qout=RT
equilC=RT

# loop over gradients in factorial manner
for(i in 1:length(WA2LAs)){
	for(j in 1:length(zbars)){
		temp=rateVfate(zbar=zbars[j],WA2LA=WA2LAs[i])	# solve for equilbria
		equilC[i,j]=temp$Cconc				# equilibrium DOC concentration [g m-3]
		fracEvap[i,j]=temp$fracEvap		# fraction of hydrologic loss that is evap 
		respiration[i,j]=temp$resp			# daily respiration [g day-1]
		respirationVol[i,j]=temp$respVol	# daily volumetric respiration [g m-3 day-1]
		Cload[i,j]=temp$Cload				# daily carbon load [g day-1]
		fracCresp[i,j]=temp$fracCresp		# fraction of loaded C that is respired
		RT[i,j]=temp$RT					# residence time [day]
		Qout[i,j]=temp$Qout				# outflow [m3 day-1]
	}
}

#### (2) simulate for same gradients as (1) above, but without any evap from lake
# matrices for storying equilibrium values
RTNE=matrix(NA,length(WA2LAs),length(zbars))
fracEvapNE=RT
fracCrespNE=RT
respirationNE=RT
respirationVolNE=RT
CloadNE=RT
fracCrespNE=RT
QoutNE=RT
equilCNE=RT

# loop over gradients in factorial manner
for(i in 1:length(WA2LAs)){
	for(j in 1:length(zbars)){
		temp=rateVfate(zbar=zbars[j],WA2LA=WA2LAs[i],evap=0)	
		equilCNE[i,j]=temp$Cconc				# equilibrium DOC concentration [g m-3]
		fracEvapNE[i,j]=temp$fracEvap		# fraction of hydrologic loss that is evap 
		respirationNE[i,j]=temp$resp			# daily respiration [g day-1]
		respirationVolNE[i,j]=temp$respVol	# daily volumetric respiration [g m-3 day-1]
		CloadNE[i,j]=temp$Cload				# daily carbon load [g day-1]
		fracCrespNE[i,j]=temp$fracCresp		# fraction of loaded C that is respired
		RTNE[i,j]=temp$RT					# residence time [day]
		QoutNE[i,j]=temp$Qout				# outflow [m3 day-1]
	}
}

#### (3) simulate over ranges of WA2LAs, zbars and inflow DOC concentrations
DOCins=seq(5,40,5)		# fluvial total DOC concentrations considered [g C m-3]

# 3-dimensional array to store equilibrium values
RT3=array(NA,dim=c(length(WA2LAs),length(DOCins),length(zbars)))
fracEvap3=RT3
fracCresp3=RT3
respiration3=RT3
respirationVol3=RT3
Cload3=RT3
fracCresp3=RT3
Qout3=RT3
equilC3=RT3

# loop over gradients in factorial manner
for(i in 1:length(WA2LAs)){
	for(j in 1:length(DOCins)){
		for(k in 1:length(zbars)){
			temp=rateVfate(zbar=zbars[k],WA2LA=WA2LAs[i],Cq=DOCins[j])	
			equilC3[i,j,k]=temp$Cconc				# equilibrium DOC concentration [g m-3]
			fracEvap3[i,j,k]=temp$fracEvap		# fraction of hydrologic loss that is evap 
			respiration3[i,j,k]=temp$resp			# daily respiration [g day-1]
			respirationVol3[i,j,k]=temp$respVol	# daily volumetric respiration [g m-3 day-1]
			Cload3[i,j,k]=temp$Cload				# daily carbon load [g day-1]
			fracCresp3[i,j,k]=temp$fracCresp		# fraction of loaded C that is respired
			RT3[i,j,k]=temp$RT					# residence time [day]
			Qout3[i,j,k]=temp$Qout				# outflow [m3 day-1]
		}
	}
}

#### (4) simulate over ranges of WA2LAs, zbars, and inflow DOC fraction labile
fLs=seq(0.5,0,-0.01)		# fractions Labile considered

# 3-dimensional arrays to store output from simulations along gradients of drainage ratio, mean depth, and fraction labile
RT3q=array(NA,dim=c(length(WA2LAs),length(fLs),length(zbars)))
fracEvap3q=RT3q
fracCresp3q=RT3q
respiration3q=RT3q
respirationVol3q=RT3q
Cload3q=RT3q
fracCresp3q=RT3q
Qout3q=RT3q
equilC3q=RT3q

# loops to simulate across all gradients in a factorial manner
for(i in 1:length(WA2LAs)){
	for(j in 1:length(fLs)){
		for(k in 1:length(zbars)){
			temp=rateVfate(zbar=zbars[k],WA2LA=WA2LAs[i],fLabile=fLs[j])	
			equilC3q[i,j,k]=temp$Cconc				# equilibrium DOC concentration [g m-3]
			fracEvap3q[i,j,k]=temp$fracEvap		# fraction of hydrologic loss that is evap 
			respiration3q[i,j,k]=temp$resp			# daily respiration [g day-1]
			respirationVol3q[i,j,k]=temp$respVol	# daily volumetric respiration [g m-3 day-1]
			Cload3q[i,j,k]=temp$Cload				# daily carbon load [g day-1]
			fracCresp3q[i,j,k]=temp$fracCresp		# fraction of loaded C that is respired
			RT3q[i,j,k]=temp$RT					# residence time [day]
			Qout3q[i,j,k]=temp$Qout				# outflow [m3 day-1]
		}
	}
}







#*********************
# Fig. 1 - proof of concept
#*********************

# a) RT - log v. log
dev.new()
plot(WA2LAs,RT[,1],type='l',xlab="log10(WA:LA)",ylab="log10(residence time) (days)",ylim=c(1,max(RT)),log="xy",lwd=2,col='grey90')
lines(WA2LAs,RT[,2],lwd=2,col='grey70')
lines(WA2LAs,RT[,5],lwd=2,col='grey55')
lines(WA2LAs,RT[,10],lwd=2,col='grey40')
lines(WA2LAs,RT[,15],lwd=2)
legend('bottomleft',legend=c(1,2,5,10,15),lty=1,col=c('grey90','grey70','grey55','grey40','black'),box.lty=0)
text(1,10,'mean depth (m)',font=4)
# b) areal C load
dev.new()
plot(WA2LAs,Cload[,1]/10000,type='l',xlab="WA:LA",ylab="carbon load (g m-2 day-1)",ylim=c(0,max(Cload/10000)))
# c) equilibrium C concentration
dev.new()
plot(WA2LAs,equilC[,1],type='l',xlab="WA:LA",ylab="equilibrium carbon concentration (g m-3)",ylim=c(0,max(equilC)),lwd=2,col='grey95')
lines(WA2LAs,equilC[,2],lwd=2,col='grey70')
lines(WA2LAs,equilC[,5],lwd=2,col='grey55')
lines(WA2LAs,equilC[,10],lwd=2,col='grey40')
lines(WA2LAs,equilC[,15],lwd=2)
legend('bottomright',legend=c(1,2,5,10,15),lty=1,col=c('grey90','grey70','grey55','grey40','black'),box.lty=0)
text(900,1.4,'mean depth (m)',font=4)

#extra plots, a la Brett et al. 2012 fig 2a & b
dev.new()
plot(log10(WA2LAs*7.917808*1000/10000),log10(1/RT[,1]),xlab="Areal hydraulic load (L m-2 d-1)",ylab="flushign rate (day-1)",xlim=c(0,3),type='l',ylim=c(-4,0),lwd=2,col='grey90')
lines(log10(WA2LAs*7.917808*1000/10000),log10(1/RT[,2]),lwd=2,col='grey70')
lines(log10(WA2LAs*7.917808*1000/10000),log10(1/RT[,5]),lwd=2,col='grey55')
lines(log10(WA2LAs*7.917808*1000/10000),log10(1/RT[,10]),lwd=2,col='grey40')
lines(log10(WA2LAs*7.917808*1000/10000),log10(1/RT[,15]))
lines(log10(seq(1,1000,10)),log10(0.0001*seq(1,1000,10)^1.084),col='red')	# model fit from Brett et al.
legend('topleft',legend=c(1,2,5,10,15,"Brett et al."),lty=1,col=c('grey90','grey70','grey55','grey40','black','red'),box.lty=0)

dev.new()
plot(log10(WA2LAs*7.917808*1000/10000),log10(Cload[,1]/10000*1000),type='l',xlab="Areal hydraulic load (L m-2 d-1)",ylab="carbon load (mg m-2 day-1)",xlim=c(0,3),lwd=2)
lines(log10(seq(1,1000,10)),log10(8.37*seq(1,1000,10)^0.982),lwd=2,col='red')
legend('topleft',c('this study','Brett et al. 2012'),lty=1,col=c('black','red'),box.lty=0)


#*********************
# Fig 2 - evaporation matters
#*********************

# a) frac evap
dev.new()
plot(WA2LAs,fracEvap[,1],type='l',xlab="WA:LA",ylab="fraction evap",ylim=c(0,1),log="x",lwd=2)		# all lines are on top of eachother

# b) Qout higher in low WA:LA lakes where evap matters when present
dev.new()
plot(Qout[,1],QoutNE[,1],type='l',log="xy",xlim=c(1,10000),ylim=c(1,10000),xlab="Qout - With evaporation",ylab="Qout - No evaporation",lwd=2)
arrows(x0=Qout[(nrow(Qout)-4),1],y0=QoutNE[(nrow(QoutNE)-4),1],x1=Qout[nrow(Qout),1],y1=QoutNE[nrow(QoutNE),1],lwd=2)
abline(a=0,b=1,lty=5)

# c) equil C conc - when evap is important in a low WA2LA lake you get a "concentrating" effect, but matters across WA2LA gradient (in a shallow lake; in a big lake water losses to evap are relatively small and don't have a concentrating effect; sort of a product of constant volume...)
dev.new()
plot(equilC[,1],equilCNE[,1],type='l',xlim=c(0,5),ylim=c(0,5),xlab="Equilibrium C - With evaporation",ylab="Equilibrium C - No evaporation",lwd=2,col='grey')
arrows(x0=equilC[(nrow(equilC)-50),1],y0=equilCNE[(nrow(equilCNE)-50),1],x1=equilC[nrow(equilC),1],y1=equilCNE[nrow(equilCNE),1],lwd=2,col='grey')
lines(equilC[,15],equilCNE[,15],lwd=2)
arrows(x0=equilC[(nrow(equilC)-50),15],y0=equilCNE[(nrow(equilCNE)-50),15],x1=equilC[nrow(equilC),15],y1=equilCNE[nrow(equilCNE),15],lwd=2)
abline(a=0,b=1,lty=5)

# d) fraction respiration - effect of evaporation
dev.new()
plot(fracCresp[,1],fracCrespNE[,1],type='l',xlim=c(0,1),ylim=c(0,1),xlab="fraction decomposed - With evaporation",ylab="fraction decomposed - No evaporation",lwd=2,col='grey')
arrows(x0=fracCresp[(nrow(fracCresp)-50),1],y0=fracCrespNE[(nrow(fracCrespNE)-50),1],x1=fracCresp[nrow(fracCresp),1],y1=fracCrespNE[nrow(fracCrespNE),1],lwd=2,col='grey')
lines(fracCresp[,15],fracCrespNE[,15],lwd=2)
arrows(x0=fracCresp[(nrow(fracCresp)-50),15],y0=fracCrespNE[(nrow(fracCrespNE)-50),15],x1=fracCresp[nrow(fracCresp),15],y1=fracCrespNE[nrow(fracCrespNE),15],lwd=2)
abline(a=0,b=1,lty=5)

#*********************
# Fig 3 - rate v. fate
#*********************
# a) volumetric respiration 
dev.new()
plot(WA2LAs,respirationVol[,1],type='l',xlab="WA:LA",ylab="volumetric respiration (mg C m-3 day-1)",ylim=c(0,max(respirationVol)),lwd=2,col='grey90')
lines(WA2LAs,respirationVol[,2],col='grey70',lwd=2)
lines(WA2LAs,respirationVol[,5],col='grey55',lwd=2)
lines(WA2LAs,respirationVol[,10],col='grey40',lwd=2)
lines(WA2LAs,respirationVol[,15],lwd=2)
legend('bottomright',legend=c(1,2,5,10,15),lty=1,col=c('grey90','grey70','grey55','grey40','black'),box.lty=0)
text(900,0.006,'mean depth (m)',font=4)

# b) frac C respired
dev.new()
plot(WA2LAs,fracCresp[,1],type='l',xlab="WA:LA",ylab="fraction C load respired",ylim=c(0,1),lwd=2,col='grey90')
lines(WA2LAs,fracCresp[,2],col='grey70',lwd=2)
lines(WA2LAs,fracCresp[,5],col='grey55',lwd=2)
lines(WA2LAs,fracCresp[,10],col='grey40',lwd=2)
lines(WA2LAs,fracCresp[,15],lwd=2)
legend('topright',legend=c(1,2,5,10,15),lty=1,col=c('grey90','grey70','grey55','grey40','black'),box.lty=0)

# c) volumetric respiration vs. fraction respired
dev.new()
plot(respirationVol[,1],fracCresp[,1],xlab="Volumetric respiration (g C m-3 day-1)",ylab="fraction C load respired",xlim=c(0,0.025),ylim=c(0,1),type='l',lwd=2,col='grey90')	#zbar=1
lines(respirationVol[,2],fracCresp[,2],lwd=2,col='grey70')	#zbar=2
lines(respirationVol[,5],fracCresp[,5],lwd=2,col='grey55')	#zbar=5
lines(respirationVol[,10],fracCresp[,10],lwd=2,col='grey40')	#zbar=10
lines(respirationVol[,15],fracCresp[,15],lwd=2)	#zbar=15
legend('bottomleft',legend=c(1,2,5,10,15),lty=1,col=c('grey90','grey70','grey55','grey40','black'),box.lty=0)
text(0.002,0.26,'mean depth (m)',font=4)


#*********************
# Fig 4 - multi-line graph - RT v. frac C resp
#*********************
dev.new()
plot(RT[1,],fracCresp[1,],type='l',xlim=c(1,10000),ylim=c(0,1),xlab="residence time (days)",ylab="fraction C load decomposed",lwd=2,log="x",col='black')	#WA2LA=0.5
lines(RT[2,],fracCresp[2,],lwd=2,col='purple')	#WA2LA=1
lines(RT[4,],fracCresp[4,],lwd=2,col='blue')	#WA2LA=2
lines(RT[10,],fracCresp[10,],lwd=2,col='darkgreen')	#WA2LA=5
lines(RT[20,],fracCresp[20,],lwd=2,col='green')	#WA2LA=10
lines(RT[100,],fracCresp[100,],lwd=2,col='yellow')	#WA2LA=50
lines(RT[200,],fracCresp[200,],lwd=2,col='orange')	#WA2LA=100
lines(RT[1000,],fracCresp[1000,],lwd=2,col='pink')	#WA2LA=500
lines(RT[2000,],fracCresp[2000,],lwd=2,col='red')	#WA2LA=1000

legend('topleft',c('DR=0.5','DR=1','DR=2','DR=5','DR=10','DR=50','DR=100','DR=500','DR=1000'),lty=1,col=c('black','purple','blue','darkgreen','green','yellow','orange','pink','red'),box.lty=0)


###############
# Figure 5 - quantity and quality effects
###############
# a & b) effect of Cq (DOC concentratio in inflow) and C quality (fraction recalcitrant)
dev.new()
par(mfrow=c(1,2))
plot(DOCins,respirationVol3[1,,1],type='l',xlab="Inlet [DOC] (g C m-3)",ylab="volumetric respiration (g C m-3 day-1)",ylim=c(0,max(respirationVol3)),col='black',main="mean depth=1 m")
lines(DOCins,respirationVol3[2,,1],col='purple')
lines(DOCins,respirationVol3[4,,1],col='blue')
lines(DOCins,respirationVol3[10,,1],col='darkgreen')
lines(DOCins,respirationVol3[20,,1],col='green')
lines(DOCins,respirationVol3[100,,1],col='yellow')
lines(DOCins,respirationVol3[200,,1],col='orange')
lines(DOCins,respirationVol3[1000,,1],col='pink')
lines(DOCins,respirationVol3[2000,,1],col='red')

plot(DOCins,respirationVol3[1,,15],type='l',xlab="Inlet [DOC] (g C m-3)",ylab="volumetric respiration (g C m-3 day-1)",ylim=c(0,max(respirationVol3)),col='black',main="mean depth=15 m")
lines(DOCins,respirationVol3[2,,15],col='purple')
lines(DOCins,respirationVol3[4,,15],col='blue')
lines(DOCins,respirationVol3[10,,15],col='darkgreen')
lines(DOCins,respirationVol3[20,,15],col='green')
lines(DOCins,respirationVol3[100,,15],col='yellow')
lines(DOCins,respirationVol3[200,,15],col='orange')
lines(DOCins,respirationVol3[1000,,15],col='pink')
lines(DOCins,respirationVol3[2000,,15],col='red')

# c & d) quality effect on volumetric respiration
dev.new()
par(mfrow=c(1,2))
plot(fLs,respirationVol3q[1,,1],type='l',xlab="fraction of inlet C that is labile",ylab="volumetric respiration (g C m-3 day-1)",ylim=c(0,max(respirationVol3q)),col='black',main="mean depth=1 m")
lines(fLs,respirationVol3q[2,,1],col='purple')
lines(fLs,respirationVol3q[4,,1],col='blue')
lines(fLs,respirationVol3q[10,,1],col='darkgreen')
lines(fLs,respirationVol3q[20,,1],col='green')
lines(fLs,respirationVol3q[100,,1],col='yellow')
lines(fLs,respirationVol3q[200,,1],col='orange')
lines(fLs,respirationVol3q[1000,,1],col='pink')
lines(fLs,respirationVol3q[2000,,1],col='red')

plot(fLs,respirationVol3q[1,,15],type='l',xlab="fraction of inlet C that is labile",ylab="volumetric respiration (g C m-3 day-1)",ylim=c(0,max(respirationVol3q)),col='black',main="mean depth=15 m")
lines(fLs,respirationVol3q[2,,15],col='purple')
lines(fLs,respirationVol3q[4,,15],col='blue')
lines(fLs,respirationVol3q[10,,15],col='darkgreen')
lines(fLs,respirationVol3q[20,,15],col='green')
lines(fLs,respirationVol3q[100,,15],col='yellow')
lines(fLs,respirationVol3q[200,,15],col='orange')
lines(fLs,respirationVol3q[1000,,15],col='pink')
lines(fLs,respirationVol3q[2000,,15],col='red')

# e & f) quality effect on fraction DOC decomposed
dev.new()
par(mfrow=c(1,2))
plot(fLs,fracCresp3q[1,,1],type='l',xlab="fraction of inlet C that is labile",ylab="fracion DOC load decomposed",ylim=c(0,1),col='black',main="mean depth=1 m")
lines(fLs,fracCresp3q[2,,1],col='purple')
lines(fLs,fracCresp3q[4,,1],col='blue')
lines(fLs,fracCresp3q[10,,1],col='darkgreen')
lines(fLs,fracCresp3q[20,,1],col='green')
lines(fLs,fracCresp3q[100,,1],col='yellow')
lines(fLs,fracCresp3q[200,,1],col='orange')
lines(fLs,fracCresp3q[1000,,1],col='pink')
lines(fLs,fracCresp3q[2000,,1],col='red')

plot(fLs,fracCresp3q[1,,15],type='l',xlab="fraction of inlet C that is labile",ylab="fracion DOC load decomposed",ylim=c(0,1),col='black',main="mean depth=15 m")
lines(fLs,fracCresp3q[2,,15],col='purple')
lines(fLs,fracCresp3q[4,,15],col='blue')
lines(fLs,fracCresp3q[10,,15],col='darkgreen')
lines(fLs,fracCresp3q[20,,15],col='green')
lines(fLs,fracCresp3q[100,,15],col='yellow')
lines(fLs,fracCresp3q[200,,15],col='orange')
lines(fLs,fracCresp3q[1000,,15],col='pink')
lines(fLs,fracCresp3q[2000,,15],col='red')

# g & h) quality effect on emergent d - volumetric R/equilC
d3q=respirationVol3q/equilC3q
dev.new()
par(mfrow=c(1,2))
plot(fLs,d3q[1,,1],type='l',xlab="fraction of inlet C that is labile",ylab="emergent d (day-1)",ylim=c(0.005,max(d3q)),col='black',main="mean depth=1 m")
lines(fLs,d3q[2,,1],col='purple')
lines(fLs,d3q[4,,1],col='blue')
lines(fLs,d3q[10,,1],col='darkgreen')
lines(fLs,d3q[20,,1],col='green')
lines(fLs,d3q[100,,1],col='yellow')
lines(fLs,d3q[200,,1],col='orange')
lines(fLs,d3q[1000,,1],col='pink')
lines(fLs,d3q[2000,,1],col='red')

plot(fLs,d3q[1,,15],type='l',xlab="fraction of inlet C that is labile",ylab="emergent d (day-1)",ylim=c(0.005,max(d3q)),col='black',main="mean depth=15 m")
lines(fLs,d3q[2,,15],col='purple')
lines(fLs,d3q[4,,15],col='blue')
lines(fLs,d3q[10,,15],col='darkgreen')
lines(fLs,d3q[20,,15],col='green')
lines(fLs,d3q[100,,15],col='yellow')
lines(fLs,d3q[200,,15],col='orange')
lines(fLs,d3q[1000,,15],col='pink')
lines(fLs,d3q[2000,,15],col='red')