'___________________________________________________________________________________________________________________________________________________________________________

'  ESTIMATE NUETRAL INTEREST RATE FOR AUSTRALIA
' email: david.stephan@gmail.com ;
' reference: https://www.rba.gov.au/publications/bulletin/2017/sep/pdf/bu-0917-2-the-neutral-interest-rate.pdf

'___________________________________________________________________________________________________________________________________________________________________________
close @all
%path = @runpath
cd %path

logmode l

!import_data = 0 'Set to 1 if you want to load the rawdata from ABS/RBA/FRED. Set to zero if you want it to load the previously imported database
include MUE_Stage1.prg 'Stage 1 Median Unbiased Estimator
include MUE_Stage2.prg 'Stage 2 Median Unbiased Estimator

'-----------------------------------------------------------------------------
if !import_data = 1 then

	logmsg "Importing Data"

	'Import Rawdata from ABS/RBA

	'Create Workfile to populate with ABS/RBA Data
	wfcreate(wf=RSTAR, page=RSTAR) q 1959q3 2021q4

	'Set Strings to Load Data using R packages readabs and readrba
	%rpath = @replace(@runpath,"\","/")
	%setwd = "setwd("+"""" + @left(%rpath,@len(%rpath)-1) + """)"

	'Open R connection
	xopen(r)

	'Set wd so can access and store R dataframes
	xrun {%setwd}

	'Download RBA and ABS Data
	xpackage readabs
	xpackage readrba
	xpackage readr
	xpackage reshape2
	xpackage lubridate
	xpackage tidyverse
	xpackage zoo
	xpackage data.table

	xon

		'---------------------------------------------------------------------------------------------------------
		'Download Most Recent ABS Data
		'---------------------------------------------------------------------------------------------------------
		'Import Data from ABS Website
		abs_5206 = read_abs(series_id = c("A2304402X", "A2302915V"))
		abs_6202 = read_abs(series_id = c("A84423043C", "A84423047L"))

		'5206.0 Australian National Accounts: National Income, Expenditure and Product
		R_5206 = abs_5206 %>%  filter(series_id %in% c("A2304402X")) %>% mutate(date = zoo::as.yearqtr(date)) %>% dplyr::select(date, series_id, value) 
		R_5206 = distinct(R_5206,date,series_id, .keep_all= TRUE)
		R_5206 = dcast(R_5206, date ~ series_id)

		'6202.0 Labour Force, Australia - Monthly
		R_6202 = abs_6202 %>% filter(series_id %in% c("A84423043C", "A84423047L")) %>% dplyr::select(date, series_id, value)
		R_6202 = distinct(R_6202,date,series_id, .keep_all= TRUE)
		R_6202 = dcast(R_6202, date ~ series_id)
		R_6202 = R_6202 %>% group_by(date=floor_date(date, "quarter")) %>% summarize(A84423043C=mean(A84423043C), A84423047L=mean(A84423047L)) %>% mutate(date = zoo::as.yearqtr(date))
		
		'---------------------------------------------------------------------------------------------------------
		'Download Most Recent RBA Data
		'---------------------------------------------------------------------------------------------------------
			'Trimmed-Mean Inflation (Year-Ended)
			rba_g1 = read_rba(series_id = c("GCPIOCPMTMYP")) 
			R_g1 = rba_g1 %>% filter(series_id %in% c("GCPIOCPMTMYP")) %>% mutate(date = zoo::as.yearqtr(date)) %>% dplyr::select(date, series_id, value)
			R_g1 = distinct(R_g1,date,series_id, .keep_all= TRUE)
			R_g1 = dcast(R_g1, date ~ series_id)
			
			'Bond-market inflation expectations
			rba_g3 = read_rba(series_id = c("GBONYLD")) 
			R_g3 = rba_g3 %>% filter(series_id %in% c("GBONYLD")) %>% mutate(date = zoo::as.yearqtr(date)) %>% dplyr::select(date, series_id, value)
			R_g3 = distinct(R_g3,date,series_id, .keep_all= TRUE)
			R_g3 = dcast(R_g3, date ~ series_id)

			'Interbank Overnight Cash Rate; monthly average
			rba_f1 = read_rba(series_id = c("FIRMMCRI")) 
			R_f1 = rba_f1 %>% arrange(date) %>% group_by(date=floor_date(date, "quarter")) %>% summarize(FIRMMCRI=mean(value)) %>% mutate(date = zoo::as.yearqtr(date))

		'---------------------------------------------------------------------------------------------------------
		'Download COVID Stringency Index
		'---------------------------------------------------------------------------------------------------------
			myfile = "https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv"
			covid_raw = read_csv(myfile)
			covid_raw$Date = ymd(covid_raw$Date)		
			covid = covid_raw %>% filter(CountryCode == "AUS") %>% rename(date=Date) %>% group_by(date=floor_date(date, "quarter")) %>% summarize(StringencyIndex=mean(StringencyIndex)) %>% mutate(date = zoo::as.yearqtr(date)) %>% select(date, StringencyIndex)
		
	xrun RSTAR_data = list(R_5206, R_6202, R_g1, R_g3, R_f1, covid) %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="date"), .)

	xoff
	pageselect RSTAR
	xget(type=series) RSTAR_data
	xclose

	'-----------------------------------------------------------------------------
	'Data Manipulations

	smpl @all

	'Rename Variables
		'Real GDP
		rename A2304402X RGDP
		'Employment
		rename a84423043c LE
		'Labor Force
		rename A84423047L LF
		'Trimmed Mean
		rename GCPIOCPMTMYP DLPTM
		'Bond Yield Price Expectations
		rename gbonyld pie_bond
		'Interbank Overnight Cash Rate
		rename FIRMMCRI CASH

		'Unemployment Rate
		series LUR = 100*(1-LE/LF)

		'Qtly Real GDP Log Level and Growth
		series lrgdp = log(rgdp)*100
		series dlrgdp = d(lrgdp)
		
		'Real Cash Rate (delfated with trimmed mean inflation)
		series rcash = cash - dlptm

		'Stringency Index of COVID
		rename StringencyIndex covid
		smpl @first 2019q4
		covid = 0
	
		smpl @all
		wfsave importeddata

		logmsg "Finished Data Import"

else

	wfopen importeddata
	
endif

smpl @all

'-----------------------------------------------------------------------------
'Creat initial series for initialzing state space

'Set Estimation Period
sample ssest 1986q3 2019q4

smpl ssest

lrgdp.hpf(lambda=1600) ypot_ini
lur.hpf(lambda=1600) nairu_ini
rcash.hpf(lambda=1600) nrate_ini
series ygap_ini = lrgdp - ypot_ini
series g_ini = d(ypot_ini)
series z_ini = nrate_ini - 4*g_ini

'-----------------------------------------------------------------------------
'Create Parameters/Starting Values for State Space Model

	' Setup coefficient vectors
	coef(10) delta
	coef(10) alpha
	coef(10) beta
	coef(10) theta
	coef(10) gamma
	coef(10) sigma
	
	'-----------------------------------------------------------------------------
	' Estimate initial coefficients
	
		'***********************************************	
		' Output Gap
		smpl ssest
		equation eq_ygap.ls ygap_ini = alpha(1)*ygap_ini(-1) + alpha(2)*ygap_ini(-2) + _
												alpha(3)/2*(rcash(-1) - nrate_ini(-1) + rcash(-2) - nrate_ini(-2))

		' Store estimates for later use
		!alpha1=alpha(1)
		!alpha2=alpha(2)
		!alpha3=alpha(3)
		!sigma1 = eq_ygap.@se 

		'***********************************************

		'***********************************************	
		' Okun's Law
		smpl ssest
		equation eq_okun.ls LUR = nairu_ini + _
											beta(1)*(0.4*ygap_ini + 0.3*ygap_ini(-1) + 0.2*ygap_ini(-2) + 0.1*ygap_ini(-3))

		' Store estimates for later use
		!beta1=beta(1)
		!sigma2 = eq_okun.@se 

		'***********************************************

		'***********************************************	
		' Phillips Curve
		smpl ssest
		equation eq_ptm.ls dlptm = (1-gamma(1))*pie_bond + _
												gamma(1)/3*(dlptm(-1) + dlptm(-2) + dlptm(-3)) + _
												gamma(2)*ygap_ini(-1)

		' Store estimates for later use
		!gamma1=gamma(1)
		!gamma2=gamma(2)
		!sigma3 = eq_ptm.@se 

		'***********************************************

		'***********************************************
		' GDP Trend
		smpl ssest
		equation eq_ypot.ls ypot_ini =delta(1)+ypot_ini(-1)
		
		' Store estimates for later use
		!delta1=delta(1)
		!sigma4=eq_ypot.@se

		'***********************************************
			
		'***********************************************
		' Trend Growth Rate
		smpl ssest
		equation eq_g.ls d(g_ini) = delta(2)
		!delta2=delta(2)
		
		' Store estimates for later use
		!sigma5=eq_g.@se

		'***********************************************

		'***********************************************
		' Unexplained Neutral Rate Movements
		smpl ssest
		equation eq_z.ls d(z_ini) = delta(3)
		!delta3=delta(3)
		
		' Store estimates for later use
		!sigma6=eq_z.@se

		'***********************************************
		
		'***********************************************
		' LUR Trend Rate
		smpl ssest
		equation eq_NRATE.ls d(nairu_ini) = c(1)
		!sigma7=eq_NRATE.@se 

	!theta1 = @mean(g_ini)

'-----------------------------------------------------------------------------------------------------------------------------------------------
'Set constraints on parameters
	'Upper bound on a_3 parameter (slope of the IS curve)
'		!a_r_constraint_stage1 = NA
'		!a_r_constraint_stage2 = -0.0025
'		!a_r_constraint_stage3 = -0.0025
		
	'Lower bound on b_2 parameter (slope of the Phillips curve)
'		!b_y_constraint_stage1 = NA
'		!b_y_constraint_stage2 = 0.025
'		!b_y_constraint_stage3 = 0.025
	
'-----------------------------------------------------------------------------------------------------------------------------------------------
'Stage 1 - No Real Rate term in IS Curve and Constant Trend GDP growth
'-----------------------------------------------------------------------------------------------------------------------------------------------
'State Space System of GDP/ UR Rate and Part Rate
	sspace ss_stage1
	ss_stage1.append @param alpha(1) !alpha1 alpha(2) !alpha2 alpha(3) !alpha3
	ss_stage1.append @param beta(1) !beta1
	ss_stage1.append @param gamma(1) !gamma1 gamma(2) !gamma2
	ss_stage1.append @param sigma(1) !sigma1 sigma(2) !sigma2 sigma(3) !sigma3 sigma(4) !sigma4 sigma(6) !sigma6 sigma(7) !sigma7 
	ss_stage1.append @param theta(1) !theta1
	
	ss_stage1.append @signal lrgdp = ygap + ypot
	
	ss_stage1.append @signal LUR = NAIRU + _
											beta(1)*(0.4*ygap + 0.3*ygapL1 + 0.2*ygapL2 + 0.1*ygapL3) + _
											[ename = e2, var = (sigma(2)^2)]

	ss_stage1.append @signal dlptm = (1-gamma(1))*pie_bond + _
													gamma(1)/3*(dlptm(-1) + dlptm(-2) + dlptm(-3)) + _
													gamma(2)*ygapL1 + _
													[ename = e3, var = (sigma(3)^2)]
	
	ss_stage1.append @state ygap = alpha(1)*ygap(-1) + alpha(2)*ygapL1(-1) + _
													[ename = e1, var = (sigma(4)^2)]
	ss_stage1.append @state ygapL1 = ygap(-1)
	ss_stage1.append @state ygapL2 = ygapL1(-1)
	ss_stage1.append @state ygapL3 = ygapL2(-1)

	ss_stage1.append @state NAIRU = NAIRU(-1) + [ename = e4, var = (sigma(4)^2)]
	ss_stage1.append @state NAIRUL1 = NAIRU(-1)

	ss_stage1.append @state ypot = ypot(-1) + theta(1) + [ename = e5, var = (sigma(5)^2)]

	vector(7) mprior = 0
	mprior(1) = ygap_ini(@ifirst(ygap_ini))
	mprior(2) = ygap_ini(@ifirst(ygap_ini))
	mprior(3) = ygap_ini(@ifirst(ygap_ini))
	mprior(4) = ygap_ini(@ifirst(ygap_ini))
	mprior(5) = 5
	mprior(6) = 5
	mprior(7) = ypot_ini(@ifirst(ypot_ini))
	
	'From RBA Paper Posterior Modes
	sym(7) vprior = 0
	vprior(1,1) = 0.38
	vprior(2,2) = 0.38
	vprior(3,3) = 0.38
	vprior(4,4) = 0.38
	vprior(5,5) = 0.15
	vprior(6,6) = 0.15
	vprior(7,7) = 0.54
	
	ss_stage1.append @mprior mprior
	ss_stage1.append @vprior vprior

	logmsg "Estimating Stage 1 State Space"
	
	smpl ssest
	ss_stage1.ml(optmethod=legacy)
	ss_stage1.ml(optmethod=opg)
	ss_stage1.makestates(t=smooth) *_stage1

'***********************************************************************************************************************************
'Median Unbiased Estimator of Lambda
'***********************************************************************************************************************************

logmsg "Estimating First Stage Median Unbiased Estimator"

call MUE_Stage1(YPOT_stage1)

'***********************************************************************************************************************************
'Stage 2 - Impose Lamdag and Add R to Output gap Equation - Constant z
'***********************************************************************************************************************************
	sspace ss_stage2
	ss_stage2.append @param alpha(1) !alpha1 alpha(2) !alpha2 alpha(3) !alpha3
	ss_stage2.append @param beta(1) !beta1
	ss_stage2.append @param gamma(1) !gamma1 gamma(2) !gamma2
	ss_stage2.append @param sigma(1) !sigma1 sigma(2) !sigma2 sigma(3) !sigma3 sigma(4) !sigma4 sigma(6) !sigma6 sigma(7) !sigma7 
	ss_stage2.append @param theta(1) !theta1
	
	ss_stage2.append @signal lrgdp = ygap + ypot
	
	ss_stage2.append @state ygap = alpha(1)*ygap(-1) + alpha(2)*ygapL1(-1) + _
													alpha(3)*((rcash(-1) + rcash(-2))/2 - 4*(g(-1)+gL1(-1))/2) - _
													alpha(3)*theta(2) + _
													[ename = e1, var = (sigma(1)^2)]
	ss_stage2.append @state ygapL1 = ygap(-1)
	ss_stage2.append @state ygapL2 = ygapL1(-1)
	ss_stage2.append @state ygapL3 = ygapL2(-1)

	ss_stage2.append @signal LUR = NAIRU + _
											beta(1)*(0.4*ygap + 0.3*ygapL1 + 0.2*ygapL2 + 0.1*ygapL3) + _
											[ename = e2, var = (sigma(2)^2)]

	ss_stage2.append @signal dlptm = (1-gamma(1))*pie_bond + _
													gamma(1)/3*(dlptm(-1) + dlptm(-2) + dlptm(-3)) + _
													gamma(2)*ygapL1 + _
													[ename = e3, var = (sigma(3)^2)]

	ss_stage2.append @state NAIRU = NAIRU(-1) + [ename = e4, var = (sigma(4)^2)]
	ss_stage2.append @state NAIRUL1 = NAIRU(-1)

	ss_stage2.append @state ypot = ypot(-1) + g(-1) + [ename = e5, var = (sigma(5)^2)]

	ss_stage2.append @state g = g(-1) + [ename = e6, var = (lambda_g*sigma(5))^2]
	ss_stage2.append @state gL1 = g(-1)

	vector(9) mprior = 0
	mprior(1) = ygap_ini(@ifirst(ygap_ini))
	mprior(2) = ygap_ini(@ifirst(ygap_ini))
	mprior(3) = ygap_ini(@ifirst(ygap_ini))
	mprior(4) = ygap_ini(@ifirst(ygap_ini))
	mprior(5) = 5
	mprior(6) = 5
	mprior(7) = ypot_ini(@ifirst(ypot_ini))
	mprior(8) = g_ini(@ifirst(g_ini))
	mprior(9) = g_ini(@ifirst(g_ini))
	
	'From RBA Paper Posterior Modes
	sym(9) vprior = 0
	vprior(1,1) = 0.38
	vprior(2,2) = 0.38
	vprior(3,3) = 0.38
	vprior(4,4) = 0.38
	vprior(5,5) = 0.15
	vprior(6,6) = 0.15
	vprior(7,7) = 0.54
	vprior(8,8) = 0.05
	vprior(9,9) = 0.05
	
	ss_stage2.append @mprior mprior
	ss_stage2.append @vprior vprior

	logmsg "Estimating Stage 2 State Space"
	
	smpl ssest
	ss_stage2.ml(optmethod=legacy)
	ss_stage2.ml(optmethod=opg)
	ss_stage2.makestates(t=smooth) *_stage2

	'Save parameter from IS Curve
	!ar = alpha(3)

  'Inputs for median unbiased estimator stage2
   group g_mue2 ygap_stage2 ygap_stage2(-1) ygap_stage2(-2) (rcash(-1)+rcash(-2))/2 4*(g_stage2(-1)+g_stage2(-2))/2 1
   stom(g_mue2,x_mue2)

'***********************************************************************************************************************************
'Median Unbiased Estimator of Lambda
'***********************************************************************************************************************************

'Estimating Median Unbiased Estimate of Lambda Z
call MUE_Stage2(x_mue2)

'***********************************************************************************************************************************
'Stage 3
'***********************************************************************************************************************************
	sspace ss_stage3
	ss_stage3.append @param alpha(1) !alpha1 alpha(2) !alpha2 alpha(3) !alpha3
	ss_stage3.append @param beta(1) !beta1
	ss_stage3.append @param gamma(1) !gamma1 gamma(2) !gamma2
	ss_stage3.append @param sigma(1) !sigma1 sigma(2) !sigma2 sigma(3) !sigma3 sigma(4) !sigma4 sigma(6) !sigma6 sigma(7) !sigma7 
	ss_stage3.append @param theta(1) !theta1
	
	ss_stage3.append @signal lrgdp = ygap + ypot
	
	ss_stage3.append @state ygap = alpha(1)*ygap(-1) + alpha(2)*ygapL1(-1) + _
													alpha(3)*((rcash(-1) + rcash(-2))/2 - 4*(g(-1)+gL1(-1))/2) - _
													alpha(3)*(z(-1) + zL1(-1))/2 + _
													[ename = e1, var = (sigma(1)^2)]
	ss_stage3.append @state ygapL1 = ygap(-1)
	ss_stage3.append @state ygapL2 = ygapL1(-1)
	ss_stage3.append @state ygapL3 = ygapL2(-1)

	ss_stage3.append @signal LUR = NAIRU + _
											beta(1)*(0.4*ygap + 0.3*ygapL1 + 0.2*ygapL2 + 0.1*ygapL3) + _
											[ename = e2, var = (sigma(2)^2)]

	ss_stage3.append @signal dlptm = (1-gamma(1))*pie_bond + _
													gamma(1)/3*(dlptm(-1) + dlptm(-2) + dlptm(-3)) + _
													gamma(2)*ygapL1 + _
													[ename = e3, var = (sigma(3)^2)]

	ss_stage3.append @state NAIRU = NAIRU(-1) + [ename = e4, var = (sigma(4)^2)]
	ss_stage3.append @state NAIRUL1 = NAIRU(-1)

	ss_stage3.append @state ypot = ypot(-1) + g(-1) + [ename = e5, var = (sigma(5)^2)]

	ss_stage3.append @state g = g(-1) + [ename = e6, var = (lambda_g*sigma(5))^2]
	ss_stage3.append @state gL1 = g(-1)

	ss_stage3.append @state z = z(-1) + [ename = e7, var = ((lambda_z*sigma(1)/!ar)^2)]
	ss_stage3.append @state zL1 = z(-1)

	ss_stage3.append @state nrate = 4*(g(-1) + e7) + (z(-1) + e6)

	vector(12) mprior = 0
	mprior(1) = ygap_ini(@ifirst(ygap_ini))
	mprior(2) = ygap_ini(@ifirst(ygap_ini))
	mprior(3) = ygap_ini(@ifirst(ygap_ini))
	mprior(4) = ygap_ini(@ifirst(ygap_ini))
	mprior(5) = 5
	mprior(6) = 5
	mprior(7) = ypot_ini(@ifirst(ypot_ini))
	mprior(8) = g_ini(@ifirst(g_ini))
	mprior(9) = g_ini(@ifirst(g_ini))
	mprior(10) = z_ini(@ifirst(z_ini))
	mprior(11) = z_ini(@ifirst(z_ini))
	mprior(12) = nrate_ini(@ifirst(nrate_ini))
	
	'From RBA Paper Posterior Modes
	sym(12) vprior = 0
	vprior(1,1) = 0.38
	vprior(2,2) = 0.38
	vprior(3,3) = 0.38
	vprior(4,4) = 0.38
	vprior(5,5) = 0.15
	vprior(6,6) = 0.15
	vprior(7,7) = 0.54
	vprior(8,8) = 0.05
	vprior(9,9) = 0.05
	vprior(10,10) = 0.22
	vprior(11,11) = 0.22
	
	ss_stage3.append @mprior mprior
	ss_stage3.append @vprior vprior

	logmsg "Estimating Stage 3 State Space"
	
	smpl ssest
	ss_stage3.ml(optmethod=legacy)
	ss_stage3.ml(optmethod=opg)
	ss_stage3.makestates(t=smooth) *_stage3
	ss_stage3.makestates(t=smoothse) *_stage3se

'-----------------------------------------------------------------------------
'Plotting Results
!bound1=@qnorm(0.85)
!bound2=@qnorm(0.95)

smpl ssest
series low70=(nrate_stage3-!bound1*nrate_stage3se)
series high70=(nrate_stage3+!bound1*nrate_stage3se)
series low90=(nrate_stage3-!bound2*nrate_stage3se)
series high90=(nrate_stage3+!bound2*nrate_stage3se)

smpl 1994q1 @last
group g_NRATE low90 high90 low70 high70 NRATE_stage3 RCASH
freeze(p_NRATE) g_NRATE.mixed band(1,2,3,4) line(5)
p_NRATE.setelem(1) fillcolor(@rgb(16, 189, 239))
p_NRATE.setelem(2) fillcolor(@rgb(14, 139, 241))
p_NRATE.setelem(1) lcolor(black)
p_NRATE.setelem(2) lcolor(red)
p_NRATE.name(1) 90 per cent confidence interval
p_NRATE.name(2) 
p_NRATE.name(3) 70 per cent confidence interval
p_NRATE.name(4) 
p_NRATE.name(5) Neutral Rate of Interest
p_NRATE.name(6) Real Cash Rate
p_NRATE.legend display position(0.5,2)
p_NRATE.options gridnone
show p_NRATE

p_NRATE.save(t=pdf) NeutralRate_Chart.pdf


