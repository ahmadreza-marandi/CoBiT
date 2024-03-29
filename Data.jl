# Julia 1.0.5
#This document contains data used in the paper
# A. Marandi, V Lurkin (2020), Static Pricing Problems under Mixed Multinomial Logit Demand
#written by Ahmadreza Marandi
# All rights are reserved.
#
#
using MAT
function random_instance_N_1(alter::Int64,R_AT::Int64)
	NUM_POINTS=alter+1; #alternatives

	a_ik=[zeros(1,R_AT);
	round.(rand(alter,R_AT)*5-ones(alter,R_AT),digits=5)];

	b_ik=[zeros(1,R_AT);
	round.(rand(alter,R_AT)*5+0.025*ones(alter,R_AT),digits=5)];
	N=1; #customers

	R=R_AT;
	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);

	price=100*ones(alter);
	UB_p=[0;price];
	LB_p=zeros(NUM_POINTS);
	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	q_parameter[:,1,:]=a_ik;
	
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R);
	Beta_parameter[:,1,:]=-b_ik;
	prob_v=rand(1,R_AT);
	w_k=prob_v/sum(prob_v);
 	# 	 file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/MixedLogitn10_r50.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
    return Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p,w_k;
end
function Intel_cooperation()
	NUM_POINTS=3+1; #alternatives
	a_ik=[4 0	0	0	0	0	0	0
		  1 -1.0334 3.2480 -0.9336 1.7094 0.4187 -0.8904 -0.9804
		  2 0.7840 4.7161 -0.3438 1.8777 2.1771 -0.4310 -0.4907
		  3 6.0054 3.8771 1.3506 2.3611 1.1723 0.8889 0.9163
				  ];

	b_ik=[4 0	0	0	0	0	0	0
		  1 0.00416 0.01840 0.00525 0.01165 0.01015 0.00325 0.00331
		  2 0.00312 0.01354 0.00394 0.00874 0.00639 0.00244 0.00248
		  3 0.00181 0.00744 0.00229 0.00508 0.00167 0.00142 0.00144];
	N=1; #customers

	R=7;
	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);

	price=[207; 299; 410];
	UB_p=[0;price*3];
	LB_p=[0
	0
	0
	0];
	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	q_parameter[:,1,:]=a_ik[:,2:end];
	
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R);
	Beta_parameter[:,1,:]=-b_ik[:,2:end];
	w_k=[0.0753 0.1126 0.1285 0.1180 0.0859 0.2842 0.1953];
 	# 	 file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/MixedLogitn10_r50.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
    return Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p,w_k;
end

function Logit_10()
	NUM_POINTS=2+1; #alternatives
	N=10; #customers

	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];

	#Parameters choice model
	ASC_PSP = 32;
	ASC_PUP = 34;
	Beta_TD = -0.612;
	Beta_Origin = -5.762;
	Beta_Age_Veh = 4.037;
	Beta_FEE_INC_PSP = -10.995;
	Beta_FEE_RES_PSP = -11.440;
	Beta_FEE_INC_PUP = -13.729;
	Beta_FEE_RES_PUP = -10.668;

	Beta_AT = -0.788;
	Beta_FEE = -32.328;
			 
	 #Variables choice model
	 AT_FSP = 10;
	 TD_FSP = 10;
	 AT_PSP = 10;
	 TD_PSP = 10;
	 AT_PUP = 5;
	 TD_PUP = 10;
	 Origin =[0
	 		1
	 		1
	 		0
	 		0
	 		1
	 		0
	 		0
	 		1
	 		0;
	];
	 Age_veh =[	0
	 		0
	 		0
	 		1
	 		0
	 		0
	 		1
	 		0
	 		0
	 		0];
	 Low_inc =[	1
	 		1
	 		1
	 		1
	 		0
	 		1
	 		1
	 		1
	 		1
	 		0];
	 Res =[ 1
	 		1
	 		1
	 		0
	 		1
	 		1
	 		0
	 		0
	 		1
	 		1];

	q_parameter=ones(NUM_POINTS,N);
	for n=1:N
		q_parameter[1,n]=		   Beta_AT * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
		q_parameter[2,n]=ASC_PSP + Beta_AT * AT_PSP +  Beta_TD * TD_PSP;
		q_parameter[3,n]=ASC_PUP + Beta_AT * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n]; 
	end

	Beta_parameter=ones(NUM_POINTS,N);

	for n=1:N
		Beta_parameter[1,n]=0;
		Beta_parameter[2,n]=Beta_FEE + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
		Beta_parameter[3,n]=Beta_FEE + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	end
	#				 
	 # file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/Logit_10.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
 	return Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p;
end


function Mixed_Logit_10(β)
	NUM_POINTS=2+1; #alternatives
	N=10; #customers

	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];

	#Parameters choice model
	ASC_PSP = 32;
	ASC_PUP = 34;
	Beta_TD = -0.612;
	Beta_Origin = -5.762;
	Beta_Age_Veh = 4.037;
	Beta_FEE_INC_PSP = -10.995;
	Beta_FEE_RES_PSP = -11.440;
	Beta_FEE_INC_PUP = -13.729;
	Beta_FEE_RES_PUP = -10.668;

	Beta_AT = -0.788;
	Beta_FEE = -32.328;
			 
	 #Variables choice model
	 AT_FSP = 10;
	 TD_FSP = 10;
	 AT_PSP = 10;
	 TD_PSP = 10;
	 AT_PUP = 5;
	 TD_PUP = 10;
	 Origin =[0
	 		1
	 		1
	 		0
	 		0
	 		1
	 		0
	 		0
	 		1
	 		0;
	];
	 Age_veh =[	0
	 		0
	 		0
	 		1
	 		0
	 		0
	 		1
	 		0
	 		0
	 		0];
	 Low_inc =[	1
	 		1
	 		1
	 		1
	 		0
	 		1
	 		1
	 		1
	 		1
	 		0];
	 Res =[ 1
	 		1
	 		1
	 		0
	 		1
	 		1
	 		0
	 		0
	 		1
	 		1];

	q_parameter=ones(NUM_POINTS,N);
	for n=1:N
		q_parameter[1,n]=		   β[1] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
		q_parameter[2,n]=ASC_PSP + β[1] * AT_PSP +  Beta_TD * TD_PSP;
		q_parameter[3,n]=ASC_PUP + β[1] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n]; 
	end

	Beta_parameter=ones(NUM_POINTS,N);

	for n=1:N
		Beta_parameter[1,n]=0;
		Beta_parameter[2,n]=β[2] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
		Beta_parameter[3,n]=β[2] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	end
	#				 
	 # file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/Logit_10.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
 	return Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p;
end
function Mixed_Logit_50(β)
	NUM_POINTS=2+1; #alternatives
	N=50; #customers

	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];

	#Parameters choice model
	ASC_PSP = 32;
	ASC_PUP = 34;
	Beta_TD = -0.612;
	Beta_Origin = -5.762;
	Beta_Age_Veh = 4.037;
	Beta_FEE_INC_PSP = -10.995;
	Beta_FEE_RES_PSP = -11.440;
	Beta_FEE_INC_PUP = -13.729;
	Beta_FEE_RES_PUP = -10.668;

	Beta_AT = -0.788;
	Beta_FEE = -32.328;
			 
	 #Variables choice model
	 AT_FSP = 10;
	 TD_FSP = 10;
	 AT_PSP = 10;
	 TD_PSP = 10;
	 AT_PUP = 5;
	 TD_PUP = 10;
	 Origin =[1	0
	2	1
	3	1
	4	0
	5	0
	6	1
	7	0
	8	0
	9	1
	10	0
	11	1
	12	1
	13	0
	14	0
	15	0
	16	1
	17	0
	18	0
	19	0
	20	1
	21	0
	22	0
	23	0
	24	0
	25	0
	26	0
	27	0
	28	0
	29	1
	30	1
	31	0
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	1
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Origin=Origin[:,2];
	 Age_veh =[1	0
	2	0
	3	0
	4	1
	5	0
	6	0
	7	1
	8	0
	9	0
	10	0
	11	0
	12	0
	13	0
	14	1
	15	0
	16	1
	17	1
	18	0
	19	1
	20	1
	21	1
	22	1
	23	0
	24	0
	25	1
	26	0
	27	0
	28	0
	29	0
	30	0
	31	1
	32	0
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	1
	40	0
	41	1
	42	0
	43	1
	44	1
	45	0
	46	1
	47	0
	48	0
	49	0
	50	0];
	Age_veh=Age_veh[:,2];
	 Low_inc =[	1	1
	2	1
	3	1
	4	1
	5	1
	6	0
	7	1
	8	1
	9	1
	10	1
	11	1
	12	0
	13	1
	14	0
	15	1
	16	1
	17	1
	18	1
	19	1
	20	0
	21	0
	22	1
	23	1
	24	1
	25	0
	26	1
	27	1
	28	1
	29	1
	30	1
	31	1
	32	1
	33	0
	34	1
	35	0
	36	1
	37	0
	38	1
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	0
	49	0
	50	0];
	Low_inc=Low_inc[:,2];
	 Res =[ 1	1
	2	1
	3	1
	4	0
	5	1
	6	1
	7	0
	8	0
	9	1
	10	1
	11	1
	12	1
	13	0
	14	0
	15	1
	16	1
	17	0
	18	0
	19	0
	20	1
	21	1
	22	0
	23	0
	24	0
	25	0
	26	1
	27	0
	28	0
	29	1
	30	1
	31	1
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Res=Res[:,2];

	q_parameter=ones(NUM_POINTS,N);
	for n=1:N
		q_parameter[1,n]=		   β[1] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
		q_parameter[2,n]=ASC_PSP + β[1] * AT_PSP +  Beta_TD * TD_PSP;
		q_parameter[3,n]=ASC_PUP + β[1] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n]; 
	end

	Beta_parameter=ones(NUM_POINTS,N);

	for n=1:N
		Beta_parameter[1,n]=0;
		Beta_parameter[2,n]=β[2] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
		Beta_parameter[3,n]=β[2] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	end
	#				 
	 # file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/Logit_10.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
 	return Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p;
end
function Logit_50()
	NUM_POINTS=2+1; #alternatives
	N=50; #customers

	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];

	#Parameters choice model
	ASC_PSP = 32;
	ASC_PUP = 34;
	Beta_TD = -0.612;
	Beta_Origin = -5.762;
	Beta_Age_Veh = 4.037;
	Beta_FEE_INC_PSP = -10.995;
	Beta_FEE_RES_PSP = -11.440;
	Beta_FEE_INC_PUP = -13.729;
	Beta_FEE_RES_PUP = -10.668;

	Beta_AT = -0.788;
	Beta_FEE = -32.328;
			 
	 #Variables choice model
	 AT_FSP = 10;
	 TD_FSP = 10;
	 AT_PSP = 10;
	 TD_PSP = 10;
	 AT_PUP = 5;
	 TD_PUP = 10;
	 Origin =[1	0
	2	1
	3	1
	4	0
	5	0
	6	1
	7	0
	8	0
	9	1
	10	0
	11	1
	12	1
	13	0
	14	0
	15	0
	16	1
	17	0
	18	0
	19	0
	20	1
	21	0
	22	0
	23	0
	24	0
	25	0
	26	0
	27	0
	28	0
	29	1
	30	1
	31	0
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	1
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Origin=Origin[:,2];
	 Age_veh =[1	0
	2	0
	3	0
	4	1
	5	0
	6	0
	7	1
	8	0
	9	0
	10	0
	11	0
	12	0
	13	0
	14	1
	15	0
	16	1
	17	1
	18	0
	19	1
	20	1
	21	1
	22	1
	23	0
	24	0
	25	1
	26	0
	27	0
	28	0
	29	0
	30	0
	31	1
	32	0
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	1
	40	0
	41	1
	42	0
	43	1
	44	1
	45	0
	46	1
	47	0
	48	0
	49	0
	50	0];
	Age_veh=Age_veh[:,2];
	 Low_inc =[	1	1
	2	1
	3	1
	4	1
	5	1
	6	0
	7	1
	8	1
	9	1
	10	1
	11	1
	12	0
	13	1
	14	0
	15	1
	16	1
	17	1
	18	1
	19	1
	20	0
	21	0
	22	1
	23	1
	24	1
	25	0
	26	1
	27	1
	28	1
	29	1
	30	1
	31	1
	32	1
	33	0
	34	1
	35	0
	36	1
	37	0
	38	1
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	0
	49	0
	50	0];
	Low_inc=Low_inc[:,2];
	 Res =[ 1	1
	2	1
	3	1
	4	0
	5	1
	6	1
	7	0
	8	0
	9	1
	10	1
	11	1
	12	1
	13	0
	14	0
	15	1
	16	1
	17	0
	18	0
	19	0
	20	1
	21	1
	22	0
	23	0
	24	0
	25	0
	26	1
	27	0
	28	0
	29	1
	30	1
	31	1
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Res=Res[:,2];

	q_parameter=ones(NUM_POINTS,N);
	for n=1:N
		q_parameter[1,n]=Beta_AT * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
		q_parameter[2,n]=ASC_PSP + Beta_AT * AT_PSP +  Beta_TD * TD_PSP;
		q_parameter[3,n]=ASC_PUP + Beta_AT * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n]; 
	end

	Beta_parameter=ones(NUM_POINTS,N);

	for n=1:N
		Beta_parameter[1,n]=0;
		Beta_parameter[2,n]=Beta_FEE + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
		Beta_parameter[3,n]=Beta_FEE + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	end
	 # file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/Logit50.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
	return Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p;
end

function Mixed_Logit_n10_random(R_AT)
	NUM_POINTS=2+1; #alternatives
	N=10; #customers
	R=R_AT;
	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];
	#Parameters choice model
	 ASC_PSP = 32;
	 ASC_PUP = 34;
	 Beta_TD = -0.612;
	 Beta_Origin = -5.762;
	 Beta_Age_Veh = 4.037;
	 Beta_FEE_INC_PSP = -10.995;
	 Beta_FEE_RES_PSP = -11.440;
	 Beta_FEE_INC_PUP = -13.729;
	 Beta_FEE_RES_PUP = -10.668;

	
	address="C:\\Users\\20176914\\surfdrive\\scientific\\University\\Codes\\Julia\\NonlinearProgramming_Virginie\\random beta\\Beta_N"*string(N)* "R_"*string(R_AT)*".mat";
	input=matread(address)
	Beta_AT=input["Beta_AT"];#AT coefficient (10 customers x R_AT draws)
	 Beta_FEE=input["Beta_FEE"];#FEE coefficient (10 customers x R_FEE draws)
	
	#Variables choice model
	  AT_FSP = 10;
	  TD_FSP = 10;
	  AT_PSP = 10;
	  TD_PSP = 10;
	  AT_PUP = 5;
	  TD_PUP = 10;
	  Origin =[
	 	1	0
	 	2	1
	 	3	1
	 	4	0
	 	5	0
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	0];
	  Origin= Origin[:,2:end]	;
	  Age_veh =[
	  	1	0
	 	2	0
	 	3	0
	 	4	1
	 	5	0
	 	6	0
	 	7	1
	 	8	0
	 	9	0
	 	10	0];
	  Age_veh= Age_veh[:,2:end];
	  Low_inc =[
	  	1	1
	 	2	1
	 	3	1
	 	4	1
	 	5	0
	 	6	1
	 	7	1
	 	8	1
	 	9	1
	 	10	0];
	  Low_inc= Low_inc[:,2:end];
	  Res =[
	 	1	1
	 	2	1
	 	3	1
	 	4	0
	 	5	1
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	1];
	   Res= Res[:,2:end];

	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	for n=1:N
	  	for r=1:R
	  		q_parameter[1,n,r] = Beta_AT[n,r] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
	  		q_parameter[2,n,r] = ASC_PSP + Beta_AT[n,r] * AT_PSP +  Beta_TD * TD_PSP;
			q_parameter[3,n,r] = ASC_PUP + Beta_AT[n,r] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n];	
	  	end 
	end
	  
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R);
	for n=1:N
	  	for r=1:R
	  		Beta_parameter[1,n,r]=0;
	  		Beta_parameter[2,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
	  		Beta_parameter[3,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	 	end 
	end
 	# 	 file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/MixedLogitn10_r50.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
    return Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p;

end
function Mixed_logit_n50_equi_distance( R )
	a=[-3.6;-68.52];# lowerbounds on the integral
	b=[1.94;3.92];#upperbounds on the integral
	NUM_POINTS=2+1; #alternatives
	N=50; #customers

	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];

	#Parameters choice model
	ASC_PSP = 32;
	ASC_PUP = 34;
	Beta_TD = -0.612;
	Beta_Origin = -5.762;
	Beta_Age_Veh = 4.037;
	Beta_FEE_INC_PSP = -10.995;
	Beta_FEE_RES_PSP = -11.440;
	Beta_FEE_INC_PUP = -13.729;
	Beta_FEE_RES_PUP = -10.668;


	#AT coefficient (10 customers x R_AT draws)

	Beta_AT=zeros(R);
	 
	#FEE coefficient (10 customers x R_FEE draws)

	μ=[-0.788;-32.3];# mean of the random varialbe
	Σ=[(1.06)^2 -12.8;-12.8 (14.2)^2];# covariance matrix of the random varialbe
	d=MvNormal( μ , Σ );
	Beta_FEE=zeros(R);
	w_k=zeros(R^2);
	eps_AT=0;
	eps_FEE=0;
	for r=1:R
		Beta_AT[r]=a[1]+(r-1)*(b[1]-a[1])/R;
		Beta_FEE[r]=a[2]+(r-1)*(b[2]-a[2])/R;
	end
	# for n=1:N
		# eps_AT=rand(1) *b[1]*0.01;
		# eps_FEE=rand(1) *b[1]*0.01;
		for r_1=1:R
			for r_2=1:R
				# for i=1:NUM_POINTS
					w_k[R*(r_1-1)+r_2]=pdf(d,[Beta_AT[r_1],Beta_FEE[r_2]])*prod(b-a)/R^2;
				# end
			end
		end
	# end
	 
	 #Variables choice model
	 AT_FSP = 10;
	 TD_FSP = 10;
	 AT_PSP = 10;
	 TD_PSP = 10;
	 AT_PUP = 5;
	 TD_PUP = 10;
	 Origin =[1	0
	2	1
	3	1
	4	0
	5	0
	6	1
	7	0
	8	0
	9	1
	10	0
	11	1
	12	1
	13	0
	14	0
	15	0
	16	1
	17	0
	18	0
	19	0
	20	1
	21	0
	22	0
	23	0
	24	0
	25	0
	26	0
	27	0
	28	0
	29	1
	30	1
	31	0
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	1
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Origin=Origin[:,2];
	 Age_veh =[1	0
	2	0
	3	0
	4	1
	5	0
	6	0
	7	1
	8	0
	9	0
	10	0
	11	0
	12	0
	13	0
	14	1
	15	0
	16	1
	17	1
	18	0
	19	1
	20	1
	21	1
	22	1
	23	0
	24	0
	25	1
	26	0
	27	0
	28	0
	29	0
	30	0
	31	1
	32	0
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	1
	40	0
	41	1
	42	0
	43	1
	44	1
	45	0
	46	1
	47	0
	48	0
	49	0
	50	0];
	Age_veh=Age_veh[:,2];
	 Low_inc =[	1	1
	2	1
	3	1
	4	1
	5	1
	6	0
	7	1
	8	1
	9	1
	10	1
	11	1
	12	0
	13	1
	14	0
	15	1
	16	1
	17	1
	18	1
	19	1
	20	0
	21	0
	22	1
	23	1
	24	1
	25	0
	26	1
	27	1
	28	1
	29	1
	30	1
	31	1
	32	1
	33	0
	34	1
	35	0
	36	1
	37	0
	38	1
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	0
	49	0
	50	0];
	Low_inc=Low_inc[:,2];
	 Res =[ 1	1
	2	1
	3	1
	4	0
	5	1
	6	1
	7	0
	8	0
	9	1
	10	1
	11	1
	12	1
	13	0
	14	0
	15	1
	16	1
	17	0
	18	0
	19	0
	20	1
	21	1
	22	0
	23	0
	24	0
	25	0
	26	1
	27	0
	28	0
	29	1
	30	1
	31	1
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Res=Res[:,2];

	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R^2);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	for n=1:N
	  	for r_1=1:R
	  		for r_2=1:R
		  		q_parameter[1,n,R*(r_1-1)+r_2] = Beta_AT[r_1] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
		  		q_parameter[2,n,R*(r_1-1)+r_2] = ASC_PSP + Beta_AT[r_1] * AT_PSP +  Beta_TD * TD_PSP;
				q_parameter[3,n,R*(r_1-1)+r_2] = ASC_PUP + Beta_AT[r_1] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n];	
			end
	  	end 
	end
	  
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R^2);
	for n=1:N
	  	for r_1=1:R
	  		for r_2=1:R
		  		Beta_parameter[1,n,R*(r_1-1)+r_2]=0;
		  		Beta_parameter[2,n,R*(r_1-1)+r_2] = Beta_FEE[r_2] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
		  		Beta_parameter[3,n,R*(r_1-1)+r_2] = Beta_FEE[r_2] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
		  	end
	 	end 
	end
 	# 	 file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/MixedLogitn10_r50.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
    return Beta_parameter,q_parameter,NUM_POINTS,N,R^2,UB_p,LB_p,w_k;
end
function Mixed_logit_n10_equi_distance( R )
	a=[-3.6;-68.52];# lowerbounds on the integral
	b=[1.94;3.92];#upperbounds on the integral
	NUM_POINTS=2+1; #alternatives
	N=10; #customers
	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];
	#Parameters choice model
	 ASC_PSP = 32;
	 ASC_PUP = 34;
	 Beta_TD = -0.612;
	 Beta_Origin = -5.762;
	 Beta_Age_Veh = 4.037;
	 Beta_FEE_INC_PSP = -10.995;
	 Beta_FEE_RES_PSP = -11.440;
	 Beta_FEE_INC_PUP = -13.729;
	 Beta_FEE_RES_PUP = -10.668;

	#AT coefficient (10 customers x R_AT draws)

	Beta_AT=zeros(R);
	 
	#FEE coefficient (10 customers x R_FEE draws)

	μ=[-0.788;-32.3];# mean of the random varialbe
	Σ=[(1.06)^2 -12.8;-12.8 (14.2)^2];# covariance matrix of the random varialbe
	d=MvNormal( μ , Σ );
	Beta_FEE=zeros(R);
	w_k=zeros(R^2);
	eps_AT=0;
	eps_FEE=0;
	for r=1:R
		Beta_AT[r]=a[1]+(r-1)*(b[1]-a[1])/R;
		Beta_FEE[r]=a[2]+(r-1)*(b[2]-a[2])/R;
	end
	# for n=1:N
		# eps_AT=rand(1) *b[1]*0.01;
		# eps_FEE=rand(1) *b[1]*0.01;
		for r_1=1:R
			for r_2=1:R
				# for i=1:NUM_POINTS
					w_k[R*(r_1-1)+r_2]=pdf(d,[Beta_AT[r_1],Beta_FEE[r_2]])*prod(b-a)/R^2;
				# end
			end
		end
	# end
	 
	#Variables choice model
	  AT_FSP = 10;
	  TD_FSP = 10;
	  AT_PSP = 10;
	  TD_PSP = 10;
	  AT_PUP = 5;
	  TD_PUP = 10;
	  Origin =[
	 	1	0
	 	2	1
	 	3	1
	 	4	0
	 	5	0
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	0];
	  Origin= Origin[:,2:end]	;
	  Age_veh =[
	  	1	0
	 	2	0
	 	3	0
	 	4	1
	 	5	0
	 	6	0
	 	7	1
	 	8	0
	 	9	0
	 	10	0];
	  Age_veh= Age_veh[:,2:end];
	  Low_inc =[
	  	1	1
	 	2	1
	 	3	1
	 	4	1
	 	5	0
	 	6	1
	 	7	1
	 	8	1
	 	9	1
	 	10	0];
	  Low_inc= Low_inc[:,2:end];
	  Res =[
	 	1	1
	 	2	1
	 	3	1
	 	4	0
	 	5	1
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	1];
	   Res= Res[:,2:end];

	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R^2);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	for n=1:N
	  	for r_1=1:R
	  		for r_2=1:R
		  		q_parameter[1,n,R*(r_1-1)+r_2] = Beta_AT[r_1] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
		  		q_parameter[2,n,R*(r_1-1)+r_2] = ASC_PSP + Beta_AT[r_1] * AT_PSP +  Beta_TD * TD_PSP;
				q_parameter[3,n,R*(r_1-1)+r_2] = ASC_PUP + Beta_AT[r_1] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n];	
			end
	  	end 
	end
	  
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R^2);
	for n=1:N
	  	for r_1=1:R
	  		for r_2=1:R
		  		Beta_parameter[1,n,R*(r_1-1)+r_2]=0;
		  		Beta_parameter[2,n,R*(r_1-1)+r_2] = Beta_FEE[r_2] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
		  		Beta_parameter[3,n,R*(r_1-1)+r_2] = Beta_FEE[r_2] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
		  	end
	 	end 
	end
 	# 	 file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/MixedLogitn10_r50.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
    return Beta_parameter,q_parameter,NUM_POINTS,N,R^2,UB_p,LB_p,w_k;
end
function Mixed_Logit_n10_r50()
	NUM_POINTS=2+1; #alternatives
	N=10; #customers
	R=50;
	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];
	#Parameters choice model
	 ASC_PSP = 32;
	 ASC_PUP = 34;
	 Beta_TD = -0.612;
	 Beta_Origin = -5.762;
	 Beta_Age_Veh = 4.037;
	 Beta_FEE_INC_PSP = -10.995;
	 Beta_FEE_RES_PSP = -11.440;
	 Beta_FEE_INC_PUP = -13.729;
	 Beta_FEE_RES_PUP = -10.668;

	#AT coefficient (10 customers x 50 draws)
	 Beta_AT=[
			 1 	 -1.2126571130  -0.8841591702   0.6042890019  -3.7638897994  -1.1085836586   0.8441170849  -1.7691852030  -0.8929604654   0.0844830025  -0.7357084861  -1.4237054695  -0.6821856429  -0.2116798356  -0.8128024620  -2.1742897620  -2.4173483150  -0.9117324286  -0.6082732346  -1.7663928392   0.2309322424   0.0261123788   1.3749771682  -1.2672851860  -0.2609197564   1.1842437850   0.2315795724   0.8503641618  -0.4169073839  -1.1187808931  -0.6225571217  -0.4421313186  -1.7181571770  -0.9869521029  -0.3283009990  -1.8096812271   1.5704833172  -2.2684050845  -0.8048123850   0.3730016949   1.8304809839   0.7934353622   0.6487468283  -2.0276615604  -1.6836837513  -1.6027640875   0.4865778123  -0.1278757411   0.0217923832  -0.4818260683  -0.1808148752 
			 2 	 -1.4577389546  -2.4718974914  -1.2370395354  -1.8665105100  -0.3763215297   0.3189424883  -0.7351874532  -0.4642126598  -2.4128207015   0.5430678057  -0.4479676279  -0.3905352154  -0.0098322953  -0.8765042474  -1.4662179397  -0.7907044320   0.4426866475  -0.8921245633  -2.1133220169  -2.5609162487  -2.5477590580  -0.8130378527   0.2044585037  -1.1006870784  -0.8384144780   0.7214443316  -0.5343669813  -0.8185175493  -1.4839967433  -0.6886220310  -2.1662545440  -0.0870683935   0.0573018778  -0.6455135089  -1.5221428218   0.2339030652  -0.9806943139   1.6786897304  -1.8828890341  -1.6995949777  -1.1717529402  -1.5593735139  -0.5995371827  -0.0010689381  -1.9371884694  -0.8331061218   1.0221748280  -1.1541507036  -1.3674343716   0.6720465198 
			 3 	 -2.3463923750  -2.6706117778  -0.4812893806  -0.1091993068   0.4357542122  -0.6518895503  -0.7403285526   1.9391803440  -1.4566853674  -1.0449337754  -1.0970918657  -0.6963978648   1.4210538565  -0.8794586704  -1.0778194169   0.7871957437  -0.6488512625  -0.7967783717   2.2212196635  -0.8549702796  -1.8137786923   1.0442253319  -1.5296312999  -0.8380453620   1.0459229055   0.6032173332  -0.9639579138  -0.6550194436  -1.2496510812  -1.7430555933   0.1972353414   0.4362082827   0.5834559398  -1.7985431121  -0.2831464752  -1.0762358961  -0.5550074734  -0.3172855708  -0.4049150726   1.4279360594  -1.5631150462   0.6268687018   0.3900374929  -2.3364308586   0.6055633314  -0.4349898056   2.8181715555  -0.9369789930  -2.3062940332   1.9381197847 
			 4 	 -3.1803579316  -0.8640949287  -0.4549893214  -1.3594283347  -0.3774009923  -2.3724145259  -1.6196189221  -2.0392518276  -1.9823374854   0.1610464251  -1.4730761161  -1.0919535723  -0.1480289662  -2.3249364495   0.5836797644  -0.0052504277  -0.0066947309  -0.9520335946  -1.5960203311  -1.0495924530  -0.2691522544  -0.4396997258   0.8511143633  -0.6780829140  -2.9346436354  -0.9408584745  -0.9135186094  -0.0646495137   2.5701116467  -0.2609124995  -1.7322374306   0.8768350395   0.3300062981  -2.1156182140  -0.9678780789  -1.9801035291  -1.4978790228  -1.3969660762   0.1080367589  -1.5322852734   0.5715467161  -0.4073468106  -0.7157214927  -0.3872886604  -0.4313422152   0.9100259210   0.3253602757  -0.5256078490   0.0843072077  -0.9737859865 
			 5 	 -0.7483639454  -2.5019833001  -1.5281193739  -2.0295497815   0.1446445267  -2.4279790994  -1.3400313218   0.2639495294  -3.3339657843  -1.3405780074  -0.7222189117  -2.1017624906  -2.2428095819  -1.1167164642   0.6477383076  -1.9064738907   0.2479395396  -2.8793573745  -0.6583004426  -0.7762274072  -1.0867179702  -1.1003416474  -0.8214084665  -1.3531758996   0.4088115674  -1.2355548950  -0.7635857319  -0.3051112070  -1.1245757139  -1.0332181895   0.1858045896   0.8323285071   0.6346296993   0.0280534617  -1.2345993515  -0.5677202345  -0.3436223418  -2.1783156532  -0.7303045108  -1.0883842723   0.0513004853   0.2073123433  -1.3394706923  -0.1403307770  -2.5338384788  -0.6169750600   0.2669972682  -0.5583763161  -1.1731847411   0.0088264887 
			 6 	  1.0218377716  -0.7222146871   0.6498573241  -0.2263375424  -1.6186191230  -0.4070714848  -1.6358508475   0.4302110939   0.0980674752  -2.3382033330  -1.9308167417  -0.6393466759   0.3215053574  -1.3125375392  -0.4676117175  -1.0647289139  -2.5264437810  -2.5096292071   0.1001757687  -1.3101873035  -1.4658382354  -0.1966366038   0.4925547237  -0.8955363491  -0.2372179055  -1.6796307976   0.4397709867  -1.1777081212  -2.0947660472  -1.5735918982  -0.3717423073  -1.0903227362   0.1339323476  -0.1203781378  -2.3088223521   0.8762729383  -1.8552363671  -1.6067035539  -0.9123927188  -2.3785068937  -1.2567976120  -0.4489950769  -0.5605472207  -0.6629087610  -1.4201301696   0.3490493325  -0.4951849690  -0.1334786888  -1.8477679967  -0.9142322539 
			 7 	 -0.4225432544  -2.2424129808  -1.0586531914  -0.6825084671  -0.6750739198   0.0998712319  -1.4153379049   0.1507653055   0.1321299319  -0.4241683802  -1.6198856622  -0.0749339237  -1.8545832929  -0.8472248351  -1.4860823605  -1.2977606893  -2.3807661280  -2.0049418946   1.7004798991  -0.3208809549  -0.6482128835  -0.2729273847  -2.5242983891  -0.8260278466  -1.1905762286  -1.5646521754  -1.6086837801   0.0303197385  -1.2650806977  -1.1766050057  -0.6562623838  -1.3696037008  -1.6628232599  -1.8084632160  -0.9035540410  -0.3378253609  -1.4785124473   0.0437476962  -0.5183867557   0.7355037021  -0.6202806803  -0.3052441233  -2.5445514433   1.0179764346   0.5505294405  -0.6334681665  -0.5505328019  -1.8583701576   0.1363408332  -0.0902977356 
			 8 	 -1.1569867128  -0.6675343469  -2.8674853293  -0.6486729580  -0.7931033246  -1.3657033775  -0.8855316525  -0.2209299540   0.7064659417  -2.0801998906   0.8582046185   0.6206078584   0.2538604577  -0.3554363075  -0.0075752063  -0.6758406490  -0.1039973248  -0.3698945684   1.7400871805  -0.5659400235  -0.1972350357   0.7967810758   0.4622365665  -0.3346160617  -3.4533452489  -1.2954679677  -1.7739340736   0.9236973633  -0.6306058700  -2.7901734944  -1.4046679621  -0.9441406426   0.8284576440  -1.5593221019   0.6121510011  -0.3992449128  -2.0870217017  -1.4345303177  -0.7799014232  -0.1602322092  -3.3173075531  -0.4439234281  -1.6898151785   0.9314972778  -1.5713351192  -2.3109355093  -0.1444224142  -0.4276553123  -0.6474598076  -0.1379332444 
			 9 	 -1.6518559996  -3.1292551115  -1.3461301745  -0.8324745405  -3.9668429279  -0.5865167743  -1.4172844288   0.4348788489  -0.7368669823  -0.9665037356  -1.0578917931  -0.4122226400  -2.0812533848  -1.6648019213   0.7932650633  -1.1020797735  -2.0558981455  -1.1602789929  -1.6120880334  -1.0980126879  -0.2302223033  -0.5113692862   0.6511768093  -0.6654856598  -2.0071251393  -0.6558445406  -0.1018655929  -1.5199985046   0.8099104034   2.5042507112  -2.6296450540   0.8045582007  -2.3774743381   0.8613302540  -0.1353087052  -0.7492882879  -0.1627499772   1.3513551684  -0.1238469352   1.2413403149  -3.5943892633  -1.1837835063  -1.1716911729  -2.5505668937  -0.4301811650  -0.4536736515  -0.4157332231  -2.3807661280   1.0952694326   0.2092207016 
			 10 	 -1.9965814456  -1.7741252395  -0.4023497705  -1.3532624216   0.7303389760  -1.0668013191  -0.5017551642   0.4204171990  -0.9848710271  -0.5820792025  -0.7483229894  -2.1606500250  -0.4173949369  -0.9594497816  -1.3028033280  -1.3182891814  -1.4084897596   0.6190935364   0.2255118328   0.1010232254  -1.1552049201  -0.0539491044  -0.5568415627  -1.4914223675  -0.4204002929  -0.8157467977  -2.1123186208   0.2353255247  -2.0602871328  -3.0915109074  -1.3715865167  -0.7951811839   1.0094468865  -0.7927113114   0.2876115448  -0.5674823886  -0.7708058157  -1.3343510843   0.1372451353   0.4354598191  -2.7909374105   1.3666932355  -1.9018599629  -0.7229759898   1.1062875299  -0.1955869515  -0.6975620924  -0.7487099379  -2.6371369392   0.4166111322] ;
	 Beta_AT= Beta_AT[:,2:end];
	#FEE coefficient (10 customers x 50 draws)
	 Beta_FEE=[ 
			 1 	-39.8356007938 -16.9411785516 -23.7614515388 -24.3408358598 -41.7506407898  -9.4776606004 -31.6035818593 -38.2003444536 -21.5839966334 -30.5944229047 -23.6417660370 -55.8958590244 -35.8630764258 -31.4499436322 -56.1508161671 -25.0709639293 -14.0214836121 -31.4519257440 -30.2687172769 -64.6611470775   1.3668558986 -56.4529994108 -22.6517361155 -40.7596962103 -10.1257334779 -16.5801833804  -0.5061922425 -65.1179219797 -23.1825251409 -30.0288305536 -40.1736508254 -25.2585359441 -25.9595985196 -54.0800075282 -12.5188756982 -23.8161532323 -11.4485462051 -44.1466767558 -29.9729533162 -16.9103805000 -39.1207528572 -11.6983066309 -17.2112289265 -40.1213111434 -36.1251483039 -34.3698935382 -32.2921795105 -45.3907834847 -16.8195475860 -17.8982243687 
			 2 	-51.3366844422 -28.9341744797 -48.5243590340 -35.4999711446 -18.7949000273 -41.4109121062 -30.6882344046 -41.0583376351 -33.2360944987 -53.8444549709 -42.9979521876 -29.9838275655 -33.0306462990 -36.9631766753 -21.1019228224 -31.2012178760 -33.3381172097 -31.1275882018 -48.6072511625 -40.0745736930 -24.1002819916 -21.1576931332 -23.4251478745 -56.0844573132 -62.0337772895 -23.6207622231 -61.0764433066 -29.7909731708 -36.6747874372 -18.9830155390 -45.0431972239 -42.7917866458 -28.2979835766 -32.9669556386 -30.9928259828 -24.3691633809 -43.3188219936 -32.6710571467 -48.2175907622 -28.5649525911 -11.1411234500 -27.1248100490 -26.5252245304 -31.0413875772 -50.1563549732 -22.9907164725 -14.4001903562 -40.6931389898 -36.5836029832 -18.2014444433 
			 3 	-29.2142931137 -33.0898800544 -35.4705886626 -23.8997123234 -24.0891135362  -9.9547211888 -51.8873514440 -30.4648627148 -41.0640068881 -37.2169937439 -21.7520578115 -46.2635399895 -51.8724392064  -7.8321969255 -29.4372770236 -55.1928844134 -11.8916377365 -50.4827243407 -42.8856107640 -29.6427973159 -35.9845231788 -31.1274490676  -7.0587724416 -15.8559434744 -28.6950129166 -49.5024436838  -6.5471726285 -28.4734921802 -43.2159405632 -21.3120079772 -23.9462570802 -33.1013675650 -48.6670962299 -27.5220307237 -29.7249301750 -23.9672102910 -25.4600312078 -46.8994687292 -25.5364497891 -42.1487255418 -25.6868742871 -32.8672775859 -53.4339028172 -36.6871880843 -45.6734549114 -24.3941212349 -22.6518182623 -20.9791627507 -25.4101569822 -36.0436918597 
			 4 	-40.0366050838 -52.2876480011 -23.1921886959 -35.0422760710 -71.5904243934 -26.3740796764 -39.3526054855  -5.1236913507 -32.1350413255 -37.2136563486 -24.8168721135 -62.5966588965  -8.2405237257 -17.6457581778 -48.3022408205  -2.5631645625 -20.7945124994 -39.5135888125 -31.0844463195 -28.7155194819 -40.8850920481 -51.4254450855 -68.4594138817 -25.1656183656 -19.8591719441 -22.7366301120 -42.4780191291 -43.1512083007 -38.0090431717 -31.6030930434 -43.4351965837 -20.6340122917 -26.0413134213 -14.0721523233 -34.9218648949 -22.8532021305 -36.0741473029 -26.9209973335 -40.7734398201 -43.6319760610 -14.4719249825 -39.1053252582 -39.7211352417 -32.0610614927 -40.1513437686 -35.6369528824 -39.6370163590 -21.8482021519 -24.6949567040 -54.4036482736 
			 5 	-37.1839377142 -49.6057518501 -32.0806405908 -29.7677911092 -47.6825967787 -15.9997113778 -51.5021064221 -45.7483992731 -40.3010451459 -31.7328585702 -24.2782716452 -37.3541377459 -35.2120328272 -45.3709909182 -24.3167130893 -42.9160366595 -48.0457102506 -17.5095414193 -35.1446841107 -13.2940766572 -43.5868082916 -13.3363309325 -41.2446356950 -10.0127649458 -42.7425703862 -10.3698273714 -32.3750737220 -54.1409079085 -22.8566657091 -21.3529419302 -28.0269932880 -57.9342423599 -12.5471060288 -38.4620442789 -26.0802038377 -12.5737109718 -20.6795580739 -45.8641417250 -24.5111482190 -25.7049498700 -38.8280386972 -36.5398569319 -60.5896202867 -36.7435309909 -30.0220830620 -12.1161618453 -33.9536690133 -38.0741129286 -21.2930676557 -52.3015259979 
			 6 	-26.1935218278 -40.2987703957 -22.7262607194 -46.4400791410 -24.9714183800 -42.8169232003 -20.1643828757 -42.0507366898 -28.5747452089 -28.0083421525 -39.7246663356 -33.7310099240 -41.1299189380 -41.2058141143 -31.4724379090 -32.5898674104 -51.6529300289 -46.9049313319 -46.4119755392 -43.6988268200 -37.2384403573 -16.9833793690  -8.3076413036 -63.2286734275 -30.8777808860 -66.6719756905 -23.3385897936 -38.5629891424 -32.1156832532 -23.1696872620 -24.3733026451 -53.4179810055 -50.7382263875 -38.1085601472 -26.2780253977 -35.0339154691 -15.3898601504 -12.9218029282 -18.6959124637  -0.2061907411 -15.7647650478 -29.5004105728 -34.7916830255 -39.7872624659 -24.5649679117 -24.1263068784 -30.1172243499 -24.7387196557 -41.8389155305 -53.0171261697 
			 7 	-44.6412099916 -24.2440743204 -30.4846556077 -22.9802534497 -25.1855236541   6.4740024615 -49.0064679768 -44.6251032282 -26.6652113819 -34.5330264209 -22.8338645565 -48.9777251477 -32.5264693214 -26.9554197480  -4.4768335627 -48.2013530490 -35.5229392545 -13.4433566680 -26.9493673874 -30.3319366441 -44.4819457755 -33.3436594187 -40.0156606348 -25.4298988520 -38.3534667557 -42.3471054699 -36.0714813812 -16.3998944839 -33.0250207251 -31.2015144995 -24.0735012932 -13.9699620884 -19.5230002054 -33.2393722725 -33.9793157942 -26.0962892888 -53.5297509335 -36.7412919298 -42.2568357404 -47.9031138512 -26.2649690423 -40.4808116565 -22.8537661315  -2.7122548488 -38.7776181859 -15.9293456965 -33.4037234957 -42.8140755874 -44.5313846065 -43.3333014465 
			 8 	-64.1961959877 -34.2926707327 -37.6341347877 -20.2463278727 -20.4596635204 -42.5219696675  12.2795310680 -52.2101973800 -52.4441987858 -15.6965556156 -22.8543997149 -28.1606632626 -38.0529596067  -8.5333985692 -43.4756367445 -55.3169322141 -50.1841801921 -28.9646513274 -10.6943166974 -35.2859918939 -38.2091313371 -31.9238947283 -26.6475013639 -29.9158048269 -35.1745487745 -22.3881671624 -22.3738533992 -35.6888841767 -22.5246886742 -33.9983168331 -60.6297284764 -32.2554423764 -18.7438414275 -20.6611557570 -28.4913583509 -43.3294785075 -58.2643732159 -12.6635864557 -30.9755940933 -22.3139173054 -20.8166902261 -37.4016723286 -13.5842584474 -12.3691553239 -39.7818818052 -38.5279753435 -35.9520688567 -26.4503837574 -28.0305511908 -29.8672689268 
			 9 	-45.3005698015 -21.2550359055 -16.2749924830 -30.9110433977 -21.6955220889 -31.7439757512 -29.4450907505 -15.2429099634 -27.1610574505 -47.9839750807 -21.5251161686 -20.6854662907 -43.9629697378 -40.6629046067 -41.5030576554 -44.1247352814 -12.6922014586 -24.8944396715 -46.2959805498 -20.3844634833 -39.0743500405  -7.2108499350 -33.5155906278 -29.2958536980 -19.9438745697 -47.5360422941 -56.7528801942 -28.9868427564 -38.9379133632 -51.3346762813 -26.2146834813 -37.8472677514 -35.7188019777 -29.0411434678 -18.4848693810 -32.9209903173 -37.6072352802 -34.7114542346 -40.1803928069  -4.7374555287  -9.2785638631  -5.7158690537 -24.0385067595 -33.0086236498 -26.5093641049 -35.5690328892 -19.4257364735 -42.6375015468 -45.6267331842 -32.1110443016 
			 10 	-34.2998345158 -32.0842875364 -49.6417443633  -0.3589067412 -24.7391033508 -37.6033459843 -11.5623758224 -50.6263587266 -45.9360589633 -26.5471577424 -37.6867621041 -52.6177139058 -12.6295288508 -25.4065561520 -28.0013472610 -38.4713617275 -47.4853101874 -16.4045918206 -45.9292958575 -30.8027627393 -25.9366345724 -28.2016599889 -19.0063133862 -36.3799124554 -25.0417908330 -32.0502490401 -32.2598577329 -38.3453715322 -39.4753808878 -29.3274766522 -40.9567918108 -25.8688002211 -51.2447909390 -52.5310377244 -40.0016450566 -22.3489503810 -25.9210858314 -28.6449429527 -31.7624626692 -50.6253623908 -15.3940769573 -36.7322995102 -42.4903249192 -27.2331450246 -29.1696009987 -35.6270737182 -38.5246449588 -35.0225869379 -29.7573533174  -2.4751127370] ;
	 Beta_FEE= Beta_FEE[:,2:end];
	 #Variables choice model
	  AT_FSP = 10;
	  TD_FSP = 10;
	  AT_PSP = 10;
	  TD_PSP = 10;
	  AT_PUP = 5;
	  TD_PUP = 10;
	  Origin =[
	 	1	0
	 	2	1
	 	3	1
	 	4	0
	 	5	0
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	0];
	  Origin= Origin[:,2:end]	;
	  Age_veh =[
	  	1	0
	 	2	0
	 	3	0
	 	4	1
	 	5	0
	 	6	0
	 	7	1
	 	8	0
	 	9	0
	 	10	0];
	  Age_veh= Age_veh[:,2:end];
	  Low_inc =[
	  	1	1
	 	2	1
	 	3	1
	 	4	1
	 	5	0
	 	6	1
	 	7	1
	 	8	1
	 	9	1
	 	10	0];
	  Low_inc= Low_inc[:,2:end];
	  Res =[
	 	1	1
	 	2	1
	 	3	1
	 	4	0
	 	5	1
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	1];
	   Res= Res[:,2:end];

	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	for n=1:N
	  	for r=1:R
	  		q_parameter[1,n,r] = Beta_AT[n,r] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
	  		q_parameter[2,n,r] = ASC_PSP + Beta_AT[n,r] * AT_PSP +  Beta_TD * TD_PSP;
			q_parameter[3,n,r] = ASC_PUP + Beta_AT[n,r] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n];	
	  	end 
	end
	  
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R);
	for n=1:N
	  	for r=1:R
	  		Beta_parameter[1,n,r]=0;
	  		Beta_parameter[2,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
	  		Beta_parameter[3,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	 	end 
	end
 	# 	 file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/MixedLogitn10_r50.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
    return Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p;
end

function Mixed_Logit_n10_r100()
	NUM_POINTS=2+1; #alternatives
	N=10; #customers
	R=100;
	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];
	#Parameters choice model
	#Parameters choice model
	 ASC_PSP = 32;
	 ASC_PUP = 34;
	 Beta_TD = -0.612;
	 Beta_Origin = -5.762;
	 Beta_Age_Veh = 4.037;
	 Beta_FEE_INC_PSP = -10.995;
	 Beta_FEE_RES_PSP = -11.440;
	 Beta_FEE_INC_PUP = -13.729;
	 Beta_FEE_RES_PUP = -10.668;

	 #AT coefficient (10 customers x R_AT draws)
	 address="C:\\Users\\20176914\\surfdrive\\scientific\\University\\Codes\\Julia\\NonlinearProgramming_Virginie\\random beta\\beta_AT_10_100_given.mat";
	 input=matread(address)
	 Beta_AT=input["Beta_AT"];
	 
	 #FEE coefficient (10 customers x R_FEE draws)
	 address="C:\\Users\\20176914\\surfdrive\\scientific\\University\\Codes\\Julia\\NonlinearProgramming_Virginie\\random beta\\beta_FEE_10_100_given.mat";
	 input=matread(address)
	 Beta_FEE=input["Beta_FEE"];


	
	
	 #Variables choice model
	  AT_FSP = 10;
	  TD_FSP = 10;
	  AT_PSP = 10;
	  TD_PSP = 10;
	  AT_PUP = 5;
	  TD_PUP = 10;
	  Origin =[
	 	1	0
	 	2	1
	 	3	1
	 	4	0
	 	5	0
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	0];
	  Origin= Origin[:,2:end]	;
	  Age_veh =[
	  	1	0
	 	2	0
	 	3	0
	 	4	1
	 	5	0
	 	6	0
	 	7	1
	 	8	0
	 	9	0
	 	10	0];
	  Age_veh= Age_veh[:,2:end];
	  Low_inc =[
	  	1	1
	 	2	1
	 	3	1
	 	4	1
	 	5	0
	 	6	1
	 	7	1
	 	8	1
	 	9	1
	 	10	0];
	  Low_inc= Low_inc[:,2:end];
	  Res =[
	 	1	1
	 	2	1
	 	3	1
	 	4	0
	 	5	1
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	1];
	   Res= Res[:,2:end];

	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	for n=1:N
	  	for r=1:R
	  		q_parameter[1,n,r] = Beta_AT[n,r] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
	  		q_parameter[2,n,r] = ASC_PSP + Beta_AT[n,r] * AT_PSP +  Beta_TD * TD_PSP;
			q_parameter[3,n,r] = ASC_PUP + Beta_AT[n,r] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n];	
	  	end 
	end
	  
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R);
	for n=1:N
	  	for r=1:R
	  		Beta_parameter[1,n,r]=0;
	  		Beta_parameter[2,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
	  		Beta_parameter[3,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	 	end 
	end
 	# file = matopen("C:/Users/20176914/surfdrive/scientific/Codes/Julia/NonlinearProgramming_Virginie/MixedLogitn10_r100.mat", "w")
	 # write(file, "Beta_parameter",Beta_parameter)
	 # write(file, "q_parameter", q_parameter)
	 # close(file)
    return Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p;
end
function Mixed_Logit_n10_r200()
	NUM_POINTS=2+1; #alternatives
	N=10; #customers
	R=200;
	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];
	#Parameters choice model
	#Parameters choice model
	 ASC_PSP = 32;
	 ASC_PUP = 34;
	 Beta_TD = -0.612;
	 Beta_Origin = -5.762;
	 Beta_Age_Veh = 4.037;
	 Beta_FEE_INC_PSP = -10.995;
	 Beta_FEE_RES_PSP = -11.440;
	 Beta_FEE_INC_PUP = -13.729;
	 Beta_FEE_RES_PUP = -10.668;

	 #AT coefficient (10 customers x R_AT draws)
	 address="C:\\Users\\20176914\\surfdrive\\scientific\\University\\Codes\\Julia\\NonlinearProgramming_Virginie\\random beta\\beta_AT_10_200_given.mat";
	 input=matread(address)
	 Beta_AT=input["Beta_AT"];
	 
	 #FEE coefficient (10 customers x R_FEE draws)
	 address="C:\\Users\\20176914\\surfdrive\\scientific\\University\\Codes\\Julia\\NonlinearProgramming_Virginie\\random beta\\beta_FEE_10_200_given.mat";
	 input=matread(address)
	 Beta_FEE=input["Beta_FEE"];


	
	
	 #Variables choice model
	  AT_FSP = 10;
	  TD_FSP = 10;
	  AT_PSP = 10;
	  TD_PSP = 10;
	  AT_PUP = 5;
	  TD_PUP = 10;
	  Origin =[
	 	1	0
	 	2	1
	 	3	1
	 	4	0
	 	5	0
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	0];
	  Origin= Origin[:,2:end]	;
	  Age_veh =[
	  	1	0
	 	2	0
	 	3	0
	 	4	1
	 	5	0
	 	6	0
	 	7	1
	 	8	0
	 	9	0
	 	10	0];
	  Age_veh= Age_veh[:,2:end];
	  Low_inc =[
	  	1	1
	 	2	1
	 	3	1
	 	4	1
	 	5	0
	 	6	1
	 	7	1
	 	8	1
	 	9	1
	 	10	0];
	  Low_inc= Low_inc[:,2:end];
	  Res =[
	 	1	1
	 	2	1
	 	3	1
	 	4	0
	 	5	1
	 	6	1
	 	7	0
	 	8	0
	 	9	1
	 	10	1];
	   Res= Res[:,2:end];

	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	for n=1:N
	  	for r=1:R
	  		q_parameter[1,n,r] = Beta_AT[n,r] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
	  		q_parameter[2,n,r] = ASC_PSP + Beta_AT[n,r] * AT_PSP +  Beta_TD * TD_PSP;
			q_parameter[3,n,r] = ASC_PUP + Beta_AT[n,r] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n];	
	  	end 
	end
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R);
	for n=1:N
	  	for r=1:R
	  		Beta_parameter[1,n,r] = 0;
	  		Beta_parameter[2,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
	  		Beta_parameter[3,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	 	end 
	end
 
    return Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p;
end
function Mixed_Logit_n50_random(R_AT)
	NUM_POINTS=2+1; #alternatives
	N=50; #customers
	R=R_AT;
	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];
	#Parameters choice model
	#Parameters choice model
	  ASC_PSP = 32;
	 ASC_PUP = 34;
	 Beta_TD = -0.612;
	 Beta_Origin = -5.762;
	 Beta_Age_Veh = 4.037;
	 Beta_FEE_INC_PSP = -10.995;
	 Beta_FEE_RES_PSP = -11.440;
	 Beta_FEE_INC_PUP = -13.729;
	 Beta_FEE_RES_PUP = -10.668;

	 #AT coefficient (10 customers x R_AT draws)
	address="C:\\Users\\20176914\\surfdrive\\scientific\\University\\Codes\\Julia\\NonlinearProgramming_Virginie\\random beta\\Beta_N"*string(N)* "R_"*string(R_AT)*".mat";
	input=matread(address)
	Beta_AT=input["Beta_AT"];#AT coefficient (10 customers x R_AT draws)
	Beta_FEE=input["Beta_FEE"];#FEE coefficient (10 customers x R_FEE draws)
	#Variables choice model
	  AT_FSP = 10;
	  TD_FSP = 10;
	  AT_PSP = 10;
	  TD_PSP = 10;
	  AT_PUP = 5;
	  TD_PUP = 10;
	  Origin =[	1	0
	2	1
	3	1
	4	0
	5	0
	6	1
	7	0
	8	0
	9	1
	10	0
	11	1
	12	1
	13	0
	14	0
	15	0
	16	1
	17	0
	18	0
	19	0
	20	1
	21	0
	22	0
	23	0
	24	0
	25	0
	26	0
	27	0
	28	0
	29	1
	30	1
	31	0
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	1
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Origin=Origin[:,2];
	  Age_veh =[	1	0
	2	0
	3	0
	4	1
	5	0
	6	0
	7	1
	8	0
	9	0
	10	0
	11	0
	12	0
	13	0
	14	1
	15	0
	16	1
	17	1
	18	0
	19	1
	20	1
	21	1
	22	1
	23	0
	24	0
	25	1
	26	0
	27	0
	28	0
	29	0
	30	0
	31	1
	32	0
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	1
	40	0
	41	1
	42	0
	43	1
	44	1
	45	0
	46	1
	47	0
	48	0
	49	0
	50	0];
	Age_veh=Age_veh[:,2];
	  Low_inc =[	1	1
	2	1
	3	1
	4	1
	5	1
	6	0
	7	1
	8	1
	9	1
	10	1
	11	1
	12	0
	13	1
	14	0
	15	1
	16	1
	17	1
	18	1
	19	1
	20	0
	21	0
	22	1
	23	1
	24	1
	25	0
	26	1
	27	1
	28	1
	29	1
	30	1
	31	1
	32	1
	33	0
	34	1
	35	0
	36	1
	37	0
	38	1
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	0
	49	0
	50	0];
	Low_inc=Low_inc[:,2];
	  Res =[	1	1
	2	1
	3	1
	4	0
	5	1
	6	1
	7	0
	8	0
	9	1
	10	1
	11	1
	12	1
	13	0
	14	0
	15	1
	16	1
	17	0
	18	0
	19	0
	20	1
	21	1
	22	0
	23	0
	24	0
	25	0
	26	1
	27	0
	28	0
	29	1
	30	1
	31	1
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Res=Res[:,2];
	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	for n=1:N
	  	for r=1:R
	  		q_parameter[1,n,r] = Beta_AT[n,r] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
	  		q_parameter[2,n,r] = ASC_PSP + Beta_AT[n,r] * AT_PSP +  Beta_TD * TD_PSP;
			q_parameter[3,n,r] = ASC_PUP + Beta_AT[n,r] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n];	
	  	end 
	end
	  
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R);
	for n=1:N
	  	for r=1:R
	  		Beta_parameter[1,n,r]=0;
	  		Beta_parameter[2,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
	  		Beta_parameter[3,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	 	end 
	end
 
    return Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p;
end
function Mixed_Logit_n50_r50()
	NUM_POINTS=2+1; #alternatives
	N=50; #customers
	R=50;
	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];
	#Parameters choice model
	#Parameters choice model
	  ASC_PSP = 32;
	 ASC_PUP = 34;
	 Beta_TD = -0.612;
	 Beta_Origin = -5.762;
	 Beta_Age_Veh = 4.037;
	 Beta_FEE_INC_PSP = -10.995;
	 Beta_FEE_RES_PSP = -11.440;
	 Beta_FEE_INC_PUP = -13.729;
	 Beta_FEE_RES_PUP = -10.668;

	 Beta_AT=[1 	  1.9798340694   0.2127092528  -1.5283686950  -0.2437148182  -0.5688439127  -0.7039031828  -1.5282484580  -2.1526582361   0.6118160383  -1.8585600676  -1.2044311044  -1.0801302756  -0.0860983468   0.2713654095  -0.2464256231  -1.6440009782   0.5993743658  -1.9624511986  -1.4161262778   0.1401277723  -0.0764609209   0.6138903628  -0.7931700693   0.8504117131  -1.1018504591  -2.1978469613  -0.3200657469  -2.3834896910  -1.1618047506  -2.5633899588  -1.0805035713  -0.7958813331  -0.6287157907   0.4834891048  -1.0880818650  -0.7678236234  -1.0367067121  -1.9408237095  -0.1024375377  -0.9612420709  -2.2259737090  -1.6091899155  -1.0057136917  -0.8437599521  -1.1671834188  -2.6709328946  -0.3363695728  -0.5565002827  -0.0871501948  -0.1922347056 
			 2 	 -2.3556698088  -0.8088819586  -1.0109090619  -1.5611762656  -1.1822996918  -0.6946579047  -0.1502042913  -0.6251406593  -2.8986926994  -1.6111013017  -1.1372894114  -0.5478099449  -0.4208527947  -2.2366052521  -0.1927872228  -0.5537387449  -1.0850911139   0.4630477134  -1.7989385526  -0.4495877617  -0.4261284597  -0.1429639224  -0.2796075580  -1.3939177382  -1.0011882328  -0.2519902732  -2.1167414240   0.1790053800  -0.6544813133  -0.7624954488  -2.3779421487  -0.5424297757   0.7373980290  -2.1151658440   0.6358219601   0.4248594669   0.6674606969  -0.4806544791   1.1591213935   0.6153218963   0.7808668836  -1.2583409245   0.7610417396   0.5521237054   0.7080541842  -2.0702313179  -1.8622042688  -1.7785298203  -1.5373455124   1.0974890264 
			 3 	 -1.4297323554  -0.9257545191  -1.0510640753  -0.6602896589   0.0645136700  -2.0092843297  -2.7743725719  -0.5581697711   0.2588200115  -0.3658589364  -0.0376397878  -2.2465281542  -0.5644174629  -1.2267004821  -2.1329637213  -0.8959483797  -2.3232406200  -0.9847446919  -1.2998139116  -1.2675468334  -0.1874122227  -0.9472618587  -0.1218955861   0.6820822511  -2.1956421743   0.0618254362  -0.6635662108  -0.5530007642   0.0048979518  -0.1078325449  -0.3681644651   0.0970192193  -1.1442471625  -2.3867993808  -0.6631360653  -0.9244828278  -0.0857905823   1.1400730904  -1.0137641367  -3.2341634500  -0.0005586133  -0.1301800165  -2.2152869546  -0.9069124541  -2.1685795769  -3.5017435082  -1.2617766664   0.8021249716   0.3678734331  -2.1308683968 
			 4 	  1.9981389437   1.4150966325  -2.2058293344  -1.5924552133  -2.0201704419  -1.0369022158  -1.8630776146   0.2441489883   0.0630885571  -0.9604412605  -2.5365314998  -1.3303570533  -3.8421446878  -1.6753737418  -0.9342752549   1.0660052523  -0.8806305888  -1.6768613706  -0.3920940688  -0.4015801107  -1.9054765762   1.2444646710  -0.9502755143  -0.1423542151  -3.6775544909   0.6940304084   0.2583757223  -1.3041316941  -1.3282160720  -0.2702587828   0.2006901598   0.5859110399   0.9648430533   0.0363825826  -3.5378567898  -0.3783452233  -0.7536627868  -0.7188630206  -0.3834414910   0.2052568572  -1.0225970721  -0.3892032507  -0.1812398248  -1.4759106882  -0.7557912228  -1.0162218260  -1.1146333616  -0.8173893788  -2.3280362876  -0.7454204880 
			 5 	  0.7534053857  -2.6060120824  -0.6333247680   1.6864766752  -0.8742041026  -2.0604162715  -0.8732838081  -0.3567654985   1.4136909052  -1.2365632704  -0.2983921588  -1.5965642294  -0.4921279652  -0.9388989130  -0.4598498244  -0.4517631431  -1.1998825697   0.1945942065  -0.5742653564   0.3309232181  -1.0477395653   0.3316922060  -0.4169869599  -0.4062763948   1.1233722017  -3.0294612644  -0.2753167071   0.2413328861  -0.2688036618  -0.9195853017  -0.5053891187  -1.3301108834  -1.2633045765  -1.3192301090  -1.1600808904  -1.8282816971  -0.6537553740  -2.0711960635  -0.4927597982  -1.2285998278  -1.0091191284   0.4167663846  -1.9620897654  -0.5407823197  -0.0675107038  -1.1781396556  -0.1049993066  -1.0986953119   0.2477564518   0.5954284955 
			 6 	 -0.1212409610  -0.9319959277  -0.3549542462  -1.2438588341  -2.3554430062  -0.4581197274  -0.9848107249   0.2099882947  -2.1170076849  -2.4437835953   0.8072475433  -0.4965463187   1.2401576212  -0.3296888141  -2.8735196562  -0.2536859567  -1.2603232590   1.2780567621   0.3851634597  -0.8018052154  -0.8909785102  -0.8991004444  -0.4032619252  -1.4403658376  -1.5616437059  -1.7993883299  -1.2503420074  -2.0646861696  -0.0780555473  -0.7661309100  -1.6172639450  -0.7171828167  -0.5588802956  -0.9622852929   0.0870661791  -1.3070313542   0.3559910870  -0.5693482957  -1.9160438149  -1.3355031612  -0.5031732741  -1.2725968290  -3.5385839439   1.1915990808   0.9721690615   1.8343661533  -1.5856710193  -0.6886220310  -2.2057514065  -0.2485781034 
			 7 	 -0.4167608084  -0.8204426637  -0.7941669916  -1.4612313886  -0.5901444241  -0.6286471900  -0.6660912972  -0.7876490072   0.6825240571  -2.5288573554  -1.1536615393  -0.7099594753  -1.0281400298   0.3566728324  -0.7829543043  -0.5824941425  -1.0671283676  -0.9031339682  -2.0913895639   1.3404458016  -1.1597969282  -0.5229377654   0.6068984125  -1.1967871218  -0.5705876268   0.9676872654  -0.6420519096   0.4647196685  -1.1572705742  -1.2793251457  -1.9558641774  -0.3528063514  -1.3526232098  -0.6276237428  -0.3728328515  -1.8053260749  -1.6447138450   0.8159762667   0.0924375829  -1.7521267552  -1.7022711327  -0.8801580404  -0.5807851693  -1.5782448701  -2.1459900402  -1.7115292891  -2.4048639171  -4.0893091493  -0.5862319412  -1.5599243338 
			 8 	 -1.7654091144  -1.5464937656  -0.8725640937  -0.3187656469  -0.6593312749   0.6986996831  -1.4685605047  -1.6484858785   0.5740235662  -0.9966774540  -1.1954901739  -1.9714460035  -1.8226920545   1.4817465020  -2.4455665819  -0.6824905558  -1.8090592843  -2.5672080526  -1.1369987827   1.3645340927  -2.2773677408  -1.5415003661  -0.6817458860  -1.5153880505   1.1776343599   0.5193750906   0.7824140968  -0.3109587461   1.1847460731  -0.0141528640  -1.2094530376  -0.4508826615   1.9814229641   0.7327456908  -0.3310267646  -1.3222636264  -0.6950400336   0.3195732233  -0.5136485261  -0.5886869605  -1.0295685414  -0.9723404263  -2.5432527829  -1.0871614672   0.4160785321  -0.9947922440  -1.6291384359  -0.5153096275  -1.7841351723   0.7109779792 
			 9 	 -0.7381401572  -0.2380163482  -2.0466452896   0.6734033588  -0.9097096681  -2.7054698979   0.5985132663   0.5933161469   0.2940842033  -0.8096150286  -1.7422849121   0.3126045682  -1.0785713353  -1.5520987587  -3.6058773061  -1.2840731117  -2.1756764506   0.4736259979   0.0683564619  -1.2540874299  -1.0283157110  -2.3983516599  -1.5537870379  -0.4445311529  -1.9656181434  -3.0753895205  -2.5568604352  -0.1732735087   0.8776422703   0.0360887672  -2.8774556081  -0.1111601070  -0.8144850743  -0.2074891086   1.3249457559  -1.0890372627  -0.4006916986   0.7509193953  -2.7112508224  -1.5554815393  -0.5387611033  -1.3696736144  -1.5551560412  -2.3729173651  -0.4095149324   1.4845510673  -0.1739032934  -1.0966569077  -0.4214882070   0.1035249234 
			 10 	 -0.6570030203  -1.8114245261  -3.2435799807   0.0290398020  -2.9084056303  -0.5451223698   1.1838278553  -0.1608428317  -0.3104115836  -1.2644713463  -1.3652121562  -1.4074969567   0.5032117223  -1.6683406507  -1.7783062004  -0.6769713063  -1.7127772658  -0.1566364402   0.3113288063  -0.6875535531  -2.4351416305  -0.7465249142  -1.7239654425   0.8592401247   0.9042202359  -0.0529651526   0.1515934447   0.2698748917   0.0552395517   0.4836946003  -2.9314592223  -0.2259842687   0.5909445735   0.3247203446  -1.3193356975  -0.2693240160  -2.5594569304  -1.2356361997  -1.1913512034  -0.1861646804   0.0047060163  -0.3886472890   0.7656516146   0.0993720539  -1.3953556381  -1.2902933766   1.3952329967  -1.3330861224   0.2080794741  -1.9189839456 
			 11 	 -0.1893213055   0.6441566216  -0.2405910677  -2.6036575905  -0.2714023589  -1.5063786639   0.1484888206   0.3532866860  -1.9413841976   0.1667840825   0.5469983648  -1.2861018308  -1.0318690318  -0.3223691353  -1.9768607004  -0.4280940774  -3.2914739343  -0.5415387363  -1.2692096126   1.4324013439   0.2651168593  -0.5712198226  -0.8622365932  -0.0159627746  -0.7104765518   1.4854932137   1.3539267198   0.4588360850  -0.9962837943  -1.7836030936  -1.3597782652  -1.5135360123  -0.5088671064  -0.1625602049  -0.3610850791  -1.5560865601  -0.5999730143  -0.7524157956   0.1844640455  -2.3730129570   1.4190544303  -1.1170481534   0.7340355730  -1.0005610259   0.6951258883   0.1623125614  -0.8577232942   0.1845098727  -1.4732010380  -0.6707403209 
			 12 	 -2.3774975910  -0.8663276572  -0.5716821608  -0.3574457104  -1.0519950268   0.1738743869  -0.5491520830  -0.7979799206  -3.5814134842  -1.8385372189   0.0004788027  -1.3392309784  -0.9004387249  -1.2137368784  -1.1024207525  -2.0217856505  -1.2088860353  -2.4884220860   2.0947133952  -0.9469244764  -1.4828148700  -0.2301533943  -0.6548957827   0.3316508951  -0.3289042885   0.6215140561  -1.1899801035   0.6289757227  -0.8041160803  -2.8420202800   0.5038095157  -0.4724971284  -0.3788730333   0.7004209099   1.4545693743  -2.2064192725  -0.3081789126  -3.0442280552   1.3580356228  -2.1467088274   0.6018286185   1.1169472175  -1.9594026004  -0.3780521699  -1.1260691556   1.9528982663  -1.1564398430  -1.0856956158  -1.4302168343  -1.1195950693 
			 13 	 -0.1615351498  -0.4398664361  -0.3622205464   0.5216938665  -1.6664777991  -0.9786967173  -1.3714178460  -1.1325326880  -0.7515466596  -0.3020388765   0.3582150042  -0.4325006766  -1.5643117702  -0.1934420007  -2.0741221231  -0.6905110238  -3.3883124905   0.2316196964  -1.0688005880  -0.5263802822  -0.1809002744  -0.1270821529   0.5808401837  -0.4076446224   0.2343030143   1.7466984997  -3.8025514635   1.8418491669  -0.5286275954  -1.4642436432  -2.2105888103  -1.2591412080  -1.8799743908  -1.0929516998  -0.6785606894  -0.1792440492  -1.9312524304   1.6977028863  -0.9487168821  -1.5293815269  -0.5296886369  -0.5517154104   0.1993910766  -1.2414957982  -1.8688855107  -0.7636457097  -1.3646973531  -0.5408997245  -1.7347470673   0.0622325457 
			 14 	 -1.2423859625  -3.9845483929  -1.2518215266  -0.6290936589  -0.5569068693   0.2536176260  -0.8907337380  -0.2600251849  -0.6434955588   0.4402166740   0.1787194182  -0.8378310460   0.2539138976   0.5970602693   0.2397773007  -1.6727391146  -1.0960472768  -0.8416543154  -1.0225914277  -0.1284166437  -1.7019532181  -2.2718409003   1.5267918390   0.2300159161  -1.9630303275  -1.0686157919  -0.0676159165  -0.8380530391  -0.1509441894  -1.2128922865  -1.1272227288  -0.7149477512   0.1854096446   1.7600203590  -0.9818591301  -1.3309602814  -0.2282879006   0.0959750682   1.1668728388  -0.1135561302   0.3048934730  -0.0996168029  -1.2137743843  -0.4614664282   0.4874050774  -2.5209975231  -2.2983982963  -0.6343235895  -0.7027625137  -0.8606020268 
			 15 	 -0.7294111312  -1.0526116603  -0.5024547323  -1.7939281734   1.4323685940  -0.9954012506  -1.3060293825  -1.6436254217  -2.1831735637  -2.3512823656  -0.7361950001  -2.3404201513  -0.8581640411   0.3209066947  -0.5257552255  -1.2892556287  -0.1975214498   0.2529545267  -0.2729919740  -2.2512611361  -0.3440395998  -2.4694686369  -1.8650106517  -1.0465681739  -0.5646059749  -0.5607185799  -2.2132156595  -2.2131183734   0.5204608096  -0.3345303607  -0.0974205804  -3.3802032563   0.8755660833   0.8201875378  -0.5657464407   0.4376895533  -2.3030607620   0.9103866341   2.0916699832  -2.1020751300   0.2881760353  -0.1560886875  -0.6737269836  -0.7451924802   0.3549038780  -0.0411042346  -0.3515522579  -0.1810505071  -0.5425443120  -0.1236800014 
			 16 	 -0.1426520610  -0.8634068783  -2.2646188789  -0.1578496330  -2.2184665505  -0.5697209002  -1.6099647775   1.0631602091  -1.7253707496  -0.5037057760  -0.8237182396  -2.4769772528  -1.9676523319  -0.3603310681  -1.2191407583  -0.9023202452  -0.3491728368  -0.7042772003  -0.3301005817   1.4459755546  -0.3689992359  -0.8637963888  -0.1560552395  -2.3805599017  -1.6093341032  -0.5942901453  -0.4825085274  -2.3545262945  -1.0224483536  -2.2687768019  -1.4717275069   0.4999308710  -0.2040596016   0.2096929359  -2.5025274041  -1.1819975848  -1.3945822431  -0.7199804324  -0.9227776517  -0.6441017531  -0.3796019532  -1.5075595168   0.2340669701   0.2543061124  -0.8333866097  -2.5033908616  -2.0938425807   0.9199792518   0.0550666722  -1.9569343181 
			 17 	 -1.4734545042  -1.8772218860  -0.5261458365  -0.6743991275  -1.4698153497  -0.3711311660  -0.5007711003  -1.6575957382   0.7700091690   0.2603340289  -1.1312800560  -1.0505798560  -0.3484787035   1.2215614432  -0.2487692440   1.2651420071  -0.5370640315   1.2137901002  -0.7291150241  -1.1564483696  -0.8432253300  -0.4122724043  -1.3915432715   1.2785272728  -1.7189897708  -1.1849097813   0.4151424081  -2.9463080143  -0.5109893488  -1.3107791565   1.1657833511  -0.8825357379   0.3012849122   0.1951574833  -1.6872104981  -2.2557538871  -0.7657223828   0.4754459386   0.7630142273  -0.8449558260   2.8533234096   0.7539287578  -0.3491408117  -0.9386206450  -1.2447258791  -1.5737008001  -0.9955485220  -1.0325619418  -1.5682194929  -0.6646911075 
			 18 	 -2.1895756893   0.1159467410  -1.8249012954  -3.4560506459   1.4867092531  -0.2123888153  -0.2552175067   1.3439772659  -0.8055015270  -2.7515122326  -1.4722238572  -0.6794681276  -0.1435232623  -1.5087949373   0.7476672060  -0.0533685727  -0.0596019617  -0.9419098392   0.6951262170  -0.0988580126  -1.9127830843  -0.9812547770  -2.4147368616  -1.7084625710   1.3512904246  -0.6139908551  -1.9167714525   0.7508599615  -1.6233376325   0.8828829375   0.0745065045   0.8483069931   0.4806445306  -0.9601306568  -2.0185907131  -2.1881430381  -0.0096893198  -0.7285338524   0.0258182383  -1.5245206836  -2.1639888204   0.7942119791  -0.2105629914  -1.0750084227  -1.7498260303  -2.5661465070  -0.6237065598  -0.1543806133   0.4212680747   1.0182687349 
			 19 	 -1.3013114095  -2.8952038558  -0.5340059001  -0.6536891555   0.5022710813  -0.8096465316   2.1631754643  -0.2904807627  -1.9216084616  -0.1801003726  -0.9190370459  -1.0056799882  -1.5276492223  -1.0293456509   0.0899208554   0.1008855019  -0.6949109286   0.3268815672  -1.7546355530  -0.2486131475  -0.2474557589  -0.4283351062  -1.8348964349   0.2162202896  -1.0763864974  -1.3127526577   0.0963544673   1.1054128853  -1.0846780018  -0.5495239754  -0.7394052435  -0.4883981713  -1.7601045994  -0.3643602063  -0.3440398030  -0.6293203040   0.2124751610  -0.1355915968   0.0048943259  -0.2296403277  -0.1830751434   0.4089206193   0.4123618886  -0.3155905057  -1.1881184562  -0.6210388951  -0.5786838492  -0.7059705800   0.9078499883  -0.0246209521 
			 20 	 -1.6232077243   0.3596873195   1.1026944321   0.2216024110  -2.3664717319   2.5137200937  -3.2234302144  -0.9421174694  -0.6625462527  -1.1036485935  -1.9872906036  -2.3151713053  -0.5965144975  -0.7963739417  -0.4363495585  -0.0825009284   0.5581757073  -2.5227205928  -1.5705161853  -3.0883377578   0.6106932664  -0.5456387273  -1.3547307044   0.8787277919  -0.9994140791  -0.1576370896  -0.3519296175  -4.4279852313  -1.0251937443   1.1907980286  -1.5948782239  -1.4547971673   0.3720954284  -0.1517603968  -0.6624845913   1.1818584125  -1.2227226198  -2.1358269350   0.4194752670   0.1882463377  -1.1624875975   0.4185005331  -1.7594817679   0.4889652777  -0.4833267223  -1.4733103358   0.2865031585   0.1077654182  -1.3825887615  -1.0218873166 
			 21 	 -1.0246909501   0.3987268944  -1.1607478536  -0.9489463967  -0.4816448774   1.5297552376  -1.7411323496  -0.0259605704  -1.3732682546   0.0531753957  -0.0446652889  -1.4057359838   1.6508127931  -0.5458419208  -0.5493636135  -0.5402339708  -1.3505123542  -0.2754704687   1.7878607492  -2.1905401173  -0.7190169699  -0.2818357852   0.6658587955  -2.8048675557  -0.1279308367  -2.7861338949  -2.2320849264  -0.5274341903   0.4959015831  -0.8131604611  -1.7601715569  -0.2394709477  -0.1840553217   1.3855682650  -0.3646085151  -1.1860614869   1.0269912662  -1.3456170430  -2.7025275671  -0.1294878291  -0.5400531687   0.2031379867  -1.4035866568  -2.6759789307   0.5894265660  -1.1819168778  -0.8012833057  -0.0906532017  -1.3819787331  -2.4203300375 
			 22 	 -0.2404273746  -0.2625909742  -0.1409176510  -0.8720190367  -0.3169890669  -1.0609660001  -0.6329338279  -0.5340860146   0.4457428585  -0.0003041821  -1.3805949199  -0.3231100452  -1.3593562314   1.9699551081  -2.3006961282  -0.5377735506  -2.4828914611  -1.8904222112  -1.9682100857  -0.4449536089  -0.0105247630  -1.3917648117  -0.7327611849  -1.0056072974   0.8304546725  -1.1927685421  -1.1071102797  -1.8285330464   0.8334664281  -2.0181249186   0.0619957352   0.2050374902   0.3670862806   0.4661356944  -1.3286909231  -0.3153164248  -0.8442568963   0.1484216536  -1.2106027539  -0.3092791815  -1.7462149383  -1.3833789122  -1.0419784519  -1.8501949557  -0.1902656127  -1.6628232599  -2.1416975344  -0.0719524078  -1.0277116141  -4.0966174045 
			 23 	 -0.1564538450   0.9189438954  -0.2154891564  -2.7987686076  -1.6335163630  -1.6374201004  -0.2239027357  -1.8519367366  -1.8198127050   0.1660269113  -1.2739000451  -0.3587753238  -0.9019616270  -0.5951896292   0.3670987703  -0.2382538197  -1.4832038799  -0.9863173169   0.0548547900   0.3274678928  -0.1274035849  -0.9971622847  -0.0438983104  -0.6334373996   0.0986730992  -1.7865579251   1.0324527198   2.1225628750   0.3426564512  -2.7005178902  -1.2142595990  -0.0954228149  -0.3867349972  -2.0849810697   0.4701931989  -2.0322257459  -0.0765460913  -0.0496440695  -1.3611141839  -0.7942737918  -1.4176432104  -0.7031114486  -0.0601223686  -0.7905881129  -2.1733379754  -1.8457887214   0.5247253948  -0.1909241966   0.3850838840  -0.4893455603 
			 24 	  0.4122912658  -0.6901169044  -0.2339873654  -0.4377524207  -1.0415048349  -1.5055251281  -1.2820980376  -2.7487829395  -0.1881020545  -0.0296620053  -1.0810371595  -0.6087483811  -2.3685566532  -1.2901963851  -0.0825439457  -1.3238105903   0.4299085828  -0.9310492511  -1.7093923580  -1.2419308165   0.8894790501  -0.3465502553  -2.8919125591  -1.8237656814  -1.7049478668  -1.7574967464   0.7765974472  -1.5736820241  -2.0305486894  -1.6192858476  -1.3409743239  -1.4616169724  -1.7254046586  -2.7569732475  -1.6397128421  -1.2588054365  -1.5327251173  -3.0875340798   2.0068949699  -0.5560385055  -0.9006709202  -0.5253696388  -0.2932927852  -1.3600022421  -1.0079333781  -2.1963229812   0.8018485055  -0.9454466827  -1.1238086333   0.3841540970 
			 25 	 -1.1687179331  -0.1381575638  -1.0408627649  -2.4520040636  -0.5711790971  -1.0793976237  -1.5846358568  -0.3899685380  -1.2333623362  -2.0622018305   0.4158868062  -0.7141181883  -2.5831602947  -1.2269984696   0.5214309855  -0.6453869809  -1.1036720628  -0.3619814713  -2.3000311806  -0.3975156382   0.0011077765  -0.4255906513   0.0364367368  -2.0842429535  -2.4372487466   0.5704618917  -1.0442806858  -0.8193309647   0.3835190988  -0.8467179175  -2.1259753965  -1.5182367145  -1.2306932250  -1.4728171070   0.5659395944  -1.4291609210  -3.4438162940  -1.0006862089  -0.3292003099  -1.3368875884  -2.0342101406  -0.9769193370   0.2009659148  -0.9093021157   0.0574974390  -0.1042058204  -2.3012176803  -3.9139324977  -1.3640626586  -0.2818357852 
			 26 	 -0.6525649677  -0.5053862687   1.2378230391  -1.1321521310  -1.8600227991  -2.3963608920   1.4115465273  -0.5451051880  -1.4758996936   0.5581418872  -1.4156244985  -0.6701546856  -1.3863524373   0.9188376684  -0.3341082264  -0.1424468695  -0.5845904968  -2.2238142595   1.6211862569  -1.8517746887  -0.2137462736  -0.9093914157  -0.8053783993  -0.4000062274   0.8440920035  -0.2538654698  -0.1014554990   0.1164651453  -0.7495021680  -0.6952177851   1.6654908128   0.2882689854  -0.7212323046  -0.2760564262  -0.5641612486  -0.3372143996  -0.5148407304  -2.0637933848  -1.9326882399  -2.2119457880   0.0184530334  -1.7233192206  -1.9015410878  -1.8217413905  -2.2462352747   0.9408972185   0.1397496977   1.1806258305  -1.6784258638  -1.8945017905 
			 27 	 -0.0645592803  -1.8616895439  -0.1182218415   0.6479366405  -0.5914995006  -2.5279522823  -1.1322698546  -0.0837185001  -2.0163141622  -0.4295833298  -2.0017701975  -2.6665972417  -1.5321818240  -0.6549825534   1.7550209533  -0.8264787544  -0.0287205989   0.1038328361  -1.4379559026   0.9728238895   1.1646475837  -2.0988323043  -1.8073703669  -1.3400901790  -0.9722365952   1.0696990823  -1.7092839375  -1.1315062135  -1.6906868187  -1.3400744735  -1.3011692232  -2.8349069317  -1.7191592429   0.0369413618  -2.6771825253  -0.7052436452   0.1740494248   0.2053249176  -2.0940999845  -1.4501027361  -1.4122350060  -1.7395467074  -0.5361667424  -0.7655810718  -0.2482787380  -1.8463553709  -0.1229294824  -2.8861384445   0.4039928108  -1.4342215918 
			 28 	 -0.2353887497  -0.6522168547  -1.4413375689  -2.2735006529   0.5550118455  -1.2779625859  -1.7662479792  -2.1271660458  -1.8070358447  -1.9202444405   0.8908690490  -0.6217014590  -1.5803572444  -1.0220913882   0.5604190807  -0.0332005935  -0.9205663492  -2.7185531968  -1.8564678854   1.0702171065  -2.6276051569  -1.2034691259  -1.1431517071   0.1615750152  -0.1509368088   0.4856127688   0.2093463689  -0.4637422399  -0.4808662857  -2.2501974597  -0.3435126815  -1.7329778371  -0.9317409125  -1.4590032267  -2.8080889605  -3.1258201042  -0.6345045943  -1.1841015286  -1.5088255688  -0.6043038102  -0.4180767657  -1.7050854517  -1.2524724034  -1.1918026159   1.1549163413   1.6976378137   0.1423070642  -2.0819157829  -0.5700605378  -0.5949567416 
			 29 	  1.6961344176  -1.3947801860  -0.6858847872  -1.4471640992  -0.5250713753   0.1775754704  -1.3392255784  -0.1169519138   0.2541244968  -2.7854807899  -0.4468427994   0.2622164292  -1.8735119415  -0.0582392540  -0.8142769469   0.2270632875   0.7742704152   1.7308402222  -0.2906958230   0.2162751998   0.2773268711  -0.9303316023  -3.1359447428   1.4903258162  -1.7555070648   0.1574463240  -2.2674254002   0.1013961723   0.6225830958  -1.1488301893  -0.0461189313  -0.9309741014  -0.1334786888  -0.6053063956  -3.7762379839  -0.2313597388  -0.9116099276  -1.2355391502  -0.2052447323   1.1581393257   0.3511488080  -0.3818548271  -0.1789743127   0.5714983451   0.2441014209  -0.7048071383  -1.1520859753  -2.5976505835  -1.8140844635  -0.5480895129 
			 30 	  0.7100458170  -0.8107355667  -1.6359527061   0.3353737390   0.5971007954  -1.6634861135  -1.1727901184  -0.5889812882  -1.3820483224  -0.0916128660  -0.4467244043  -1.6210289849  -2.5924479352  -2.0184409377  -0.0913079196  -0.9858816317  -0.5261010911  -1.3945064076  -1.1910739375  -0.9733828589  -1.3438503686   0.5651899038  -3.0204641221  -1.5968074365   0.5712235337  -0.1377101386  -0.6658463750  -0.5861030831  -0.6762387857  -1.2939013658   0.0414520051  -0.3209029811  -1.3335105710  -0.5805848057  -0.0790209590  -1.2639557696   0.3524840478  -0.7169297342   0.3527051372  -1.1461658036  -0.0267034536  -2.7779062724  -0.2344864916  -0.9580345303  -0.0748554332  -1.0552112085  -1.2282229182   1.9430693304  -0.3651659710   0.7296290847 
			 31 	 -1.1204289396  -1.6346078277  -2.0897705217   1.3771957551  -0.2504352842  -0.7115637930  -0.4686502527   0.3910606689  -1.6872949104  -1.1346647084  -0.8408968865  -0.6110039702  -0.3654750816   0.0003104585  -1.7256203213  -1.5944387527   0.3012985274  -0.4695609522  -2.6049568358  -1.2960564397  -1.1598386817  -0.4078540386  -2.9504705234  -1.3939664834  -1.5089785237   0.1605055335  -0.7998603991  -0.2053084766   0.0950752982  -1.7125842304   0.5630208664   0.0696173215  -0.6172679476  -1.8760200002  -0.8194090463  -0.8406494007  -0.2262817333  -0.5539136260  -2.3301726821  -0.0298120055  -1.8556408076  -0.9670868503  -0.5869136418  -0.2880641481  -1.1813070438   0.7924592932  -0.7830082045  -1.3591194351   0.9654089797   0.0448008493 
			 32 	  0.7103509197   0.1633287684  -0.5418007185  -0.6137708303  -0.4275239216   0.3655984197  -1.2194691174  -0.8556659647   0.0549060253  -0.8369957801  -0.4975811703  -0.5714672721  -0.8068598457   0.2223382920  -0.6983349555   0.0393649078  -1.7099164482  -2.3971186823   1.0782360539   1.2401411698  -1.2884983182  -0.9914803620  -2.3306330990   0.6369289150  -0.9002172438   0.6280166743  -2.1478290165   0.0595022539  -1.5366206303  -2.0104176559  -0.5475756738  -0.4916043322  -0.4810510805   1.0707375158   0.3489285198  -1.3870586909  -0.8473103594   0.3024150786  -1.3594428780  -1.8528480104  -0.5199977528  -0.1717096471  -2.9645403312  -1.6338027922  -1.2849234077   1.9612756375  -0.8101639963  -1.0055361557  -1.8902173037  -2.0562163765 
			 33 	 -1.5050470421  -0.7126272731  -1.3048119506  -0.8162489943  -0.0113751374  -2.0280615772   0.7920077216  -0.3959973654  -1.3819691353  -0.8453103339  -1.1198017135  -1.0900048816  -1.1712862452   0.7370847297  -0.3609335958  -1.0284050519  -1.4431338988  -1.0780177050  -0.3727280058  -1.5252409317  -1.4786140507   0.0708757218  -1.6523682257  -2.0816737188  -2.2549158578   1.1345427861   0.0510901880  -1.3294782382  -0.0970783594   0.7134127742  -0.7716361145  -1.1775040950  -1.7530444071  -1.2148482238  -1.8502713920  -1.4697339747   1.4766400136  -0.5585042554  -2.8037010472  -0.9991290600   0.9098325217  -0.8877309366  -1.6931739406   0.8357126514  -1.0217003229  -1.1105360918  -0.7641125745   0.7764573380  -0.9302944749   1.0334728908 
			 34 	 -1.7414403295   0.4746599915  -0.9582714225  -1.2010743861  -0.7367848107  -0.0115380048  -0.7148266257  -0.0802006916  -1.2848849440   0.4707303202  -1.1457800366  -1.8706428017  -1.7071238354  -0.1111302149  -0.2693966418  -1.8715701381  -3.2110843686  -1.2353907929   1.2117683568   0.0986251002  -1.5886130714  -0.1058353894  -1.1529514743  -0.7953264956  -0.9071630952  -0.1408400004  -0.8064618241  -1.9981490832  -1.1326320268  -1.2253007428  -1.4492731533  -1.2952223467  -0.1425002276  -2.0509890827  -1.3296379932  -0.1191217471  -0.2056090075   0.6182764218  -1.3817068721  -1.8905433024   0.1792929042  -0.5764410051   0.5271883616  -1.4978874044  -2.0047247831  -1.6412751232   0.7842296884  -0.8392078760  -1.2149929946  -2.5460760473 
			 35 	 -1.0937001886  -1.9563207002  -2.0390981845  -0.7527239866  -1.4442385107  -1.6798834229  -1.9525707927  -0.8153809410  -0.9511466807  -0.9452904202  -0.3657718456  -0.2001365405  -1.7527888113  -0.5819494657  -0.6775491630  -1.5294002662  -1.9486561552  -0.6355696788  -0.1286615696  -1.0815440926   0.0813103356  -0.1481706513   0.3731162780  -1.7384696352  -0.0287071760  -0.9248830742  -1.8283239656  -1.2966836589  -1.4738742629   1.2174304897   0.3207158564   0.3662059718  -0.9236618597   0.5983888303  -1.2042123054   0.8358677088  -1.0322676883  -0.3107231398  -1.3709224578  -1.6715966805  -4.1391049820  -2.6759467924   0.0184336266  -1.5096499149  -0.3562862266  -0.3563633112  -1.9453036136  -0.5790322386   0.4657626047   0.1611524177 
			 36 	 -1.1879267042  -1.2376734714  -2.2435415542  -0.6504510778  -1.8415094195  -1.4177736847  -1.7643263137   0.7971584161  -0.9237168013  -1.9137591134  -1.2552916703  -1.3556284793  -0.5030180755  -0.7235190492   0.1639013858   0.3096857077  -0.9745158521   1.3913256268  -0.7954995708  -0.2321157905  -0.8609468078   0.3514940293  -1.4795385099  -0.1313440311  -1.0130355519   0.0286060210  -1.6211751508  -0.2447902938   0.0178756328  -2.7731143009  -0.7859831335  -1.9525458945  -1.4660104385   0.4321775256  -1.3472440045   0.7927618877   0.1635008948  -1.4396338463  -0.1825665524   1.7831159603   0.6468081310  -0.7227372340  -1.1974320360  -1.2796083461  -0.8947317561  -0.8012518678  -1.3953019099  -1.0641329710  -1.3067636403   0.6044878389 
			 37 	 -0.6464276246  -0.9937876818  -2.1488731522  -1.8640186084  -0.1303520168  -0.0788527131  -0.3472742375  -0.3391687003   0.0380895566  -2.0269068479  -0.6553594985  -2.9181653631  -0.9387128032  -0.2922772233  -0.7575977157  -3.2127523339   0.8218948535   0.1717145367   1.1381568927  -0.2231965985  -0.9525668390  -1.7035104243  -1.0976509025  -1.7113474810   0.5662856467  -1.7897021090  -1.9809063101  -0.2973533518  -0.4054759116  -2.6078761100   0.6480488848  -0.4396931500  -0.5479609766  -2.1890800676  -0.7678411720  -2.3045873014  -2.9158986561  -1.6774193893  -2.0380771489  -0.3623486110  -1.7951505144  -0.0912688635  -0.3682867460  -0.7515755474   0.9466202487  -3.5895653850  -0.2257538571  -1.8007538903  -1.8456759655   0.0255802651 
			 38 	 -2.9997348339  -1.2434832643  -0.3944197546  -2.3828303427   0.5429232362  -1.3586301418  -0.5090705679  -0.9141974609  -0.1341233584  -0.1121865519   0.3795690800   0.5662424796  -0.7424276675  -1.0810007037  -0.4273403845   1.3185119681   1.1700582078  -1.2148211758  -1.3913145416   0.2603919808  -2.1212434044  -2.1271208159  -1.7582959654  -0.2127783966   0.6704095484  -1.6931179981  -0.2226825149  -3.3608289642  -0.0454360872   0.0872480802  -1.1826336402  -0.4374543299  -0.7542024365  -0.7502064050  -1.1472547972  -1.5013413349  -1.8402706643  -1.3013134710  -1.1801591690   2.3889016286  -2.3560483603  -0.6341430712  -0.7170488177  -0.4276757182  -1.6610831677  -0.5213533240  -2.0480448939  -0.7641782970  -2.3936148851  -2.5683032854 
			 39 	 -0.8679194305  -0.5836424921  -2.0832109925   0.1272447538  -1.2748368821  -1.9829842111  -0.1434596352  -0.5752368392  -0.8673282875   0.3481726005  -0.6266428827  -1.3869815787  -2.5697304926   0.4520326609   0.3873866824  -2.8163898702  -1.1392759167  -0.5680550584  -1.0831608772  -0.5924919209  -2.6360607266  -0.6115885794  -0.0230698470  -1.3830962939  -0.8006208178  -1.2592694991  -1.1656662224  -0.5831422996  -2.9907592673  -0.7291227978  -0.5921504259  -0.2259833945   1.0375813134   0.7918960447  -1.4568479863  -0.6644154725  -0.7495762953   1.8940074185  -0.9206399464  -0.6990093314  -2.2112915489  -1.7840300110  -1.2497521044  -1.7014124395   0.2154423145   1.0100787664  -0.6430985069  -1.0500424645  -0.0834807335  -1.5719804662 
			 40 	 -0.6029480041  -1.7072156016  -0.2646747027   0.9027220552  -1.1318578463  -3.1015139497  -0.2703161528  -1.9060022889  -0.2217236128  -2.0862397967   0.1220651308  -0.3937027613  -1.1677311036  -1.0455626073  -1.4656264656  -2.0167242976  -0.9181578754  -1.6334120235  -0.6996260395  -2.0540705071  -0.3199675334  -2.8001462899  -0.6560900906  -1.3084000466  -2.1099124793  -1.3718564397  -1.5560093981  -1.3743700816  -2.5424990616  -0.8424679270  -0.5502853409  -0.7174186658  -1.1754171107  -0.4266973307  -0.2460096330   0.0642071400  -1.8303743028   0.4925547237   1.9451641194   0.3260646429   0.6741005576  -1.4912946804   0.8770483150   0.4302069680  -2.3237159121  -0.2420810106  -0.3621873144  -0.0075737998  -0.7147122761   0.0393649078 
			 41 	 -0.9461831276  -1.5604641731   0.3473132196  -1.4035737306  -0.1243494041  -1.0697189297  -0.0136833562  -1.2573859310  -0.8129585461  -1.2945881136  -0.2343499043   0.3235367999  -0.7543514182  -2.5884606534  -1.1507639056  -2.7904176764   0.4183426603  -1.8790018104   0.6344602547  -2.5567406137  -2.7129143947   0.8921370293  -1.7591947562   0.0724385070  -0.6545007010  -2.8659604251  -0.8101510443   0.6576578075  -2.7820845185  -0.1450244570   0.3627365221  -0.9975722028  -1.6828727365  -2.8518690624   1.6310644373  -1.4281668280  -1.2329340396  -1.0355883402  -1.6463568230  -0.5560939027  -0.4842293268  -0.9394334483  -1.4451997870  -3.2744130974  -0.3962825710  -2.6346287853   1.6331564060   0.5124563635  -0.1590745349  -0.5159720991 
			 42 	  0.1713325349  -1.1159075568  -1.1137700490  -0.7088851527   0.5359242378   1.1066783180  -0.1179002571  -0.9907485573   1.0074743586   0.8919502652   1.5088758348  -0.2255964673  -0.6655763451  -1.1476577444  -0.8698285102  -0.9827579890  -2.1512188946  -1.7398847616  -0.4104158898  -2.0839162973  -0.3171145266  -1.9938263299   0.1868610166  -1.5080310808  -0.4231086693  -1.3742033793   0.3167809022  -0.2639143713  -0.9294790103  -0.9916331887  -1.6373346539  -1.4403759614  -2.1721246855  -1.3721613364   0.2232352008  -1.1017037438   0.4318202484  -0.2686466519  -0.3837544259  -2.0672755621   0.6020310402   0.5597172579   1.1330320469  -0.4320337853  -0.4096100732  -1.2376508796  -1.4103268168  -0.3941145777  -0.1697686981  -1.4040132806 
			 43 	 -0.4039106628  -1.0783548154   0.5745974281  -0.3610310214  -2.0544111078  -1.6306470802  -1.4535896847   1.5479635144  -2.0835195547   1.0627903624  -2.2114968652  -0.5767426158  -0.1300421625  -1.3669932525  -0.8849824160  -1.3586674402  -1.1183345151  -1.4213372433  -0.4841009401  -0.2093671026  -0.8216516300   1.9044585452  -0.8309987437  -0.4557617389   0.0235482754  -1.9162524103  -0.3058155030  -1.1649685098  -4.7989915151   1.5550221670   0.7216900292  -2.3269004722  -2.4973046366  -1.1025992732  -1.5542234043  -0.8494102954  -0.1855338543   0.6164692317  -0.1698978307  -0.9474267493  -1.1023598995  -0.2401633088  -1.9047788756  -1.5564498626   0.4553508588  -1.4012451926  -0.7213499183  -1.8168338414  -0.2604528056   0.1832659104 
			 44 	 -0.1426219493  -2.0010334246   0.0392980359  -2.9376071105  -0.8031613789   0.8793478055  -0.5933594588  -1.3811148865  -1.4990117813  -0.8269331108   1.3180924358  -0.8436533827  -3.1239009444   0.7529416776   0.7337354657  -0.5989917670  -1.5759577269  -1.2986304165  -1.7279784795  -2.1457025763   0.2329196191  -2.7373094860  -0.7440219117   0.8250526632   0.0036682371  -1.1052927790  -0.3484971582  -1.9085771602   1.3901251018  -0.8977561614   0.4442774136  -1.3469446662  -0.6462208404  -0.0950299033  -0.7068288935  -0.3623486110  -0.8027192170   0.1196644781  -0.8864747013  -0.0277296290   0.4683389648   0.0049558541  -0.0872307991  -0.7843812948   0.3058571707  -3.0257149037  -1.8737708340   0.4723966608  -0.3440302955  -0.0033100099 
			 45 	 -0.4263590110  -0.5203738130  -0.1124267706  -0.9855730455   1.1833298931  -0.9180835839  -2.1457913833  -2.1408095137  -0.1840030341  -2.7314974900  -1.8240337383  -0.5600577086  -1.0239573461  -0.5385510675   0.0048761354  -1.1182744295   0.4793025934  -0.8260461846  -0.3753575054  -0.7402499962  -3.1533523865  -1.3023148273   0.0637366600  -1.0362454812   0.2659719178  -1.2134678023  -0.8680970666   0.6508424674   0.0743154969  -1.3593885198  -1.9232751883  -1.7849878583  -1.1633372387   0.4576601911  -3.1006728309  -0.1049079239  -0.5854857570   0.1092664572   0.3089292952  -0.0224702165  -2.7063359710  -1.5816251428   0.3558774338  -0.2234206830  -0.5267873516  -0.0339494468  -0.4584782140  -1.1060683695  -1.7702063669  -1.3523963254 
			 46 	 -2.4555746128  -0.1067254623  -1.3396962415  -1.0461870070  -1.6616040399   0.7866768038  -0.8728711715  -0.4046883900   0.1300323801   1.5359570941  -0.8442309967  -1.4410512592   1.4052691954  -2.7162517383  -1.5792277756  -0.8027005011  -0.6237691398   0.0275508802  -2.0807927839   1.2995670752  -0.0453768441  -0.2288648983  -2.0403949398  -0.8027964479  -1.5245995146   0.3259525110  -0.4611161877  -0.1771317864  -1.5621077969  -2.0903416230   0.0759824042  -0.0081722279  -0.9497524944  -0.5148339769  -1.8942601380   1.6373300679   0.5779117431  -0.7347294786  -1.2987048570  -0.1773362303  -0.7458969860  -1.6284204635  -1.3193308077   0.5847113262  -2.0039660963  -0.2690851838   0.3004809212  -1.2816047465   0.8471539402  -1.0367815244 
			 47 	  1.5145165367  -2.2587596461  -1.1109262264  -2.3216817369  -0.3765245251  -0.3841431124  -3.4735236170  -0.2486506060  -0.7197610927   1.2002098652   1.1674403385  -2.2364669571  -1.5754129414   0.1118233886  -0.3614629190  -0.3871236021  -1.2129096100  -0.7144701058  -2.1302194182   0.8153004769  -1.8519293926   0.1422918994  -0.9122692915  -1.4773978213  -1.8859381856  -2.9330504620  -2.6811148701  -0.3260248219  -2.7441223025   0.3116176003  -1.4786345215  -1.6632327231   0.2314170157  -0.1861457527   0.2469657121  -1.7220942414  -1.5011114316  -0.7837709392   0.2753572242   1.2165189399  -0.3535610412  -1.6356312030  -2.0313033574   0.1117846521  -0.9212649729  -1.8006457304  -0.3102830104  -1.3277353380  -0.4813468089  -1.6357289824 
			 48 	 -1.1995718623  -0.3794097548  -1.1637768532  -1.9337479183  -0.1565734294  -2.7622323394  -2.0928901009  -0.8313836462   1.1321579241  -1.2198649740  -0.3780037935   0.5566085538  -1.8587322938  -1.9747571051  -0.6239114238  -0.1441290467  -2.2924499537  -1.5300750039  -2.5679469784  -1.9266547317  -0.8848774216  -1.6540277553   1.1401475536  -0.7029485850  -1.7270553652  -0.7211613251  -1.6538934299  -0.7787278954   1.1836047454   0.0244309138  -0.9698676663  -1.0473915852  -1.4989257842  -0.9528142856  -1.5347336413  -1.3932156909  -1.3265798900  -1.9701204980  -1.6371516857  -1.5370920433  -0.3963208272   0.4925850550  -0.6744706309   0.1610365731  -0.8100819545  -0.6234954660  -1.2999855169   1.0425202468  -1.4112370684  -0.5335718183 
			 49 	 -1.4453736796   0.6534120186   1.1179344397  -0.4529975350  -1.9416667407  -1.1603782320  -0.5228772232  -2.5692097915   0.1465908691   0.5480031098  -0.9304401588  -0.7703464537  -0.3867700911  -0.4420505566  -3.4007652336   0.0704147262  -0.4883330363  -0.5914589598   0.2007647527   0.2022499699  -1.0209187356  -0.8504421226  -1.7333644705  -1.4504277747  -1.7058038753  -1.2965576944  -1.5319475126  -1.2125863610  -1.0361325996   0.1490835365  -2.1342138209   0.3532244769  -0.8017761118  -0.2554856939  -0.4191576832  -1.6826456072  -2.4500963148   1.5531799962  -0.6334826809  -2.0341475016  -0.5634914346  -0.7553023926  -0.1874126344  -0.3792275794  -0.1812385810  -1.0517894223  -2.5488466262  -2.2179498257  -0.4496278419  -0.4138290395 
			 50 	 -0.6749768598   1.4190544303  -0.1691050715  -0.2017168488  -1.5116694161  -0.1449884954  -0.9617789640  -0.6241828591  -1.1064015124   0.0905024001   0.6921195253  -1.8992588578  -0.5440806268   0.7638858163  -0.2613779533  -1.0781232265   0.1955596632  -2.3964704434   0.3081212041  -1.1177959090  -2.1435813049  -0.6759798609  -1.0697111340   0.5853198053  -2.2498204287  -1.6177463459  -1.8263651334   0.1244553644  -0.0545939953   0.5856803404  -0.7454907725   0.7204608185  -1.7945253163   1.2330363007   0.7992678869  -0.2570599930  -0.7018522792  -4.3159742996  -1.0959736827  -1.1932630483  -0.3621307316  -1.6931122403  -0.8052563503  -0.9607304514  -0.5155413888  -1.2315896624  -0.4617482400  -0.4015344491   1.2133100504  -1.8082536661 ];
	Beta_AT=Beta_AT[:,2:end];
	 Beta_FEE=[1 	-27.9059084152 -38.1562227016 -18.9718150980 -33.6177558837 -39.1468832389 -14.8575984580 -14.7217042415 -42.7627555555 -23.7737697396 -12.4880575838 -36.9782353101 -67.2300161495 -28.4431964415 -37.4171240120 -50.7201024092 -47.4622116246 -29.0216207012 -42.3567207000 -47.1374511873 -47.8849005458 -50.7588758281 -35.8503406322 -30.4665016495 -12.9929273637 -22.4822870887 -41.5872121106 -51.7896764812 -48.4245999055 -23.2583098260 -64.1910697080 -39.8782381440 -21.1751645298 -19.1258277193 -38.7110213806 -41.0896309719 -24.8733469930 -55.9776482998 -42.7942327446 -50.3715370856 -49.2171952660  -9.7235920795 -51.9894313575 -39.5889960180 -35.2998021173 -32.0578005898 -34.8099416509 -34.5228220712 -35.9272737805   2.8807321556 -23.9692017781 
			 2 	-23.5126897028 -15.2961923304 -16.8623116624 -57.6110853108 -31.0466499235   0.4977695870 -28.1980123938 -51.4645389659 -39.3926872794 -57.4382076804 -31.9334860107 -42.9090541716 -26.7510983475 -54.3143594661 -42.2541851970 -30.1498961506 -34.8361453845 -25.3091740574 -46.4816777629 -55.0347749658 -30.5297256430 -29.5610840247 -49.0719347039 -30.1813142925 -27.2949658248 -32.3965480804 -24.9841867232 -25.7711695659 -18.1649186625 -42.4021201350 -19.9961516898 -28.1049172599 -28.8703516583 -55.4929519008 -15.5995992960 -41.4808435332 -37.9421936623 -28.0830669264 -31.2264279615 -34.2913781597 -55.0246344691 -33.3350236141 -39.3096809598 -29.2098472344 -59.0131373966 -53.5655895488 -37.2684585095 -40.7928781110 -28.0823063799 -37.2940753889 
			 3 	-28.3420862160 -33.4641023624 -28.3045957318 -42.5363065373 -43.1261382002 -21.8891536398 -16.5021508560 -65.1286449106 -33.1810756189  -4.3594187659 -44.1025569449 -70.6123557194 -37.2734273577 -30.6669892616 -16.2373750593 -50.4571107080 -20.0031502991 -68.2331802682 -14.2862836522 -19.6964222670 -38.2534861245 -23.0149632941 -35.8415689841 -37.5410769621 -20.1094179756 -55.5240581755 -13.3944543108 -20.8792784770 -27.4399565948 -23.2731945259 -33.6597008754 -38.1478614242 -45.8709075242 -42.3384293499  -8.2735378301 -13.9344045757 -32.9386555165 -49.8108137432 -28.3407700882 -19.8163969436  -2.2391907195 -53.7211900905 -21.3365269308 -12.9143517084 -38.7549253161 -70.7528453798 -60.4388203017 -25.4097524076 -27.0372907702 -32.4219789548 
			 4 	-43.3906454196 -26.2293948367  -2.5358540149 -52.2697158808 -46.6623024257 -23.2147683435 -11.4321470663 -50.0362517050 -17.7977243148 -22.7723333225 -18.0050223097 -25.0119954581 -27.6754530605 -22.5340259168 -12.6568806269 -35.5503001267 -47.0033559712 -36.6011467794 -37.9492123857 -14.5276975730 -37.7125675065 -26.4992633946  -3.5619081196  -6.4228485100 -40.3829737350 -12.3180510497 -30.7795484697 -46.1941499105 -49.8870802941 -13.2684062527 -31.1964580531 -35.3802333741 -45.0438375113 -26.8592833557 -31.7480413141 -41.3321921437 -48.4355917165 -10.8653476997 -26.2020806314 -37.3975763792 -33.6996752379 -62.0117109943 -22.4762643944 -28.0174744681 -11.1472035181 -39.1739061719 -21.0618515188 -37.1417510754 -35.7374071688 -22.9289168231 
			 5 	-45.5953736509  -3.1376498126 -38.5113843731  -7.4957472256 -24.0959630982 -26.8773472200 -36.8411815110 -31.4947511570 -23.2917727024 -40.2506571194 -35.7021424579 -20.0953243732 -32.1996449613 -23.7651456654 -46.4854336251 -51.8499870937 -23.4747652297 -21.1885947599 -36.7411611750 -33.0105279769 -17.7621137707 -46.1415217749 -10.6564337363 -22.2893173836 -37.2696714956 -36.7001454893 -46.6664898255 -20.5659854220 -35.2407025060 -11.2505157057 -24.2459660343 -28.8531198484 -45.1806466244 -32.5010806272 -13.5500889468 -22.8376507275 -31.3561523902 -29.7040212565 -24.2440743204  -8.2791204336 -37.9155456868 -20.3897502338 -21.9959809959 -15.7102838231 -16.2452888656 -20.7633004493 -32.7646947377 -28.3646176688 -13.6071957441 -41.9657441900 
			 6 	-33.6879575920 -32.5484946023 -34.2011840449 -36.3063898724 -20.9768899955 -53.8927016677 -49.6160697799 -11.1062409909 -24.5116184647 -21.8449800808 -54.3254711143 -43.3816482949 -58.4941859388  -8.9626200784 -16.1138319805 -20.7102109646 -15.7469622221 -29.7420136053 -14.9192987366 -45.9678145710 -32.5384878056 -31.0646670284 -34.7623704603 -29.6435401817 -51.9427690916 -23.0035069931 -39.1259767187 -30.5707249072 -37.7006230145 -29.0181573101  -4.6469707617 -40.4973825983 -53.5751389876 -33.9090488786 -62.9326527282 -28.2450564106 -60.2805388318 -25.8080178239 -58.5183443265 -37.3553899008 -42.5559444711 -52.1412570373 -66.1048808109 -34.5931341090 -52.8391193403 -52.2856518361 -36.2382214100 -17.7500162004 -45.6258803735 -26.5785942470 
			 7 	-21.0808027162 -39.1930482193 -66.0229568879 -56.5190352736 -30.7392087410 -55.6337661872 -41.0287798279 -25.1098380087 -61.0395997913 -44.3480474819 -27.1035401300 -30.2217945041 -24.3043667325 -50.4586508679 -22.2457240858 -10.7636280844 -49.5211416929 -40.8208752391 -21.3752294236 -27.7624396720 -27.2050315271 -42.8000302306 -17.9210488340 -39.1360956826 -38.7671933853 -38.1626011668 -16.7026823236 -28.4604978499 -42.2683995411 -29.5015329837 -22.0106548127 -20.3654795995 -30.5581672727 -15.2287577311 -31.9288496953 -31.8845946875 -18.3444610984 -21.4119305396 -33.9738803866 -44.1756465837  -4.9335714733 -20.3115267103 -48.2010769867 -26.7507230631 -36.1246396047 -22.5447952400 -46.8467211211 -27.7652796947 -45.4017197007 -25.9997502506 
			 8 	-39.5858142090 -51.7567565064 -11.9931037869 -37.2508009571 -47.5818675299 -20.9844861000 -45.3941228065 -20.6973253169 -37.5266658986 -16.9517346948 -34.0101268879 -47.3112739960 -48.5201758932 -45.5572741842 -11.0827877859 -39.8222028540 -15.9236935972 -26.4542937320 -38.6801564142 -21.6952558320 -36.7316198412 -20.1799257712 -58.6685686013 -36.0091948518  -5.8575089813 -40.1133319236 -27.7583206417 -48.0208751366 -36.4768615511 -31.0524052698 -16.5709157023  -5.0057334197 -33.9281655786 -55.5950652152 -22.8796391278 -30.8568674521 -53.2462800618 -20.0334563406 -12.2842458295 -35.1613507866   1.6799993282 -35.8832169886 -34.6941101447 -40.8850920481 -17.6611718203 -54.4220431296 -49.9999273255 -62.1841352324 -25.0565130644 -31.8026281509 
			 9 	-44.9759084629 -45.2100402131 -41.9729810510 -33.2278441503 -35.2082130610 -35.2545414284 -17.1401689715 -29.1991711199 -33.4087970323 -31.7561514439 -37.3113324942 -18.1198603921 -20.4130677894 -47.9025295244 -25.8581849135 -30.5204907378 -22.8954179121 -17.3166824109 -23.3496459129 -34.3222378272 -32.4256682779 -28.7427188485 -26.6724234527 -25.9080620742 -42.8609708787 -38.9037238413 -35.9541983793  -5.0796260578 -47.2466881907 -62.4706122326 -27.0592320781 -39.7110738290 -50.7037879324 -26.4459419773 -71.3733436450 -15.0018790369 -23.6606846575 -31.3276329622 -34.8250001643 -15.4062796300 -27.5123438912 -36.0113276280 -31.8806817592 -36.1840114252 -50.5736520263 -45.4640793544 -45.2309990258 -22.6314310421 -24.3609547722 -16.2053995429 
			 10 	 -8.7938133693 -27.0808728775 -29.9279799224 -41.0895133409 -36.4100075754 -30.3881929778 -28.5724878856 -43.3464065601 -45.5493500768 -14.9196003252 -37.5513908044 -37.3710021864 -25.8068017296 -26.1977739319 -43.7987852547 -30.9387018403  -7.2814795132 -34.1786303689 -27.9635413141  -5.8106024327 -36.6691165410 -27.3326393951  -5.7183498283 -36.1736425644 -33.4255493497 -49.3780293274 -28.2961208860   1.3218273034 -31.5075554751   9.7255419604 -29.5314265054 -45.1626185961  -1.4137883012 -49.4700193716 -28.9086240544 -15.2358330205 -37.6007412445  -8.6740682326 -27.6886633738 -41.4206815040 -23.0545015083 -26.3344322119 -32.8381225650 -37.3353108816 -33.3476442847 -22.1646611750 -49.9217697523 -18.5078922357 -14.5711456281 -54.1421502554 
			 11 	-44.4123609439 -32.7209537972 -38.2417705275 -37.6849594238 -19.7214571219 -27.9126294593 -27.1228609262 -66.8380431535 -29.4945677395 -27.0709673538 -51.8715507161 -20.4903121851 -19.3706280941 -17.7027319356 -24.6052763306 -21.0329027835 -42.9777264692 -52.8463138697 -29.4917690573 -21.1645802512 -11.6727563378 -40.2440013862 -24.3556000535 -31.2181295554 -20.2571035346 -54.2470354410 -29.0554592985 -14.8276006002 -17.7685057829 -33.4874826085 -27.2057438397 -36.5737713660 -25.5172199826 -62.8024693824 -29.2535311029 -63.5205325174 -34.0214932173 -17.0630015310 -32.7786430627 -25.5565282426 -45.5095896250 -15.8103280371 -34.3380524909 -38.9473428091 -32.1979420049 -21.7818645088   1.1532401446 -19.1668033727 -39.8523814483 -40.3808327102 
			 12 	-47.5304668791 -48.3510315692 -37.8073622574 -35.2270873491 -32.0596826868 -27.0808728775 -52.2541710207 -48.6124288092 -34.2771713068 -17.1426566415 -36.8044949807 -40.6836821355 -20.0337199795 -29.9089610852 -26.9866868725 -31.7364154165 -19.9910686242 -38.7440826539 -24.9488945980 -52.0353723920 -65.4548034084 -49.3397000142  -8.4191868799 -42.5251088187 -15.2980867427 -24.4657623091 -25.1327457555 -16.4810244233  -5.5904908171 -28.8834579476 -39.8342742689 -32.1790087079 -39.5754305641 -22.9177749163 -44.0319589767 -25.7516857523 -24.0765281103 -58.4252784844 -15.5993393592 -35.5355493648 -40.9718328161 -27.9801373264 -55.3821071083 -21.7690718140 -23.6306475620 -31.9877890993 -19.3577593923 -17.2265440267 -11.5447949193 -22.4175753508 
			 13 	-31.3626640439 -24.7589766046 -36.3203243767 -17.1221179601 -38.4909136852 -40.3547879958 -37.8502366021 -46.2852276304 -54.8536798578 -32.5433667595 -12.1701967657 -31.3061635311 -40.0622386446 -54.1038910005 -28.9794920580 -21.8347722877 -12.7582495568 -25.1400669300 -40.9348454697 -41.9600385250 -53.0750947605 -32.9756826229 -59.6714228358 -17.7166398921 -18.5997222257 -41.0113676515 -47.2982143710 -45.6674197162 -23.2189148092 -26.1525860957  -7.8061228135 -42.3488687101 -23.2446869540 -33.0150736844 -24.3881706188 -43.0181171790 -25.4738942467 -40.2433571463 -22.9533593903 -35.8991293117 -62.1413263673 -30.1740030739 -37.5852750285  -7.3866496596 -36.5494442821 -41.3505020765 -13.8101717695 -11.1903567997 -48.7752163857 -14.2307234795 
			 14 	-28.3790657369 -36.6639767724 -48.1592465620 -61.1888008194 -38.9295815466 -20.0853720581 -22.7885741697 -53.5540683503 -30.8187582819 -52.4824514850 -13.5171242422 -43.8795796598 -11.4948305902 -35.2089880368 -14.2565054811 -18.5223156905 -36.0411852226 -37.1011541371 -15.9608345282 -32.3438410319 -20.2725119800 -38.9950749857 -21.8779772921 -30.8288512008 -38.5880843507 -41.2701108650 -18.9942513013 -37.0062655243 -43.2012281365 -37.2765996564 -47.9327303114 -32.6066404594 -44.9106557004 -44.9874648069 -30.5873002962 -27.8230424686 -51.6395941101 -38.6934578281 -56.3811563097 -42.6064976960 -43.7950381174   5.5881985372 -47.2127072449 -29.8205994192 -19.1430349342 -20.8107861845 -42.4129627648 -51.3533692584 -47.8168291200 -34.8148555149 
			 15 	-33.2727589109   3.0363632290 -24.9668240197 -39.0665769979 -19.0466018256 -31.7666058547 -41.2344255467 -41.9311320785 -25.9789235204 -12.6783392514 -28.5043645773 -17.2857816070 -34.9918303140 -26.1220981195 -24.8477698636 -21.3460629926 -29.2465649892 -19.1699172565 -55.1531128767 -31.9955895241 -26.4596929015 -13.2345700158 -51.8306985345 -50.7771803594 -25.3692159657 -31.7582322097 -26.7770729314 -53.6573406041 -10.0012797300 -46.2425666220 -40.5184240250 -25.1241567892 -34.0738514422  -9.6379116699  -6.5164741552 -53.6237661541 -38.5835819664 -24.4810269848 -37.8491089027 -29.8450040393 -43.2292631519 -15.1273270356 -30.9751459789 -22.8843612629 -46.1000453742 -48.7291722845 -53.3317068123 -39.9610107560 -17.0072591758 -41.5058428291 
			 16 	-23.1211273604 -49.5856770446 -41.5957944045 -47.3114243259  -6.9732663854 -19.8359670293 -32.9725692715 -39.2745407723 -29.8447835401 -47.1911234878 -45.9022297950 -50.6282392404 -12.3311130907 -30.3161163780 -41.7881472119 -46.3163454340 -14.2381732935 -34.8751272297 -13.1397743839 -25.4960860711 -38.5429402592 -59.9114066778 -40.8830262845 -20.8153532743 -18.7924954981 -25.0745164882  -4.9066929423 -36.3334239521 -30.6908433102 -13.5128896452 -38.5167302851 -31.3138447963 -36.3286953171 -38.0134436728 -15.1484451510 -27.6269416219 -23.2624164826 -47.5612017079 -23.3259110973 -57.3228870855 -33.4580698168 -23.2475236058 -20.6377584789 -24.5928977170 -39.2200841033 -54.9106303076 -44.8192574646 -14.7688721973 -19.9354981561 -17.2350292476 
			 17 	-21.0672851395 -12.3841944735 -19.4540611250 -42.5023786676 -21.3626730207 -33.2772209923 -66.3577813542 -21.2121132858  -4.5178978888 -17.4431006566 -39.4796834523 -29.4920029286 -33.6031452674 -40.7909364434 -25.2722410955  -5.2988258968 -42.9699675250 -30.5943152160 -16.3295933586 -10.3795953288 -40.8912365470 -33.4905810355 -45.1148145359 -13.3007956040 -55.3988314853 -13.6838706324 -33.8315931683 -17.2159971657 -16.1860759736 -48.2694316105 -38.8000580134 -38.6695017596 -41.0997361702 -58.1174691200 -47.4539575385 -16.6284478585 -32.4602584226 -28.6664755166 -21.8796890029 -55.0135567381 -14.6061984364 -23.3847478366 -35.9170956772 -24.3618047168 -24.2366967684 -48.9787299182 -31.8447753193 -42.7629629277 -22.2223831100 -34.2118331098 
			 18 	-31.2834526166 -29.4790093579 -30.0977452730 -30.8259017327 -33.5206645929 -37.2408242998 -11.9246726340 -20.4679941888 -35.7976180823 -36.5333091609 -43.7703583637 -50.8245583475 -44.3444066696 -26.6862544627 -32.1807581524 -38.4956274748  -5.5179673510 -42.7291232669 -17.8530497615  -5.8694001740 -18.5340639543 -27.4245645343 -26.0234134396 -48.5011627469 -28.6315921341 -13.5839623585 -17.7987619941 -34.7185924328 -29.1434265146 -37.4119116103  -8.0659474693 -35.4971126985 -15.2681049171 -33.5173887291 -40.5620077040 -38.0721310914 -30.4660849111 -16.8590957385 -21.1947273711 -37.3302841415  -2.8071767234 -61.6131071061 -28.1172537619 -27.3141820704 -44.8091169277 -27.2527668478 -29.1587232289 -47.0620661909 -40.0333960783 -31.6006219655 
			 19 	 -5.9448897098 -42.2349321884 -34.9019517983 -43.2080943333  -6.8238982887 -28.3380114601 -11.0400596704 -11.6309132903 -18.0737627119  -8.0796004157 -49.7560669528 -29.4705074188 -50.7811856420 -24.3982106844 -21.5385137377 -36.3504484154 -42.7554194573 -25.9726553360 -27.8463738678 -34.7613362138 -45.8373161780   0.1359931004 -18.3165686231   0.4401240117 -31.2268525989 -43.2601965444 -44.4579556462 -33.9970300536 -33.3012290627 -26.3929146622 -27.4989759015 -10.8178498906 -31.0523071526 -28.1971285263 -61.5817556573 -43.3914927457 -40.7997339033  -6.9997766809 -32.2554770382 -38.2046365313 -49.3011165146 -29.8057812826 -23.3344714012 -45.4627641489 -53.7247535807  -0.0243005095 -30.1071980361 -36.0405186090  -7.5812151449 -48.9712971642 
			 20 	-38.4908519101 -38.1791793317 -23.8666765592 -21.0129949562 -46.7015151974   0.6475555129 -17.0107079386 -30.4040649302 -37.4382373585 -24.9178241316 -52.5082041988  15.9661319063 -27.5204881965 -35.1962531101 -34.2006297923 -17.9043719018 -42.6657406841 -24.1823384533 -45.2869701944 -32.4033979459 -46.1611702366 -30.4891095952 -33.8975461034 -24.6239461001 -45.3038667991 -30.9754043214 -20.3033857605 -27.8050186362 -30.9189164283 -43.9519329026 -62.3735830269 -41.9918673515 -15.3713210247 -30.8544915829 -13.4739305150 -52.1417640302 -43.2089086989 -15.5598733381 -12.1317575809 -29.0132255455 -58.2428549434 -37.7531026623 -22.1464910694 -54.0758013544 -38.7216308326 -43.6279778172 -26.0581429721 -26.3523729669 -32.7182823889 -50.7545025758 
			 21 	-38.3079486830 -31.9392181515 -31.4641276604 -29.3310291130 -41.4864488819 -11.3119530246 -38.6149468794 -36.7214417273 -25.9307846410 -18.2227738993 -44.2382206617 -18.0979505927 -45.9164111630 -16.8929431869 -26.3798140258 -29.6274239887 -38.3020643704 -43.5669525445 -40.5046392610 -34.8895613169 -41.6050547400 -27.4701621388 -27.4104391115 -20.3538187175 -40.8407566380 -18.0266861165 -36.3274942836 -10.6832069231 -37.9843072821 -37.2361246258 -17.0941260503 -22.3646363766 -28.7320491851 -21.5087063483 -36.4335872985 -30.2461809518 -25.2163258505 -15.7596381760 -38.2043214400 -45.8351913851 -23.5131904564 -44.5767406877 -44.6807848539 -55.2985259509 -44.1815005641 -30.6484501864 -36.2200958783 -30.9648654202 -11.3401050373 -40.3511951420 
			 22 	-19.4324926674 -30.4555967926 -25.2499251100 -40.2395879815 -41.9943082091 -48.7183666703 -24.3754586271 -27.1666423087 -40.8198798739 -38.6421460015 -33.6374364580 -61.3061371728 -39.5589005639 -24.2767157083 -32.4714130316 -45.9958763985 -43.0221622379 -67.4810796170 -61.4303504554 -40.5655644299 -50.2889284710  -7.0840597477  -4.6669234657 -18.2231162568 -36.0749909714 -15.8330827477 -30.4168848436 -36.8892334821 -19.2388490871 -37.7470503736 -33.7288717540 -32.1356917204 -36.7687242798 -39.1335442991 -26.0409896408 -34.3279136609 -40.7112142781 -15.6734646986 -45.8257911520 -32.7155136056 -17.6444240580 -34.9824783466 -30.7427666599  -7.8880083898  -8.7508557107 -24.7747515670 -34.9375521000 -58.1854601912 -30.9801589044 -37.5826815696 
			 23 	-35.5427787283 -46.6879020779 -17.1158877020 -44.2575857495 -62.2164593360 -45.6319571764 -23.8708396489 -28.9799444111 -40.6163339223 -60.9604991292 -20.6354618338 -46.6002327082  -4.8000286680 -36.4721884163 -23.5393045517 -12.6054320250  -8.7396601649 -18.1666486063 -26.8690073448 -23.0732645493 -41.5918906424 -35.4979629128 -41.2978121072 -16.6352860585 -19.0371151981 -21.7509972576   2.1415190634 -14.8017507322 -26.5292553926 -45.6736370296 -17.8879184816 -58.1189036522 -23.0461112537 -19.8291729823 -49.5723882357 -20.6527261754 -37.6493948496 -37.6916104216 -30.8419802362 -49.8386564663 -28.3763033281 -44.6216194769 -30.4831709495 -42.1615314083 -42.4263301624   2.5159930450 -43.5742552680  -9.4997447940 -17.7532938293 -27.0461648662 
			 24 	 -4.8215198767 -39.5491183391 -25.6414517516 -26.0086309316 -36.5243189781 -28.9377979186 -46.6506373772 -29.9618742769 -31.8209322311 -52.5330324343 -45.9914659458 -35.9614641297 -26.3164897731 -40.0580030868 -25.6217945654 -57.1359721798 -27.9135074002 -20.2353575717 -27.9500186505 -46.7777293889 -34.7301604473 -34.5867124498 -46.8554268287 -34.5214714298 -50.7966403462 -20.8213073927 -20.5891840281 -21.6325718531 -20.6300674544 -32.0630481863   2.1716925703 -46.3865847623  -8.0916980696 -37.6136551046 -40.7279252758 -49.5815800113  -5.2913842491 -23.9358461064  -0.0931298901 -30.6347088056 -39.7976174958 -25.7046692071 -13.4910057725 -15.7831866027 -36.9651470489 -34.2996823860 -35.1885594728 -25.5323813873 -31.1499944453  -9.9641261353 
			 25 	-51.0325854543 -12.0647147243 -31.7179990278 -10.0140686800 -19.9375785172 -39.2058650056 -38.0815777746 -39.7277696074 -44.4295665533 -20.5275983520 -19.8103923681 -27.6361929492 -27.8799641763 -26.4502160898 -35.6897829434 -20.2799532316 -33.9949469592 -25.9801764077 -30.2713408719 -16.2490366618 -21.1139417859 -47.4218188111 -29.3499184759 -23.2260321919 -41.1066703153 -27.3561993569 -29.2853916395 -21.6957944891 -41.2464438506 -34.5275709375 -45.5748841629 -39.4368352431 -16.0012386669 -25.6147694127 -38.9037075100 -38.0197709122 -34.7401102975 -37.8872815780 -55.2391186685 -47.5096440108 -37.3921495300 -32.1841125566 -21.1941890915 -32.4382048903 -19.6758675738 -24.8118038376 -54.5273626185  -3.8759377689 -24.7674765203 -37.9154468603 
			 26 	-19.9133362691 -49.8197692943  -7.0243508788 -17.7451915220 -37.5221063889 -38.8880273273 -38.3282249758 -46.2097032620 -64.6439604257 -23.3552336142 -35.2843559566 -27.8427987893 -18.1029608268 -28.6398373535 -45.3962336430 -25.7937744610 -35.8009116200 -18.9398983821 -10.5314333043 -21.1652743085 -40.1421976673 -15.5102875290 -54.4667439281 -66.6923414138 -18.7793214988 -15.3463128913 -19.1548522137 -36.0867885411 -25.5297503718 -44.9816023913 -43.6696891948 -17.0829057920 -24.0010580229 -40.5146669150 -19.5958023164 -38.7701557459 -30.9809685031 -48.2732596993 -27.5324692195 -51.4758827182  -0.6136511182 -28.0280111969 -12.1773437599 -11.3326214097 -10.2222198892 -51.7268267074 -10.6883110321 -36.1904238562 -26.7096633726 -18.9268693990 
			 27 	-34.3306656468 -42.9197674497 -54.1992112048 -48.8669821583 -33.4476264469 -50.7289586655 -36.6075247322 -28.1851057159 -51.9208616902 -13.8242930168 -41.0557137070 -21.2816307986 -22.6383077930 -35.7763400742 -34.7098531327 -25.6880807635 -58.8996104070 -39.5546434019 -23.1173000786  -8.4427821818 -25.2642597700 -35.4674242642 -21.8349050597 -22.6158880748 -17.2419406216 -39.4009446344 -13.7688954295  -1.8829068910 -55.3699278665 -24.1623877957  -5.9833591016 -42.3368225661 -33.8497521822 -31.8825760365 -11.3550652775  -6.3580397949 -11.1066775742 -12.8757640011 -34.6661463134 -27.8620894127 -37.0860028311 -56.8341056999 -29.6982064712 -59.4408160807 -21.0614079437 -59.5285694205 -43.0913210258  -8.0122625358 -37.1890709220 -10.6664721183 
			 28 	 -5.9965245825  -9.9993925479  -9.7104709686 -38.9866875812 -38.1597384657 -11.3576600887 -34.4230932669 -42.6291285605 -35.9542207924 -22.7895569503 -38.2303410502 -26.2542502344 -24.7165494673 -44.1008382909 -34.2052550048 -34.2726600843 -46.6180381973   0.5892596273 -17.0455559484 -43.6032517217 -70.1695007103 -52.2115470260 -19.2871955632 -22.4731774923 -26.8934131453  -5.3626177081 -42.6296975173 -41.0769150043 -40.7023681520 -35.5916032410 -31.3416790206  -5.0057334197 -29.6447564986  -9.1366470259 -28.5996160172 -38.6043683236 -28.1122461837 -27.3277282014 -10.2029594199 -29.7228547364 -31.5782031452 -32.6950627839 -27.5726716696 -29.3124908907   7.1481301823 -36.3622620943 -49.0252314509 -48.8601175515 -43.3974463949 -36.8246510785 
			 29 	-35.6881103204 -23.8721288978 -35.1433036229 -43.0808499266 -20.9423023540 -24.0051220971 -20.1251686423 -26.9278959756 -22.8879671506 -52.7325983051 -62.4474373543 -35.9384640601 -51.9833061401 -30.5318836221 -52.0043338036 -34.1498488519 -39.4418485236 -53.0752700980 -40.6318858710 -26.1407118101 -44.8075210625 -65.7848566528 -32.3455077792 -48.0304085261 -20.3292216827 -45.5123286785 -43.1528861680 -10.8568362622 -39.5094297325 -39.6083585435 -32.3194414416 -44.8804461593 -43.4411016751  13.5645384846 -43.9183681817 -55.8265320360 -16.8425115809 -23.5839303610 -10.1837793827 -68.1815530188 -52.8578361039  -7.1605447502 -25.8471258978 -47.3084178580 -30.3279195476 -59.4565733190 -40.0456162759 -24.7975771366  -4.0604965132 -29.9181995596 
			 30 	-26.2365035896 -41.3599944868 -36.0490245163 -21.4677037993 -36.1440315533 -41.9719005651 -24.7689526702 -22.1514002125 -34.1425868304 -38.0371965642 -48.6872286832 -32.1419539168 -30.9487439794 -52.0670224183  -7.4368936349 -19.9553345737 -61.5343626143 -31.8269369886 -33.4627557821 -24.5665235085 -33.8672304287 -47.2158751876 -39.6266045692 -40.0287234693 -41.0315214703 -37.4672666861 -27.2941710178 -37.1987183624 -25.0774751606 -13.5440510306 -56.3676594200 -43.2242789495 -33.3538730287 -31.1997951318 -20.4016402763 -63.2685990874 -34.4707021247 -24.0624635799 -39.1152952908 -39.4468481730 -47.2661695708 -38.4435898925 -23.2942204567 -39.0741905683 -38.8548039612 -31.8757671236 -21.7564752516 -36.8665434737 -50.2003388959 -50.8937977503 
			 31 	-31.5038240810 -47.8224138988 -41.2548781040 -39.9399494357 -32.2570721923 -29.2602041708 -37.2228180101 -19.8317006186 -20.0899751948 -15.8239040679 -62.8805490530 -38.8883478679 -42.4110338209 -29.8230875205 -67.7712265581 -37.6623045404 -18.8017045538 -16.6617372117 -31.5144523258 -24.5505111704 -25.2773837981 -34.6859474553 -33.4565840749 -39.2206792818 -45.7003016731 -13.1836337002 -19.2901475419 -52.9823350547 -27.6075313111  -6.6449725576 -37.7683364424 -30.8345248515 -19.8792589271 -28.0192617391 -18.2774549223 -30.8081604958 -42.7293208886 -16.9103805000 -33.0198740338  -7.5215201101 -29.8252662770 -54.7404975782 -43.8722792081 -56.4601008728 -25.9355419060 -35.1268734383 -30.9347951123 -25.3745949580 -22.9687514139 -46.0322869730 
			 32 	-25.3032212589 -35.6241155617 -40.8393811466 -40.1102414945  -6.2133472257 -47.0098562494 -51.2837895769 -55.5799925131 -21.2754674351 -33.8586681271 -13.4212441280 -23.9203041415 -20.6831141645 -35.6058501234 -27.1871687354 -29.1405553345 -22.8433387562 -49.5976394967  -3.2471609999 -35.1410223450 -25.7903290594 -22.8465025034 -48.2763253011 -28.9217741449 -22.1247884178 -33.6890617085 -31.1005337016 -25.7128647967 -22.4206222535 -49.2151386622 -37.9839730604 -64.0078118772 -14.2912145657 -11.9494686253 -19.2083320609 -30.7273426902 -22.2106805673 -39.6976450236 -24.5248384764 -23.1746255166 -16.6402829905 -46.8127227616 -26.5738278423 -42.4108748581 -25.5565097524 -36.4790328660 -33.8128306169 -27.1568136711 -15.2621875727 -33.2639044641 
			 33 	-26.5174512697 -39.5098502627 -28.5744819250 -17.7152705333 -35.6166123353 -62.6356412093 -32.1212851165 -21.4405235770 -40.1757564043 -37.1411201247 -27.2165460931 -32.5186416961 -28.6171588608 -18.5911525612 -23.2117555698 -29.0762547684 -27.2487340221 -44.5651219824 -36.7593511419 -43.0445287634 -58.9523597460 -47.9064438205 -45.9307320122 -23.3395699167 -45.8843055795 -41.0879624840 -26.7159080749 -35.4956626010 -66.1097038123 -43.0013781230 -16.7899645381 -26.5434784741 -31.8945317133 -13.4814492328 -32.8234072307 -61.1426605678 -41.0836813218 -21.0183350949 -20.1282043094 -23.9079200897 -33.9734182366 -28.0371899244 -31.9883912102 -48.1236330546 -21.5138380195 -47.8327016950 -63.4726725980 -37.2464481813 -44.2875964635 -41.2484813526 
			 34 	-30.0492543700 -45.9533419359 -35.5480594707 -43.2994057158  -8.6018309812 -47.3210127794 -61.0958336357 -32.4969724100 -21.7153205091 -39.1726802372 -36.5613472769 -40.0062080172 -54.7571281214 -37.1096548417   7.4802244787 -37.4838792896 -34.5553062030 -49.7056633047 -37.6327693125 -41.2588246213 -49.1912459830 -46.2177731385 -26.4666005026 -57.7162130510 -23.8185772413 -21.0608838655 -20.1183157017 -40.1370011028 -62.0347995127 -42.4927216425 -28.9207754751 -48.2863952505 -34.0882846943 -13.2808759375 -18.5569506334 -52.4814824163 -37.1474242976 -30.4980544141 -28.0338320951 -33.5202943277 -25.9559830726 -42.0516268426 -33.2498263485  -9.5679820178 -19.5412139794 -46.0738696260 -28.2567398194 -22.8549845761 -42.1968203942 -18.5773191435 
			 35 	-30.4141705538 -24.3752523173 -32.8643163473 -29.5919472936 -37.4588295344 -34.9943948734  -6.0406715142 -43.1584656852 -34.1920383482 -18.8241439156 -31.6603004502 -29.0196576193 -11.9543882123 -32.9660217033 -39.9416013468 -64.7495745257 -28.2334748258 -45.7029937492  -2.8132427287 -17.8509771185 -40.5134442488 -37.3317868590 -44.1917217709  -9.4972988756 -26.3819621053 -20.4278731856 -41.2172713482   0.7315665890 -12.6759355153 -44.9526714657 -26.9944284942 -25.4322188284 -42.3631356806 -61.8322922073 -18.3989489779 -58.7027604131 -17.7962543009  -8.8259891192 -15.5897003047 -20.8966492348 -35.2379236317  -1.6015539873 -37.5067008553 -24.1590690382 -37.1483145830 -48.2881623886 -24.6005141472 -29.4654020940 -34.3821978144 -13.1920470427 
			 36 	-32.9176797673 -41.8771442663 -30.0170331609 -34.2904415696 -32.0955152720 -22.0396235632 -28.2424120618 -23.1464374797 -28.8569538051 -40.6687944357 -43.8878001976 -33.3988948794 -12.6006171491 -64.0076835941 -39.0593508696 -52.9634047716 -16.8805541510 -27.7195335582 -35.2239532362 -41.8453248748 -11.5192825579  -9.0637751058 -33.2128419798 -25.2469173566 -20.5719336992 -36.7400595606 -23.4222178595 -33.0178676428 -44.3949260642 -37.6997180286 -37.5324810723 -23.9796238492 -19.6522973152 -52.6589931624 -34.0843423787 -26.1570674718  -0.1922484953 -27.8968501258 -38.1440302874 -49.9631332994  -4.8057209086 -48.7514359137 -23.2541500493 -48.2946925432 -42.1596056151 -47.0986677961 -54.2927294856 -50.9692579859 -27.6706493806 -36.5328371604 
			 37 	-33.4828626218 -41.0519123871 -25.1521133875 -33.2051227264 -30.1496830110 -38.5068543819 -35.2435250892 -65.6542887376 -38.2700931830 -61.6057791954 -22.9942217289 -23.9091306901 -49.7599537136 -18.9394410633 -39.2695098205 -62.5378069205 -29.4770892222  -7.4076357401 -41.3722452049 -12.3860305971 -21.9537269564 -24.6767863017 -22.7988297841 -32.7318757237 -23.0458314997 -40.3374581531 -11.2381164386 -16.4099095774 -19.5742139343 -30.0437321711 -34.4327980027 -55.3831764732 -79.1168083611 -30.9649782220 -26.0292861392 -44.7008403169 -13.3584386184 -61.6250296736 -48.7751888233 -46.4430321992 -39.3618516204 -47.6452834098 -21.4645367133 -29.9606532297 -39.4346363770 -34.3609295246 -49.7189208111 -12.9974180323 -21.6327440120   0.2030972504 
			 38 	 -9.4302037219 -35.0752054730 -50.5898597153 -36.7540739885 -11.2623503972 -15.7863830391 -59.8457954323 -32.4935669562 -45.7409893212 -49.5113753571 -33.0371345991 -35.7249773779 -35.1613649146 -37.9911245090 -10.5320556742 -32.4859510994 -24.8247142940 -20.5255896980 -31.0235157092  -3.5504087612 -45.8366412897 -22.8012961444 -35.8761706572 -47.6676938781 -22.4980175137 -41.9799193699 -55.5892420049  -7.6763702762 -15.3550679425  -2.7090579556 -31.6726553541 -39.6902800195  -4.1439418631 -28.1514304806 -36.0591310663 -45.9022297950 -37.3635738952 -39.8682355672 -35.2384206400  -0.8103856350 -24.3945392400 -13.5826066936 -20.7232274868 -29.6676844865  -8.6073194903 -42.5827431574 -48.1055627952 -26.3028320647 -33.0413349352 -41.5846854441 
			 39 	-40.2992187476 -54.8782076973  -5.2991427170 -19.0347615984 -37.9701991087 -40.4644209648 -51.7537836561 -49.7551171567 -20.4676402377 -33.6718499997 -34.1425256582 -74.6986142205 -21.1493368527 -17.2007226290 -22.2517771683 -17.7701418254 -30.9276383876 -35.3905720177 -61.5387845957 -12.9885749901 -49.1480881159 -11.2909399158 -54.8051986496 -21.5318535751 -42.9194407227 -38.9059390037 -64.0797953844 -21.2961091992 -22.3440663503 -31.2613576535 -28.4872480135 -14.2920731139 -18.3413597790 -59.7433019494 -35.6679503903 -40.4243938233 -53.1535986925 -33.8113918174 -40.3706824241 -41.3754314301 -28.2140823201 -24.7349104916 -38.7242686916 -15.2860915670 -32.3973584280 -53.5766359913 -35.0642114531 -30.5308418289 -30.1779615275 -44.0517825162 
			 40 	 -7.2809901389 -30.0772604609 -39.5912715701 -30.1849158924 -30.0611872974 -49.1531159320 -45.2428290777 -71.4078635940 -40.1483582404 -44.9621193356 -24.8959272182 -11.9746322809 -20.0252975468 -40.0370870496 -30.2817106299 -34.5028009796 -28.6188549853 -25.1017316928 -34.2674026716 -49.0970681398 -31.7079730879 -23.6976756333 -38.6424856028 -27.4213428370 -46.7930692374 -36.1921559709 -24.8008402188  -3.3609536573 -39.1989598307 -26.8721805353 -45.5548192147 -58.1492757888 -31.2205838490 -22.2621703800 -51.0047539051 -27.5030594924 -30.5995862049 -28.9837433930 -11.3644453828 -31.2881888849 -18.6972122844 -26.0209879448 -30.4736227976 -40.9193278389  -4.8495895525 -28.9798545541 -54.1013312410 -26.8861225347 -29.5003264932  -9.7995502850 
			 41 	-49.5483151373 -20.1615730415 -40.4249418341 -44.7919052372 -34.9839369096 -30.2874050748 -41.6576875153 -16.6409886716 -27.9784418925 -19.8746949430 -38.8114695364 -15.6434076993  -6.7336370659 -48.7012432347 -35.4749752100 -49.3116714253 -10.8821550059 -39.6914061937 -27.7069541681 -50.5931000122 -19.1336498653 -31.9882441233 -61.1366542227 -40.0368933584 -15.7272088526 -35.7704430236 -45.7950387708 -21.6295046933   3.9199359705 -43.4969057833 -26.6545017719 -43.2298364527 -40.3627606468 -27.4340631308 -44.8447163639 -67.1683264082 -17.9505386557 -43.1123228461 -19.7970809563 -15.3048558611 -44.0408821999 -25.9317854135 -46.1852787882 -24.2364464233 -58.5480201557 -22.5916150438 -27.3132023595 -15.7747843679 -18.6835678062 -12.1803595776 
			 42 	-40.1281323854 -36.5222389102 -50.3977780670 -31.6874728658 -58.3214390010 -23.1033544228 -44.6080582069 -40.9613184119 -29.7222185993 -21.8945835879 -32.7382540925 -39.5136874985 -56.6955905537 -32.3420303095 -13.7460063812 -41.3215214045 -71.4119885587 -53.4187944800 -41.7227013127 -45.2253790211 -30.1541306487 -31.8747715316 -24.3801618575 -36.8503019798   1.0465231553 -36.9256019771 -24.7708260887 -22.1450802647 -63.2460905850 -23.1816443322 -32.5872585910 -36.7535612044 -36.9404464630 -38.8609944354 -38.0460913467 -34.2796305961 -48.3539165483 -18.1126194423 -27.5712469352 -34.2379694097 -32.9323386860 -21.6120521138 -29.3123789744 -62.7570473168 -46.1818753456 -40.9230909648 -44.4746043973 -40.0110797476 -15.6212398866 -19.9984781027 
			 43 	-42.4682305610 -60.1698924841 -34.4902958404 -50.2214658064 -57.0546098489 -32.5396666938 -30.6171609055 -39.8571079735 -36.3788278919 -37.3104074547  -3.6420633285  -8.3884700858 -42.8271380394 -21.7239973482 -24.4159317435 -42.5804292682 -39.3254457145 -33.5974150500 -23.3102982780 -47.8417389166 -39.6634774691 -26.6605832513 -36.1025369750 -22.3540670719 -34.8334243183 -69.6648028890 -10.7430796467 -37.2280490175 -36.3512521053 -38.0828999183 -35.0308064476 -36.9830364931 -59.3896511983 -37.1916985491 -35.6283736047 -18.2230665713 -50.4004516212 -37.6366862215 -11.0509888247 -23.3970187454 -16.1478362292 -16.0417579118 -22.9614379180 -34.1759349330 -51.1936152571 -28.3802785717 -34.1560785129 -33.5160951270 -20.6889024899 -41.8005151828 
			 44 	-37.6596520536 -52.8599791913 -28.2615458189 -26.7951456942 -16.3134611450 -43.5048576282 -35.1562641134 -16.8471351018 -40.7607942470 -17.5207956700 -14.2714429779 -49.8158726539 -50.2839314375 -56.5942091807 -24.7761541085 -20.2958847588 -53.6490092527 -21.8779772921 -22.5925861770 -20.2868222560 -53.5678561417 -34.8877263600 -55.9729105908 -42.0465826980 -40.1768440241 -25.4561977428 -38.0931128853 -22.2435003596 -41.9183170537 -46.9356357323 -37.2435614252 -38.5985769870 -37.1641781292 -12.7883435835 -35.7079647855 -28.5603679145 -43.3221763601 -42.4357858332 -52.1657469017 -25.3495074098 -39.9341630114 -54.9910587357 -29.0077650469 -27.5115747774 -37.4535792028 -54.3793069675 -29.0677631694 -28.2284751645 -26.3239150715 -21.5243815998 
			 45 	-13.1232951766 -25.9412822144  -7.2190735803 -41.6927473545 -23.1211273604 -53.9578662525 -55.5087514976 -34.9835264797 -40.3803259529 -22.9743022664 -20.9048099586 -49.0992085144  -6.0975870661   2.1082226132 -57.6320348142 -25.4614052214 -29.2071459500 -29.0903292848 -33.2262097963 -60.7304916618 -42.5556755831 -19.1540948324 -26.5313498186 -23.2961829512 -38.2638763111 -30.2183817650 -54.1096179002 -16.0133236106 -12.2251846508 -26.6875457346 -43.6467336727 -28.6439570613 -11.7640321619 -26.8968184820 -23.7301271499 -37.7637519203 -60.3581377143 -26.3045908503 -34.0540916317 -25.8118962829 -13.2424821845 -46.9676431144 -33.7116805977 -55.8983303402 -68.4825875902 -22.4526636467 -34.1928433499 -41.5748257577 -28.2088532110 -27.4689380325 
			 46 	-49.9016755917 -53.5696407945 -19.1378473518 -48.9533715902 -61.9985530878 -11.5889598568 -17.4078493432 -25.1451051985 -44.1770092190 -10.9377708305 -35.7201131119 -31.2631771422 -39.0510847948 -33.6953671716 -57.3811467893 -43.4400221098 -51.6347907630 -51.7322419149 -58.5909566255 -28.7570328710 -24.2043101123 -26.2275113060  -1.4044731777 -37.6609388184 -35.0424801727 -64.8245564924 -38.3511524211 -24.8296529319 -30.8497856646 -41.5118839462 -14.1990051292 -10.6394435834 -36.3645538684 -39.7115037557 -36.6908579347 -34.3946927472 -32.6860129263 -24.8841419043  -5.2263667452 -36.1133529524 -44.3360581257 -38.3848233539 -10.5920264137 -29.5397410677 -38.5996369873 -38.1090035945 -26.2879874577 -32.2606800932 -37.8493086517 -31.6685987022 
			 47 	-24.2020552070 -31.3611764214 -55.0522338295 -49.3572707313 -20.9897370172 -38.3580590860 -22.0599382743 -45.1509132169 -30.9525288587 -45.2967705910 -50.4749736225 -22.3252210919 -32.5143736816 -32.0128849693 -22.2223831100 -58.5117019375 -14.6042329384 -32.7026929188 -44.7779614937 -30.3876994330 -38.6136965893 -36.6196107034 -29.9929112472 -55.8759826049 -56.1569564094 -47.5515256141 -44.3312515223 -35.1580751647 -37.0631190553 -31.4603758367 -38.6727919344 -42.6992591793 -24.7998072035 -38.6924981522 -50.0174688687 -20.6437979022 -54.5902980004 -42.4201426185 -55.6384194231 -36.5388376642 -31.2146560731 -37.5283140112 -32.4086968640 -33.6019578711 -35.4439059561 -32.6835072800 -38.5623606165 -26.0817701589 -49.2193086132 -48.4684796649 
			 48 	-54.4961174182 -24.2059608875 -59.5629201374 -41.0665281086 -24.4474545472 -37.1439489941  11.0194728353 -25.4881265668 -19.5273604007 -28.1411373574 -60.2179551767 -18.8884425122 -41.0308445675 -45.3061757452 -50.5939742360  -7.5758757932  -8.9696600397 -52.2950318323  -0.4979875028 -11.8889833134 -28.8299041133 -34.1788355770 -58.0372610715 -55.9790338380 -32.2173634925 -36.0507528950 -48.2053175507 -28.4292572043 -30.0786282543 -21.9974255777 -35.8059052424 -41.4578748674 -39.0060805856 -16.4489896982 -27.8072729799 -19.1603145503 -32.7600189881 -30.0261139851 -39.0029770988 -27.5140483468 -31.9984614095 -25.9167053038 -38.4865969120 -31.4940666584 -24.8310369879 -16.2742573005 -16.8162737290 -39.4688767418 -49.4517107810  -7.3669980786 
			 49 	-49.3486108307 -13.4453886374 -42.1464985568 -41.9613081129   3.5239293948 -24.1645681024 -54.5839067674 -36.9369615314 -43.9558650597 -19.4036007883 -40.4193386533 -46.4885133402 -43.3719800037 -37.2304287514 -24.8788034630 -42.5957014401 -26.7916333312 -37.4543446753 -30.6128162322 -44.1429474271 -42.0196122656 -33.5469268737 -16.6170767082 -39.0590168390 -25.8366810740 -27.6545462466 -19.4926432045 -51.2379305513 -24.8374361345 -36.4203062539 -58.8978644686 -19.5495552746   3.3139830199 -15.1843188390 -32.1982334105  -1.3090023863 -26.0019248610  -8.8995736175 -43.6500010851 -26.9735483977 -51.6237954799 -44.9073652860 -45.6589483574 -35.2286636724 -40.4406438685 -19.7124996322 -22.8098498215 -32.8284459616 -48.2566518937 -28.3472201887 
			 50 	-36.4475172259 -31.2647993369 -49.9024826799 -34.6013881072 -13.9752100124 -36.6045937996  -5.0839078603 -37.6798210334 -33.8008363199 -55.3853496022 -40.6343508760 -39.6919717079 -25.1298577990 -35.0824978690 -25.7918272193 -46.0958789490 -35.8860667055 -21.3421946620 -42.5277380813 -38.3176220432 -30.1314815805 -35.6010971611 -28.9547149123 -42.0313418003  -5.7366856147 -28.5408383812 -24.4481510913 -32.9947145976 -32.7356585264 -55.7034606148 -29.2789430996 -17.8952385787 -35.4709346889 -44.5117747752 -71.1249456163 -59.4088596132 -33.9148803866 -44.3203576484 -26.2785824911 -27.5489048682 -30.8659941544  -8.2544287208 -51.3280927306 -56.3298361607 -40.3551213744 -57.4649068152 -44.0046555145 -27.3894293236 -49.2644881800 -35.9210271733 ];
		Beta_FEE=Beta_FEE[:,2:end];
	#Variables choice model
	  AT_FSP = 10;
	  TD_FSP = 10;
	  AT_PSP = 10;
	  TD_PSP = 10;
	  AT_PUP = 5;
	  TD_PUP = 10;
	  Origin =[	1	0
	2	1
	3	1
	4	0
	5	0
	6	1
	7	0
	8	0
	9	1
	10	0
	11	1
	12	1
	13	0
	14	0
	15	0
	16	1
	17	0
	18	0
	19	0
	20	1
	21	0
	22	0
	23	0
	24	0
	25	0
	26	0
	27	0
	28	0
	29	1
	30	1
	31	0
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	1
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Origin=Origin[:,2];
	  Age_veh =[	1	0
	2	0
	3	0
	4	1
	5	0
	6	0
	7	1
	8	0
	9	0
	10	0
	11	0
	12	0
	13	0
	14	1
	15	0
	16	1
	17	1
	18	0
	19	1
	20	1
	21	1
	22	1
	23	0
	24	0
	25	1
	26	0
	27	0
	28	0
	29	0
	30	0
	31	1
	32	0
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	1
	40	0
	41	1
	42	0
	43	1
	44	1
	45	0
	46	1
	47	0
	48	0
	49	0
	50	0];
	Age_veh=Age_veh[:,2];
	  Low_inc =[	1	1
	2	1
	3	1
	4	1
	5	1
	6	0
	7	1
	8	1
	9	1
	10	1
	11	1
	12	0
	13	1
	14	0
	15	1
	16	1
	17	1
	18	1
	19	1
	20	0
	21	0
	22	1
	23	1
	24	1
	25	0
	26	1
	27	1
	28	1
	29	1
	30	1
	31	1
	32	1
	33	0
	34	1
	35	0
	36	1
	37	0
	38	1
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	0
	49	0
	50	0];
	Low_inc=Low_inc[:,2];
	  Res =[	1	1
	2	1
	3	1
	4	0
	5	1
	6	1
	7	0
	8	0
	9	1
	10	1
	11	1
	12	1
	13	0
	14	0
	15	1
	16	1
	17	0
	18	0
	19	0
	20	1
	21	1
	22	0
	23	0
	24	0
	25	0
	26	1
	27	0
	28	0
	29	1
	30	1
	31	1
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Res=Res[:,2];
	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	for n=1:N
	  	for r=1:R
	  		q_parameter[1,n,r] = Beta_AT[n,r] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
	  		q_parameter[2,n,r] = ASC_PSP + Beta_AT[n,r] * AT_PSP +  Beta_TD * TD_PSP;
			q_parameter[3,n,r] = ASC_PUP + Beta_AT[n,r] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n];	
	  	end 
	end
	  
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R);
	for n=1:N
	  	for r=1:R
	  		Beta_parameter[1,n,r]=0;
	  		Beta_parameter[2,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
	  		Beta_parameter[3,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	 	end 
	end
 
    return Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p;
end

function Mixed_Logit_n50_r100()
	NUM_POINTS=2+1; #alternatives
	N=50; #customers
	R=100;
	# UB_p=1*ones(NUM_POINTS);
	# LB_p=0*ones(NUM_POINTS);
	UB_p=[0
	1
	1];
	LB_p=[0
	0
	0];
	#Parameters choice model
	#Parameters choice model
	  ASC_PSP = 32;
	 ASC_PUP = 34;
	 Beta_TD = -0.612;
	 Beta_Origin = -5.762;
	 Beta_Age_Veh = 4.037;
	 Beta_FEE_INC_PSP = -10.995;
	 Beta_FEE_RES_PSP = -11.440;
	 Beta_FEE_INC_PUP = -13.729;
	 Beta_FEE_RES_PUP = -10.668;

	 Beta_AT=[	 1 	  1.7752919700   1.4650269027   0.0246131522  -1.4301606721  -1.9055857350   0.4811625732  -1.5978920933  -0.9651872802  -0.5067220717  -1.3048485394  -1.5127537899  -0.5814813489  -2.6579411756  -3.1893796161  -0.6533916966  -1.1148105685   0.0055555855  -1.1066780617  -1.6229402151  -1.1715269464  -1.4723275899  -1.8340293914   1.9454925107   0.5812019402  -0.5501841797  -0.5313047386  -0.8642175548  -0.8537063677  -0.0926835375  -1.3629332996  -0.3561510457  -1.7664966462  -0.3654385046   0.0616334755   0.1242572115   0.8331292928  -2.2002131724  -2.1341958953  -1.2971951785  -0.6134864239  -0.2157006774  -1.1943214425  -3.5335042040  -2.0690310300  -3.6442945260  -1.4019047195   0.7678352308  -0.0294795928  -0.4169210716  -1.5672660848  -0.2451172250   0.1605868883   0.1013684668   1.2731219223  -1.4245499889  -0.9693240216  -1.0533766655  -0.3585142487  -0.1295563115   0.4972655449  -0.8964660409  -2.1686306903  -0.3725119604   0.3765444306   0.9171494877   0.4988318274  -0.5510708754  -1.1657806974  -1.9551330714  -0.4234021364  -2.0563753464  -2.1183319355  -1.2342298012  -1.6403775564  -0.0139827668   0.7859428961  -1.1961477774   1.2762088621  -0.0942329609  -0.0265036761   0.4936647480   0.4160713321   0.1258442404  -1.3513304681   0.0254556422   0.5123888677  -1.0634230574  -1.3416243614  -0.1229679054  -1.5531604753   0.1902170772  -2.2480981097  -2.2352401134   1.3511899620  -0.1921716570   0.8967216230  -2.6486258751  -0.8323448165   0.6548641715  -2.6135568382 
		 2 	 -0.1665611017  -1.6792303905  -0.1236847389  -2.3799194309  -1.2430896211  -1.7689909124  -0.6670147193  -1.9451424885  -1.3226361535  -0.1751176995   0.7638447015   0.1125907520   1.3489554921  -1.4705915789  -0.7267716045   0.2389045489  -1.1208168594  -2.1670280353  -2.1483861394  -0.4008369918  -1.4208665673   0.9362745028  -1.4093370776  -1.4356176850  -1.9945021672  -0.3189490941  -1.5102201010  -0.4382108725  -0.6649443221  -2.0919439482  -0.3220882640   0.7972679442  -0.7277302965  -0.9792136528  -1.2710414168   1.1115544679   0.1642111785  -2.1568929396   0.7059714611  -2.1834196497   0.2409163200  -4.2118676561  -2.4988400177   0.1831655157  -1.3699671248  -0.8228961423  -1.5241593981  -1.1013011478  -1.7958589579  -1.7764385280  -1.1192945619  -0.8315663295  -1.3049616150  -1.1613304222  -0.7217126187  -1.6699178785  -1.0391277559  -0.4285531459  -1.2458788974  -0.6199415186  -0.0360911423  -0.2962641928  -0.6536684025  -1.8984329315   0.6780222840  -0.1117856382  -0.2128285295  -1.1400686721  -0.2841327631   0.6222118608   0.4043549850  -2.5455980652  -0.5013802765  -0.3441319772  -2.3166426494  -2.0016449963  -2.0762140012   1.0111585610  -0.5377160211  -0.8772527173   0.0785511041   0.3595550366  -1.2552407508  -0.7423693376  -2.1758541228  -1.7559255013  -2.3634817711   1.0205038225   1.3606015841  -0.1180203364  -1.1374768406  -4.6208891990  -0.2466317214  -0.7206669746  -1.8183382624  -1.8287826072  -1.0627914310  -1.5618543914  -1.1480963760  -0.4785176501 
		 3 	  0.6889694468   1.5885689528   1.6449202530   0.2884371091  -1.1577524272  -1.0131528300   0.0335618538  -1.0275786944  -0.9601722341  -2.0795049934  -0.9896731755  -0.3849807534  -0.0335136930  -1.6855887061  -1.7695242102   0.4555265264  -1.4448660866  -0.0095590895  -0.6669754517  -2.0435981540  -0.3919504602   0.6903998190  -1.7777312276  -1.4012652084  -1.8801321118   0.1293788526  -1.9962316428   0.8572504375  -1.3042726606   0.1982080273  -2.6172044800  -2.3638455568  -0.8886806977  -0.6542966972  -1.0238534807  -1.9814176707   0.3884522979   0.7279630148   0.3542958232  -0.2382522848  -1.8646471060  -2.5783214127  -1.7461384246  -0.3562310193  -1.1101991563  -1.3384504207  -2.6652222198  -0.8735087364  -0.7356511993  -1.6276397197   0.1097008463  -1.5284735531  -2.0614179884  -0.9517122921  -1.9199742579  -2.0401324041  -2.4611799486  -0.8176587786  -0.1392086586  -1.9272113719  -1.6096728673  -0.9858799518  -0.3547579469  -0.7234327764   0.0776833584  -0.2838255648  -1.3458696674  -0.1647702075  -0.8082893808  -1.0647289139  -0.9551300918  -1.2360390693   1.4077559320  -1.5095409834  -2.0336567206  -1.7937079357  -1.4646691800  -0.4320471994  -2.5754443818  -1.5918001049  -1.8862270225  -2.3786403104  -0.7324783310  -0.2412424776  -0.6379367289  -0.9054490914  -1.6351122265  -2.7170955063   0.7991595056  -1.3564901042   0.8050672909  -1.3529906802  -1.2875314742  -1.6378493260  -1.2220642982  -1.8719979213   1.1669495012  -2.8070147543   0.0405936018  -1.2554174281 
		 4 	  0.5590778661   0.3595847598   0.3765831166  -2.4194850036  -1.8020266831  -0.1592448350   1.0953031787  -0.3686176288  -0.9356236372  -0.6280491488  -1.1552836720   0.2277207231  -1.1320004829  -0.6860848098  -1.1065327994  -0.1387152087   1.0007857243  -2.5587837443  -0.1248228572  -0.4890150798  -1.0811116913  -0.4178553072  -1.0692012093   2.5951898105  -0.1814256298  -1.0654009970  -0.2510456890  -2.2947189326  -1.5541480049   0.0160581483  -0.4123006145   0.4402898049  -0.0095874828  -0.2023419299  -0.9013967640  -4.5605991943  -1.1721071731  -1.3221300047   0.3401858415   0.4975260101   0.5293347305  -1.3512144190  -0.7945607306  -0.1784784030  -2.7177505697  -0.3245931564  -0.1786077820  -0.5734441279  -1.6034070601  -2.2873378601  -2.3059616795  -1.6794851094  -1.1967554685  -1.9407116877   0.1961749287  -0.9186520115  -1.2122821728  -0.8181212657   0.2233000226  -1.2665517040   0.9858232296  -0.8685051639  -1.6696269124  -1.8685702423  -0.6258377887  -1.1267735968   0.4201937685  -1.3997269337  -0.4568641517  -1.2429026826   0.0505116574  -1.2530948254   0.5949382256  -0.2733812696   0.7746972495  -0.8914382932  -0.6887455476  -0.3450680764  -0.2167868395  -2.7355124602  -0.8741069907  -0.5919773420  -0.7024416259  -3.1345053244  -1.4661644472  -1.4933479458  -1.8371895219  -1.8102786043  -0.0850194206   0.2741951858   0.3151853144  -1.5104781728  -0.3005495367  -1.7810039078   0.3282746446  -0.0507886086  -1.3176181112  -1.9643821154   0.1857245935  -3.0303702837 
		 5 	 -0.7885303739   0.6498573241   0.9489960754  -0.4027836142  -1.8778870376  -0.2047089093  -0.8796792193   0.7622536241  -1.2158627485  -1.3075029519  -0.3218672293   0.0129066456  -0.6060707830   0.4341161177  -0.1631887278  -1.3774196438  -1.9370720513   0.1095748683  -0.5092796281  -1.5226422529  -0.7943588974   0.0931198800  -0.0679072021  -1.1725023298   0.0294549024  -1.0070827952  -1.1691878370  -1.1028013940  -2.9406927770  -2.4703893757   0.3755206615  -2.0538792561  -1.0087090547  -1.9177901427  -1.2556343204  -0.2023116695   0.0235903596   0.0946551383  -1.1420066995   0.2768100182  -0.5159602072  -2.1199023989   0.0220029228  -1.5080039233   2.2280612600  -1.7530286403  -0.8768950442  -1.8583576084  -2.9101968965   0.0292996733   1.0275150333   1.6299793111  -1.0731942601  -2.1115199647  -0.9207038033  -3.5994764662   0.1950533332  -0.6573241658  -0.8831139979   0.1394842294   0.9339445339  -0.6687788306  -0.8112511662   0.0228340630  -0.2074714973  -2.1372078312  -1.1085309275  -0.4827917325  -1.3754623011  -1.5474573235  -1.7226635897  -1.0162366134  -1.2154966625  -0.1889765652  -0.8468398222  -1.1318486957  -2.3716148667  -1.4432087097   1.4294819091  -1.3093550627  -1.5752813939  -1.0787219971   0.2864722695  -1.8929110594  -1.9319618433   0.0844367663  -1.6779979141  -1.7725341761  -2.4707843350  -1.7699921099  -1.1632166572   0.2048022632   0.1186391164  -4.1173628679   0.7237947399  -1.2453098083   0.6164500045   0.2815759336  -0.0840904893  -2.1999299069 
		 6 	 -1.6320933259  -1.4545303070  -0.4761579025  -1.6567037678  -0.3283903312   0.2855053268  -1.1501849150  -0.3782638848  -2.3810387369  -0.7543513767  -3.5032208043  -1.8909635924  -0.5835509140  -0.3398707492   0.7947946473  -2.9385430120  -1.8802671386  -1.8160754053  -2.5617416876  -0.4706089738  -0.0978079947  -0.1085941187  -1.6876389814  -1.1377141208  -1.5197592375  -2.4987074045  -2.0630743225  -1.9360413987  -1.0508312976  -0.9566025545  -1.9720093857   0.1692597149  -1.0544345647  -0.9445459696  -0.0762642050  -0.2215676913  -0.7062813625  -1.5247054278  -1.4446675646  -3.4900486940  -1.6195453403  -0.6460365146   0.2629523832   0.2527658918   0.7907660384  -2.9860152507  -0.1937885842  -1.1043812148  -1.7341513241  -2.3970385784  -0.3830384147  -0.3037291962  -0.8653280696  -0.1973129671  -1.1278553818  -3.1777346488   0.2168206263  -1.2415358615  -2.2010723938  -0.4302570742  -0.4636423115  -0.2491351000  -2.4571994624  -1.6678196467   1.2113051797  -0.5138193757  -1.3386170036   1.1033106905  -0.2364793227  -3.8004940816  -1.6097605450  -0.0356264703  -0.5379548888  -0.1793226482  -0.6855343195   0.4135970878  -2.3451694270  -3.0773380861  -0.8834614898  -0.3070112492  -0.1758323259  -1.7133812814   0.4178665319   0.3471647778  -0.6689208023   0.2397970023  -1.2047361673  -3.0293441492  -1.0272561872   1.6031115538  -2.2992198373   0.6048174834  -0.9023261720  -0.0340465447  -0.3392526171   0.3683186492  -1.6975166369  -0.7428596882  -0.4876923007  -0.1097619176 
		 7 	  0.5917818306   0.5733252361  -2.1056962429   0.6960841221  -0.0612020458   0.7633195074  -0.8788769258   0.3872741173  -0.7379460637  -1.0936625615   0.1086016599  -3.4301440051  -1.0097882778   0.1334595426  -2.0040164705   0.1763158959  -0.2755101461  -0.5263802822  -0.1756083558  -0.3149389923  -1.4746989431  -0.6303290134  -0.7854657503  -0.0578076988  -1.4001741526  -2.5001987481   0.1943662779  -2.4243819002  -1.7766010939  -0.7359353997  -1.5778662889  -0.2094274834  -0.1547430710  -4.1757443045  -0.6419073437   1.6061187333  -1.6043280557  -1.9370428681   1.5473843882  -1.6036832328  -0.8704604252  -2.2718518292  -0.2298414060  -1.1450667018  -1.7503743946  -2.6355051241   1.0236855380   0.1912705775  -1.5088043279  -1.3620228205  -1.5129930552  -0.4676593946  -2.7179111658  -0.1208474378  -1.0002372093  -2.4928868455  -2.0644579213  -1.4256459569  -0.6785075843  -2.4303250042  -1.0674478582  -2.0307491583  -2.2134908847   0.7080767650   0.1205157978  -0.6380274390  -0.0518730102  -2.1710165413   0.7940416871  -0.4005047285  -0.9808038552  -1.1581273422  -2.4827494102  -0.8236820075   0.0241611877  -2.5785083654  -1.1929229655  -2.9747524027  -0.6073576774   0.7243692627  -0.6496874730  -0.1845160369  -0.9818429612   0.6365501213  -2.8190665347   0.9514363417  -1.5346773226  -1.6815859228  -1.4760654718   0.6668385246  -2.5156625976   0.5219384830  -0.0039509719  -1.1066442974  -0.7280774045  -0.7898217258  -0.8511492331  -2.2399775445  -1.9900373872  -2.5815220749 
		 8 	 -2.2074138653   0.8488946249  -1.2316907759  -1.3355398242   1.3928852811   0.4820472821  -1.7267526276  -1.7469104660  -0.8855071746   0.2146236122  -1.7343902343  -0.7589041264  -0.3632467878  -1.0070916273   0.3921484654  -1.0139577428  -1.2484987676  -1.4274828392  -0.9855757278  -0.0231631506  -0.7427110990   0.7380611175   0.2683877902  -0.5650748531  -2.3206719980  -2.5179245165  -0.9462225976   1.1527862474  -2.4401358832  -0.4100885147  -0.0007437788   0.2092049431  -1.8081673903  -2.6575744087  -0.4807521997  -0.0421166085  -1.2243626242   0.7018668811  -0.4082576472   0.7826998469  -1.9949706935   0.4931499091  -2.9277771285  -0.6874292410  -2.7802744576  -0.7911417599  -0.9657509730  -0.8206418894  -0.3590096670  -0.6571199717  -0.3795864869  -0.7684463563   0.6760811487   0.4328604488   1.3535012405  -0.7372856990   0.5179782604  -0.3606829192  -1.2023732117  -1.1495850901  -0.7217732629  -0.7710234435  -0.8388369984  -1.4311022322  -1.5952556529   0.3305977948  -0.6553608337  -0.6054640197   0.9869355498  -1.7036497164  -0.4482715609  -1.2689835300  -1.9302122547   0.1592930123  -1.8660575299  -0.9294464987   0.9279513820   0.5185960308   0.0294771690  -1.6781948340  -1.6104043365  -2.2144794119  -1.4097806561  -1.1043023705  -1.2127022221  -0.1562140440  -1.0276195869  -1.1033656158  -1.9813883548  -0.3432793433   0.8479029378  -2.0949702263  -1.2380722182  -1.3518897359  -2.5675460233   0.5738122151  -1.5250466427   0.9717825135  -1.3916355733  -1.1853802391 
		 9 	  2.0599268681   0.0547640282  -1.5968547897  -1.4533655973   0.4710850475   0.0913059464  -1.3671251691  -0.1066192443   0.6416264843   1.9836986922  -2.4249808921  -0.4575944450  -1.0698397769   0.4244056106  -0.5444435340  -1.8883449015  -2.2841167365   1.3421973766  -0.2199070210   1.3313314315   0.3367238596  -1.5277158951  -1.2777818181  -2.9521154835  -1.9033000693  -0.7862907607  -1.1714835181  -2.8164609785  -0.8014385689  -0.1956694590  -0.0315898827  -0.4279779895  -2.2911167592   0.5932958540  -1.5391753947  -1.8744484498  -0.3531776245  -1.6690638619  -1.1067484025  -0.9421849219  -0.5801559216  -1.5009641437  -0.4234816243  -1.3383602469  -1.4812107523  -1.9689120550   0.3369700879  -2.8521879503   0.5261364473  -0.0908430457  -0.3214895717  -1.1054976459  -0.5190769644   0.9046144739   0.5406084670  -0.5012504442  -0.0826438783  -1.6194357456  -1.0372724285  -0.8261067640  -0.9203787286   0.1909185238  -2.2462759044  -2.1379863508  -2.2370107835  -1.3131470638  -0.8604246019   0.0863119655  -0.2642364936  -1.3562508022  -2.0910764973  -1.2393216995  -0.1585177117  -1.5215752637  -1.3987538320  -2.1445044504   0.3759747135  -1.5331321191  -1.5095535313  -1.0435621833  -1.3060300420  -0.3167840443  -2.0913888296   0.5124112566   0.1559086651  -0.8912612633   0.2472716772   0.2407223853  -0.9138652879  -2.1133179289   0.6201015425  -3.0267299997   1.1344825586  -1.3380608873   1.1890255813  -0.3991355929  -1.3509181326   1.0581638462   0.3227738471  -1.8541727284 
		 10 	 -2.9716083135  -0.5863210032  -1.3996176863   2.6147484168  -0.3861586781  -0.6273011469  -2.6578417083  -1.6828129205  -2.6869731341   0.2872454657  -1.0680466242   1.0234539901   0.5216863355  -1.3127794995  -0.3405630902   1.2497617391  -1.5469572070  -2.0364664517  -0.2189851183  -1.1212163354  -1.1348998194  -0.6868545595  -1.5647585826   0.3077606058  -0.7180592638  -0.9018081032   0.0531925086  -0.0684613726  -2.0864647519  -0.0673076004  -0.7334032580  -1.0510076977  -1.5816044955   0.5004168781  -1.2860785091  -0.7086789000  -0.7783581667  -0.4385369753   0.1149128976   0.1735982834  -0.2775915910  -0.7206893574   0.2911327887  -0.6343687896   0.1061521717  -1.6914518765  -0.6163609504  -1.2816040617  -0.3144283952  -0.9867487748   0.9295231614   0.6110856837  -0.1402787559  -0.9628531714  -1.4390041682  -1.8036708457  -0.7501677623  -1.5103091995   0.5295313605  -0.0206212509  -0.3428291720  -0.3333546663  -1.9660230006  -1.5454485234  -0.7381801493  -2.6387057280  -2.2244297260  -0.8762071203  -0.6994752494  -1.1057561432  -0.7496500223   0.7223549239  -1.0276233906  -1.5463756038  -1.1458342519  -1.3942886458   1.2097644104  -0.2657757204   0.5377277073   0.1969350966   0.7405724826  -1.1915877971  -0.3916631494  -0.3594206331  -1.2431591699  -0.5078808841  -0.7030148888  -1.1954276293  -0.3512588138  -1.2202222363  -1.1431911202  -1.1105819371  -0.5832841806  -0.5796827988  -0.3799780117  -1.3822547556  -2.3018633534  -1.6043135320   0.3431194197  -2.5655949167 
		 11 	  0.8096359679   1.7069109678  -1.3024809423  -0.9342736092  -0.2745452061  -1.1056897550  -1.2698560344  -1.2850045470  -2.8885642954  -0.8262250924  -2.0118697137  -0.2429351021  -1.3615093559  -0.0374140076   0.8141190893  -0.5136992674  -2.7553303446   1.2150987587   2.1936521985  -4.1652650521  -2.1333377881  -0.1653235154   0.2477189776  -2.6087093757  -0.1587068403   0.5845149340  -0.1353673542   1.4921849357  -0.4242957035   0.2755617246  -1.8359901141   0.8167667602  -1.9202041384  -0.3036521168   0.5514758126  -1.4043628353  -2.0223122208  -2.2476409701   0.3362455893  -0.8868388187   0.2494535058  -1.0477395653  -2.2862683532  -0.9792206668  -0.3821686170  -1.1173694779   0.3760133009  -0.5221285629   0.1411503363   0.9320810262   0.9803773719   0.5537614291  -1.2699034451  -0.8218215507   0.4349156716  -1.4685451205  -1.2597269780  -1.3021492835  -0.9975550062   1.0961750913  -1.4721652184   1.0655510423   0.0371187911  -0.1830663716   1.1722755569  -0.9141910913   0.0299735842  -1.8594486881  -0.2110934131  -0.2517033876  -0.7954526547  -0.6691629175  -1.8974623893  -1.5784770174  -1.3590137004   0.4073152585  -1.0253868652  -0.5076842552  -0.4427310268  -1.9079117971  -0.7902262985  -1.0507158144   0.6699309698  -1.3094390401  -1.9493630644   0.0855817850  -1.6945870009  -1.9019879258   0.8656017207   0.0920757590  -1.1326929135  -0.8029330887   0.4344684909  -0.4400312541  -0.0417798083  -1.2564049329  -0.3809729464  -0.6576690455   0.1694582826   1.2741570549 
		 12 	 -0.8949880707  -0.5163651818   1.5554964814   0.8112434911   1.1025996083  -0.4151120663   0.2115766791  -0.8991223235   0.3073934502   1.6118673626  -0.7698273449   0.2085789561  -1.6259302372  -0.6249003123   1.3243227290  -2.3827159845  -2.0796774082  -0.9761415790  -1.0221029478  -1.4348486698  -1.9061521366   0.4192297892  -1.6243427237   0.6047613010  -0.0619271292  -0.9394133746   0.1652087367   0.9039469535  -1.6462803775  -2.4742720233  -0.6080848830  -0.8158662410   0.2460446908  -1.4491653130   0.5960430878  -1.0076615822  -1.9471209833  -0.7246908270  -2.2202458992  -0.0848819796  -0.7452947689  -0.1285299954  -0.0590107519   0.7295342324   0.4891128435  -2.5624422499  -0.8324635335   0.5324392513  -0.1434066716   0.8018421278  -0.4345922146  -0.8414316775   0.2530896234  -0.5303068451  -0.8298468262  -1.8858312350  -0.4997332709  -1.6184521654  -1.0240512295   0.8190747639   0.1356729257  -3.3889018403   1.8232705271  -0.9315780080  -3.5799389558   0.1024751314  -1.9389674892  -3.6010548037  -1.8384781774  -0.5279755299  -0.3918926782   1.4638028773   1.4119594132  -0.3778044349  -1.4822332292   0.7904130098  -0.7492756793  -1.5063241073   0.7424537316  -2.0835075991  -1.5110350436  -0.4211036445   1.5877219495   0.0944772647   0.4492554538  -0.3954732198  -0.6874139329   0.2650846622  -1.6393552424  -0.8708833165   0.5756947725  -1.3151600460  -0.8081059001  -0.0122483305  -1.4581406750  -0.4522246219  -1.9051129970   0.2045854315   0.1929435181  -0.5674857242 
		 13 	 -0.6961968385  -0.3854487273  -1.8142651993  -0.0880424420   0.4134482442  -1.3485922454  -0.1751082370  -1.5659055212   0.8302192284  -0.3413741136  -2.7055812394   0.3963968537  -1.5290118778  -0.3746304869  -0.9417481913  -1.4235189933  -2.5203459047   0.1312057517   0.3928742068  -1.2692053380  -0.7604366749  -0.6843215136  -2.7179680425  -1.5813865529   0.0883507542  -0.9930037433  -1.4861695731  -1.0085572640  -3.3633315728  -1.8900113351  -0.8893925244  -1.8669582895   0.4435795666  -2.1850698879   1.2099368261   1.6136296063  -1.5944002239  -2.4294431597  -1.5262760896  -2.8965712813  -0.5171188158  -0.4064378484  -2.5575127747   1.5106672878  -2.4801644760  -2.3167672285  -0.1371752120   0.0110784770  -0.6878762239  -2.8366083895   2.2392614432  -2.9022799567  -2.3254723615  -2.5105822640  -1.1041822104   0.0220503199  -0.5975943880  -0.6466651777   0.0914300426  -1.2801054329  -0.2479890014  -1.1225164809   0.6996132930   2.1631754643  -1.4902742349   0.0090420381   0.2953469012  -2.1580616375   0.6722214774  -1.0679713153  -1.8546633915  -2.0281452761  -1.1326382563  -1.5420160071  -1.2626102985  -1.8647478421   0.0982245764  -0.8108222535   0.1063575008  -1.6671113677   0.5408806778  -1.4332366587  -1.3138222435  -0.6421851632   0.7663243189  -0.7159713669  -0.4686531837   0.6074390495  -0.4379985560  -1.7758725479  -1.5202060753  -0.1874012286  -0.7731176620  -1.3051744750   0.5462664043  -1.9194831558  -1.7605397450   0.7201460357   1.2885103749  -0.1854561210 
		 14 	 -3.3463648925  -3.1954554176  -0.8533236408  -0.5549581285  -1.1955055333  -0.0380071446  -2.5000736946  -1.5152031629   0.2615880730  -1.8822429313  -0.0618344987  -2.5657854900   0.5261083993  -2.2568470253  -2.5747329554  -2.0128014417   0.8329458261  -1.3681011449  -1.9235130905  -0.0868643300  -1.6094411382   0.2517723220   0.8511702329   0.1459440216  -0.0458565479  -0.3985669752  -1.7920273563   0.0268927449  -0.7321957096   1.4647611344  -0.9779865769  -0.8742591088   0.1726092813  -2.9165848698   0.4987557477  -1.2509193690  -1.1307990972   0.7924080470   0.0318096535  -3.5549462634  -0.0691648525   0.4477410224   0.1224127703   0.1451538618  -0.8612572708   0.9249131633  -0.1172554337   0.1132542572  -1.2589195627  -1.5917008012  -1.6020238821  -0.5877496029  -1.9420626210  -0.6116073072  -1.9160386372   0.1454352980  -1.1708858776  -2.0347500190  -0.3889347710  -1.6706887150  -2.4164777476  -0.7059772156  -1.2685695412   0.0513886682  -1.9279745869  -3.5677596738  -0.3689955353  -0.5959488033  -1.3397918291  -1.2279533161  -0.6308717590   0.9482570191   0.2696193318  -2.5850761750  -1.1647264663  -0.6308782090  -0.0580446946  -1.0528322067  -0.9081574673  -0.2482348499  -0.7836584992  -0.5771879804  -1.2855676617  -1.5485198440  -0.0966912886  -1.2069585437   0.1272583002  -0.8329823677  -1.7146276345  -1.1418108798  -0.4453504073  -1.3205569979  -0.8383416924  -1.9995523114  -1.3984870286  -1.8466836202  -1.1265576963  -0.6318910351  -0.5946482681  -1.8692016778 
		 15 	  1.1435776009   1.0721582506  -0.5925450531  -2.2076452501  -1.1073389169   0.2256111376  -0.4176693843  -0.8143762689  -0.1708654973  -0.7466663598  -1.8065488147  -0.2133780008   0.9372442085  -1.4241780133  -0.9803232494  -1.4205274696  -0.2644753447  -2.1496195153  -0.1310877041  -2.2980068216   0.4933756144  -0.2412127704  -3.2277301540   0.2461330997  -0.2608288647  -0.6625135717  -2.0397932195  -0.0295656896   0.6454818779  -1.8669204150  -0.6107878697  -0.7679999375  -0.9873800778  -0.9329131173  -1.5453429718   0.3710062790  -0.9160078053  -2.0920691074  -0.3277800963   0.3066319390   0.6852874514  -1.0064099134  -1.4340809412  -2.1727016658  -2.2057514065  -3.1119420302   1.3810307806  -1.3806332956   0.4463618189  -3.0737693688  -0.4128699820   0.7672547037   0.6334934060   0.4083453793   0.2745681686  -0.5373299248   0.7152429498  -1.0225811038  -0.5418409650   0.4204071039  -1.0065134884  -2.1446173428   0.3522953923   0.3208034321  -1.5527154263  -1.9587883368   0.5673937987  -0.1690342125  -1.0487720439  -2.3988825709  -2.1296175619   0.4869491000  -0.1649768055   0.8801395243  -0.6388642495   0.1990252864   0.1407840105  -1.8385993352   0.5188028418  -0.4844517425  -0.9181232150  -2.2127527640   1.2342624813  -1.2275063626  -0.3289240294  -0.8010524592  -0.6573151001  -0.0139241488   0.2062523935   0.5733892644  -0.1273000664  -1.3073952332  -1.9802706185  -0.2407946699  -0.8851300718   0.2075530319  -1.9834994232  -1.5108528421   0.4096211156  -0.8943824120 
		 16 	 -0.4715016598  -1.3495729730   0.8800098029  -0.6330007137  -1.1755821962  -2.4078579542  -1.9494861050  -0.4650454973  -0.2904852982  -1.1235964586  -1.3033515643  -0.6897682549   0.3927432287  -0.2574453839  -2.1110918198  -2.8903983897  -0.2724913357   0.3279282454  -1.8096453090  -1.8274027806   0.9676750768  -0.4310218445  -0.6042054121  -1.0543292382   1.7747690398  -0.8548388021  -1.5700895545  -0.2360230454  -0.0927511838  -1.3151438085  -3.3456296056   0.8531770639  -1.8860297740   0.5995738523   0.0489217330  -1.2444289699  -1.4649888254  -1.6321983165  -0.1237056689  -1.3013134710  -1.9034698266  -0.3120251718  -1.0188518035  -0.2157253384  -1.1490870514  -1.2970850325  -0.7063966453  -1.9671079283  -0.1520225680  -0.6711203118  -1.2055154842  -2.3061344130  -0.0397200890  -0.3603963124   0.4191443402  -0.2768884885  -0.9280162018  -1.0032615659  -0.3526590218  -1.3531441267   0.3539471512   0.0867619995  -1.3106889440  -1.7711472901   1.2082550923  -0.9308406489   0.9787522217   0.1524969523  -0.0335921396   0.1288699570  -1.5959173855  -1.5783648282  -2.9779051919  -1.9299207194  -1.4790537599  -0.7217617419  -0.2327573903  -0.0861873087   0.1536380254  -1.1073965879  -0.1639612968  -1.3966316405   0.9850091756  -1.2345109268  -2.9882279251  -0.4793296991  -0.0475899389  -1.2178604668  -1.9493313172  -0.6972370033   0.0604318411  -1.5883172397  -2.2665775624  -1.6014204103  -0.6757552444   0.5094357825  -1.5542422198   0.0768526012  -1.7277076111  -0.1947478097 
		 17 	 -1.0674680256  -1.9813357668  -2.5692751866  -1.6747852282  -0.4783392089  -0.8148413552  -0.8098401177   1.9559287367   0.2703572571  -2.1068776663   0.3085286130  -0.2086677465   0.1405701597   0.8808918359  -2.0104767155  -2.1097879668  -2.3112368091   0.4434021116  -0.8018064806  -2.2841480303  -1.2039681097  -1.2332160988  -2.3064288542   1.8647652201  -0.2791617526  -0.3005494618  -1.9182646883   2.0338835106   0.1109683821  -0.3052837662  -1.6524942218   0.8054408574  -1.0258883718   0.0669374060  -1.9995749807  -0.0506135710  -0.7195790266  -1.0309550915  -0.1907909713  -0.1854136379   0.0189387809  -1.7707788171  -1.2106811472   0.6199755380  -0.3824957906   0.5814875994  -1.5626225002   0.3789128236  -0.1580224335  -0.9391311773  -0.2002525570  -0.9022919290  -1.0270245295  -1.0355698552  -0.8135103503  -0.3555033546  -2.5889833173  -1.6114686286  -0.0356264703  -2.3559085954   0.1080733607  -1.2213854037   0.9995197130  -0.2218636732  -1.0247774717   0.0666187381   0.4747095803  -0.7170169940  -0.9839044810  -1.3206136355  -1.5644661337  -1.8716079995  -0.4908342728   0.6629016469  -1.4142341784  -1.8237446592  -1.1975942306   1.1857863624   0.0936948901  -1.7048999266  -1.4306680103   0.9894599979   1.1932132593  -2.8776189551  -0.8075471739  -0.2566986527  -2.7489699030  -2.2180864842   0.0901299859  -0.8869674823  -0.9855862554  -1.8098790050  -0.5234997279  -0.6001346158  -1.0680775587  -1.9349018548  -1.7505233480  -1.6072851881  -0.7278833486  -1.9085053821 
		 18 	 -1.6778149761  -0.7392079433   0.2665329612  -0.2046374149  -0.5901465485   0.3825543202  -1.4361985164  -0.7819751695   0.3624620677  -0.9801680138  -1.4383695484  -2.2564249837  -1.1005756594   0.5936768409   0.0902016181  -2.3212479244  -0.3575400316  -0.8235842453  -1.3177402115   0.2313697188   0.0717460717  -1.9571073387   0.5424297140   0.3449854552  -0.4996592491  -2.1097765852  -3.2977963806  -2.2850345222   0.6822833844   0.2649232046  -0.8614890146  -0.4221246511  -0.3448320271  -0.2227435054  -1.9883729683  -2.3395195859  -1.0037964665   0.7428775619  -2.9350112168  -1.5788055619  -1.2447150336  -1.2468726616  -1.2915721842  -0.9210809777   0.0880527250  -0.4506016482  -1.4994782925  -0.4367862289  -2.2701336905  -1.5477068593  -0.1718325763   1.0565971223  -0.9723968994  -3.8580039661  -1.2603089022  -1.6732733146  -0.5845926956  -0.2158868715  -1.8721839071  -1.3478459437  -0.6942255118   0.2849515814  -1.4343794131  -0.4095161940  -0.4130596049   0.1592579726  -0.0087662386  -1.8162061602  -0.0148603992   0.8646314767  -3.5697984449   0.7451959273  -1.3282273814  -0.2795366631  -1.5956443446   0.7501322770  -2.1694803156  -0.6914916866  -0.7975528619  -1.8957943395   0.2549122602  -1.3372530269  -0.9845197090   0.3012849122  -0.8865744161  -2.8440800157  -0.3916084024   0.1370847107  -1.5176429893  -0.5939115766   0.2534973007  -1.8826849708  -0.4434385646  -1.8501011087   0.2903615729   0.4457008707   0.7228136796  -1.0620765094  -0.4297535515  -0.8835128543 
		 19 	 -0.7330204068   0.5015383432  -1.6309805872  -1.8708175464   0.0654389657  -2.3678780719  -0.2404230545   0.7674263001  -2.0389509035  -1.6283118956  -1.3004788902  -1.3631969634  -1.9663740956  -0.1786719292  -0.7229907308  -0.7606472758   0.5178614644  -1.4793280983  -0.4502902537  -0.7826449518   1.2112874668  -2.8953281713  -1.4940685282  -3.0313846521  -2.0394551376  -0.1753544858  -1.4259609235   0.1451438580  -1.7404097968  -1.7462141579  -1.0734042202  -0.3414159886  -1.4612590760   0.3785767278  -1.9647436581  -0.9636544042  -1.2279168731   1.4500534533  -1.5036893339  -1.3000623739  -0.2421597136  -0.8715427058   0.4673344488  -1.1854428639  -0.7236879408  -0.7816491916  -2.3429160432   0.6707924355  -1.1524238753  -1.2115020449  -1.0217111356   0.0121889320  -0.4736713464  -1.1055487544  -0.1490335577   0.4783572433  -0.9741189821  -1.5865340022  -2.0721259308   0.3800297179  -0.0524312920  -1.2317789885   2.1934233458   0.2335096745  -0.5999459181  -0.5492904932  -0.8572892584  -0.6853601162  -2.1386711512  -0.8970774740  -0.6747677473   0.4364642039  -1.3952935343  -0.4704992977   0.0271305550  -1.9533987135  -1.0383006474  -2.2065988759  -0.9238467962   0.6608411829  -1.3199920866  -2.9084056303  -0.2138984125  -1.8395438611  -1.3730028368  -0.8108237907   0.6493223407   0.3263970858  -1.7186257847  -1.4023685508  -2.2647331719  -1.6764361492  -0.0580929237   0.1464054583  -1.0243103503  -4.0426477890   0.3588390044   0.3392893796  -0.7104878783  -1.4307885882 
		 20 	  0.8177123912  -1.0789085708  -0.2036574254  -1.0902069329  -1.1924394029   0.5175909434   0.9244471208  -2.4569948510  -1.3697055843   0.4592926099  -1.3239087952  -0.5478077082  -0.5977183998  -0.9419927247  -1.0270671673   1.6936324751   0.2022469949  -1.3476680646  -0.0637434600  -1.6448333871  -1.8700042406  -0.5382198849   0.9321532939   0.7078486855  -0.4962386464  -0.2562567818  -0.8247258727  -0.1792440492  -1.6976825825   1.4547113439  -1.0768971685  -0.7682806106  -0.0664955111  -1.0174176347  -1.9452806881  -1.8633193625  -2.2698029166  -2.0213969082  -1.7319492423  -1.0406797896   0.0715261590  -0.0885177747  -1.2267030550  -0.4508897388  -0.8223472392  -0.6627995750  -0.1922921525  -1.4130187945  -1.2283667197  -0.6653069133  -0.2237376797  -0.7708748197  -1.2940726659  -0.3429716735   0.4548451623  -1.6952328088  -1.4072714201  -0.6666337732  -2.1165924603   1.1521493871   0.7853553421  -0.8009559154  -0.3153018218   0.7117048743  -1.6811081785  -0.0937410011   0.2856868206  -1.9625244100  -1.2361753244  -0.2131790558  -1.8605832123  -1.0904389021  -2.6294872730  -0.7443423464  -1.3644879668  -2.2620917753  -1.7254610625  -0.1479346126  -1.2589575885  -1.0377666071  -1.1918261005  -2.8452719455  -0.7430656369  -1.4146395967  -0.2914572393  -1.1627787542   0.0982979504  -0.7630515866  -1.4181090864  -1.4966104435   0.6318947262  -1.2063507650  -3.0970134997   0.1308462816  -0.6743629105  -1.6714528313  -2.2626568456  -0.4023162800   0.0747900898  -0.7050809956 
		 21 	 -0.7903865074  -0.5659983985  -0.6351581955  -0.2076279228  -2.6475465852   1.7830948128  -0.2341621894  -0.1024375377  -0.7268076125  -1.7924475970  -1.5742291698  -0.3428157307  -0.9447139372  -1.1690168284  -0.4444834891   0.3406426043  -1.4879291934  -0.4321797411  -0.6270079423  -2.5003419818  -0.6680917943  -0.3334806774  -0.1958107828  -1.7062650024  -0.5427828804  -0.5170858086  -1.3256924133  -0.6760519470  -0.3721948502  -0.4863948841   0.5433993160  -2.9466188539   1.6230882583  -2.5759520657  -1.7510369364  -0.6441017531  -1.3687481575  -2.7196504147  -0.8463628314   0.5155312277  -0.5247169069  -1.4054557368  -3.0237215550   0.4281477994  -2.9724976331  -0.5891571156  -1.8597877286  -1.3816171511  -0.3461793100   0.1272309289  -0.5961576646  -1.0289004446  -1.5261085049   0.0161885125  -0.7646604554   0.1833479093  -1.4620082823   0.5826036480  -1.4724146384  -0.1381683454  -2.2307461316  -0.2631357107  -0.2002687468   0.7124886283   0.5418711272  -1.7836433298  -1.6797979064   0.7759645646  -0.2941195472   0.0307119879  -1.9869108274  -2.0670331927   0.0229220034  -0.8037353647  -2.1581777075  -1.5061966291  -2.0708489546  -1.5874294730  -0.0755831194  -1.3471294142  -1.2066074784  -0.8443551539   0.2262907493  -0.5061407516  -1.9592032651  -1.5807207919  -2.4474339952  -0.0983976464   0.0610843396  -0.4948002444   0.4622173247  -1.4104540533  -2.3045593345  -0.4273403845  -3.2346671854   0.6404166552  -0.3063737296   0.1219833300  -0.9342978358  -2.6672061274 
		 22 	 -0.5377330139  -1.5443417478  -0.9942520098  -0.4374752754  -0.7745444438  -1.0789317209  -0.2924011586  -1.3499338432  -1.8482788245  -1.5589498983  -2.0801454142  -0.2988356437  -0.7820824015  -1.1497140631   0.7933875178  -2.2576290725  -0.3211210285  -0.7775802811  -1.7451136579  -1.9079908502   0.7512764310   1.3200412839   0.6307456627  -1.2294642109  -2.7048102160  -2.1575159043   0.8854751543  -0.0983976464  -0.1771880567  -1.0780415316  -0.8170612100  -0.5274060923  -0.6647325693   0.0910897775  -1.0063495708  -0.8949373157   0.5138747286  -0.6439660965  -0.8779180286  -1.4337735829  -2.7348765373   0.7695032736  -1.5878875737  -0.5865916243   1.5263075938  -1.9217699418  -0.4260394781   0.0876516456   0.0816654129   0.3520429657  -1.1768795942  -1.1460392104  -2.6326573394  -3.3448677536  -1.3962246653   0.2448851014  -2.7494375102  -1.3359345699  -0.2608732261  -0.4481713640  -0.7872298157  -0.1989141882  -2.1580295910  -0.0689414795   0.3306240871  -0.4663613506   0.2050374902  -1.5684524419  -0.9894469755  -0.6225148686   0.0359467354  -1.9622295711  -1.2204195196  -0.9589053421  -0.3449452068  -0.9088914239   1.5019625428   0.1403721707  -1.5513195477  -1.1988775161  -2.2571895231  -0.2985301427  -1.5933762070  -0.0812140096   0.1952038085   0.2018655229  -0.7512216041  -0.7254782464  -0.1347067157   0.3600536809  -2.6908872452  -1.2724547767   0.3875457910  -0.5970078110  -1.6888241621   0.3489266348  -1.6535377487  -1.8079719505  -0.7271719460  -2.0892434189 
		 23 	 -2.3724767294   0.5596582620  -1.0269532330   0.6571446431   1.7627291732   0.4752758403  -0.8602530145   0.0143926847  -0.0651550313  -1.9062664055  -0.2532206951  -2.2105888103  -1.5059014175   0.7326355222   0.0171657592   2.6983986943  -0.9198908510  -1.6782413931   0.3280497228  -2.4642388285  -0.5268866925  -0.8882224104  -1.5875135077  -2.5152607106   0.4227550142   0.5114737297  -0.4327649500  -2.6722586556  -1.0218963651  -0.5162206078   0.9500648156  -0.7823842355   0.9767367501  -0.9676958130  -1.6262514991  -2.2861502407  -0.7842822150  -1.0531876480  -1.2564110756   0.0184335868  -1.7142118582  -0.9181578754  -0.4358588126  -1.0848320778   0.7732078212  -0.7307874261  -0.9347470752  -1.3389625295  -0.5319502249  -1.9795102578  -1.5418851222   1.1569524009  -4.2952339356   0.0533341007   0.6166179037   1.4190137436  -2.0401885044  -1.7052815961  -1.8183464596  -2.1499618517  -1.9190188893   0.2308127789  -0.3740041656  -0.0754598228  -1.5281401395  -0.3740825128  -0.1874105898  -2.7448888969  -0.5037872159  -1.5093042031  -1.2468973023  -1.9104931017   0.8507551714  -1.2970850325  -0.6224808306  -0.8636022006  -1.0179802125  -0.9287099402   0.6274264318  -3.0522055818   1.3102273082  -0.4549856201  -0.9726025476   0.0105044268   0.4194448721  -2.5271352677  -0.4221183310   0.8171751932  -2.3274265836  -1.5328186972  -0.1688684394  -2.1538185026  -1.0829192867  -2.2020137321  -2.1010661940  -2.0170497350  -0.1790664662   0.6498715893  -0.9369454754  -2.5623358868 
		 24 	 -2.4302410801  -0.5985946404  -1.0174300109  -0.4855422340  -1.2964779950   0.0577419527  -0.1008881885  -1.7334355062  -2.2887840128  -1.2765847077  -0.8163596580   0.2739507132   0.0517779181  -1.1812771375  -3.2282477665  -1.5820732468  -2.2557004966  -1.2581430294  -2.7903864305  -1.1058735001  -1.7075718263  -1.4202660201  -2.9794296484   0.3007106393  -0.8703553152  -1.8137786923  -2.5600242378  -1.7412195348   0.0071244197  -1.6006122919  -0.2419251464  -1.1509866162  -0.8582503245   0.8411451620  -0.7363312729  -1.1856259320  -0.6359595248  -1.0655288039  -1.1947760629  -0.1811128255  -0.4675911636   0.3434574698  -0.0744255383  -3.3636843096  -1.7874179451  -1.0435534377   0.8041060916  -1.6322515225  -0.3900289900  -1.0870381997   0.2888298535  -2.1538716436  -2.1099605376   0.4557211285   0.3619094158  -1.2228008014  -0.3437156497   0.2710023747   1.2677593472  -0.2048393687  -0.2227493046   0.7676764686   1.5514096867  -1.7038207983  -1.4321913522   0.6511768093  -0.9990882264  -0.3667235069  -0.7915102912   1.3563587019  -1.4951922972  -1.4022303164  -1.7621134692  -1.3005467296  -1.4632135428  -1.3033515643  -1.3810708930  -1.9063490924  -0.4288508359  -0.0646918256  -1.0037642554  -0.6007012173  -0.5656513691  -0.9825516895  -1.8953006999  -0.9363758504  -1.7046573707  -0.2588901393  -3.4959904063  -0.1839754729  -0.2237376797  -2.6379640650   0.3572462436  -1.2234610991  -0.9839205587  -2.5062844698  -2.8182307487   0.6430967225  -1.7441158485  -1.2708979618 
		 25 	  0.0353116210  -2.2229380863   0.9327669383   0.2267089294  -2.2194295425  -0.3773821605   1.4521072153  -2.2360597392  -0.2088351069  -0.5033536627  -1.0052989378  -1.4696991825  -0.2862035921  -2.5181917776   0.0809810263  -0.3275940551  -1.4173258158   0.7806340869  -1.0398345280  -0.5077606146  -2.7925422995   0.5261116182  -0.6535793343  -0.2708234811  -0.4264121757  -0.5191178473  -2.2160384798  -2.8581387463  -1.4431338988  -1.3324213199  -0.4838180149  -0.5645431169  -0.5840009633   0.8042251810  -1.8368181961  -0.2139124752   0.3717927478  -1.3098722972  -0.5634356789  -0.0165924173  -1.4166072816  -0.0076037446  -0.4214219263  -0.6879105335  -1.2183954737  -1.1671124240  -0.4730197072  -0.5424197518   1.9323419606  -1.1941149674   0.0723771846  -1.5108659693   0.7543464148  -1.4920918536  -1.0746797306   0.1568361715   0.6400470359  -1.5685772578  -1.8067679403  -1.0578917931  -0.5595677622  -0.5536848788  -2.1226999088  -0.3254766963  -1.0729591411  -0.1584260100  -0.8949907322  -0.1055513981  -2.3196517271  -1.2438048937  -1.1331909812  -0.8389201740  -1.8419737461  -1.9930979007   0.0739907831  -1.7209933679  -1.9015410878  -1.8438470869  -2.0332689545  -1.2373622206  -0.1578480150   0.1958553700   0.1015522684   0.4938022455  -1.8619307329  -0.8708901020  -1.3170069918  -2.2160850808   1.0782453848  -0.3107317437  -0.3940714144  -0.7270606532  -2.6161794303  -1.3926075791   0.4186532438   0.6884437007  -0.9823621269  -1.2159937178  -2.1907636711  -0.6024919319 
		 26 	 -2.2430033159  -0.8421142383  -2.6324426282  -2.8633983739   0.2884371091  -1.0420100775   0.3490980347  -1.0387094305  -1.7341175166  -1.0608484481  -1.1892796245   0.3886604638   0.3169864773  -1.6846792260  -0.8267294250   0.2443363229  -1.4593927671  -0.0693496129  -0.9429878022  -2.5595226112   0.0670864351  -1.0319181961  -1.2839795611   0.9705758350  -1.1651099972  -0.0958567885  -0.8234872337   0.7742716120  -1.0435926415  -1.0394843427  -0.8529734691  -0.7706019292  -1.2657234321  -1.5869247765   0.0893614757  -0.0539058524   1.3134831071  -0.9709594757   0.2769172544  -2.0553892452  -0.0852726380  -0.4063219382   0.4342907166  -0.8224099811   0.3310142220  -0.8422373641  -2.2734815023  -1.2786890171  -0.9312937534  -1.3187931657  -0.5803759571  -0.4446168888  -0.8353266833  -2.1322547156  -0.3081160690  -2.0451150979  -0.9320299954  -1.2827269455  -0.1989702389   0.8500880030   0.4254176874  -2.3018090756  -2.0949546156  -0.9675734251  -0.2671625128   0.0827177823  -0.1150867547   0.3401993754  -1.8951142625  -0.1695358639  -0.8163596580  -1.0988539024   0.1019350866  -1.7421179474  -1.0959334488  -1.1797274192  -0.8428427901  -1.5534950136  -0.7596425045  -1.3357727151  -1.0303047500  -0.1772343062  -2.4585777179  -0.6339776353  -0.6049061725  -1.0387427400  -1.9739853382  -0.8623019529  -0.8182131711  -2.8742854651  -0.4992937723  -0.0188286291  -0.1329962964  -2.2523065928  -1.2326929367   0.4474839180   0.8444815709  -1.0190532980  -0.7860608377   0.8604412389 
		 27 	 -0.1468442083  -0.6225832410  -1.4816737345   0.3473045202  -0.5863601770  -1.3252464837  -1.0128756771  -0.0235779672  -0.6716992248  -1.1391347860  -0.2515660260  -1.4153379049  -1.3418488197   0.1643642056  -1.1153988030  -0.2696101466  -0.5851782125  -1.1571627488  -2.2339823671   0.1351637740  -1.4054585942  -2.5587283169   0.2287354093  -1.7240072485  -1.6234413255  -1.0056199836  -2.1371892566   0.7820011948  -2.1185812498  -0.0813654683  -0.6800509127  -0.4840377714  -1.4825355950   0.4673114842  -1.3072739252   0.1887048401   1.0469844735  -0.0633005238  -1.1512346092   0.9206811692   0.0282496196  -0.3145418441  -2.2478095878   0.1748515305  -1.0720750797  -2.0233445418   0.1272780886  -1.1445801170  -0.2546663266  -0.9992590449  -1.6714600932  -0.4678354257   0.1426399907  -1.3747591281  -2.4472948914   0.1453007359   0.8330466568   0.1061176650  -0.6127436853   1.6034155178  -1.7747425986  -1.5117697909   0.5359047338  -0.5799507929   0.4436778481  -1.2051750260  -0.3873889892  -0.6471890704  -3.0708540612  -0.8532005503  -1.0546967807  -0.3354909856   0.4233984839  -0.6931751680   0.5807279112  -1.0483151627  -1.4693903377  -1.6018940644  -2.2525177541  -1.3678234588  -0.9595261947  -0.1885728626  -0.7248312205  -0.9400421401   0.1843251494  -0.0276423807  -2.1406969936   1.4235338512   0.4058315548  -1.7549853234  -2.2126699721  -1.3992120551  -1.7712904504  -2.0849535484  -2.1819288780   1.6572517707  -0.8881421098  -1.2270323682  -1.0090822555  -1.6371516857 
		 28 	 -3.1831099555  -0.8741203601  -0.8269584093   0.0513886682  -1.4713902029   1.0301467132  -0.8035604687  -0.7717635660   0.5235847256  -1.3808982472   0.8187280068  -2.3302566269  -2.4430370759  -0.9462661216  -0.1747382751  -0.2970571009   0.0114594489  -0.9122167954  -0.9989909008   0.1754751874   0.4568593512  -0.7792890653  -2.1430878729   2.9401528722  -0.6283105735  -2.3164701562  -0.6756478337   0.1206541176  -1.8859239407  -1.0720128340  -0.6824130804  -0.6593531890  -0.9138366298   0.9406185157  -0.4252440280  -2.1227913388  -0.4635717610  -2.6129388620  -0.6653576359   0.5782425449  -0.1444222562   0.0830126669  -0.3927684309  -0.8101808548  -0.0717506834  -1.1231463153  -2.2985962628  -1.0597748003   1.6736307171  -1.3319253262  -0.4158666643  -1.1862636134  -0.8984698908  -0.3678011587  -0.3246741483  -0.4885937628  -2.9947661491  -2.0223467578  -2.8556236803  -0.2530109372  -0.1166949980  -2.8474450203  -0.7154782135  -2.0037326271  -1.9781510332  -1.4795617270   0.7171787466  -1.9800372364  -0.5032810389  -2.2408759681   0.0826179238   0.3664542311  -0.3096213731  -2.2376873783  -0.3573507473   0.4482250439   0.3324980323   0.5520793848  -0.7041202962  -1.9939822009  -0.9150616218  -2.5821082311  -0.7824025825  -0.9191660924  -0.6917715461  -1.5374126147  -1.3026270039  -2.4394792576   0.5378397654  -2.2363779540  -0.6354819745  -1.7119852574  -1.1026523245  -1.8882311425  -1.8438586705   0.9126777447  -1.3193488028  -0.1420039145  -0.9684962801  -2.5741945472 
		 29 	 -2.3822667979  -1.5213877934  -0.2475827005  -0.3485682334   0.4802355826   1.1912763563  -1.9239642539  -1.4447157216  -1.1327672087   0.0093272618   0.8102341138  -0.7272331865  -0.8805510285   0.1568901975  -1.7549890383  -1.4084606027  -1.2279408162   0.4516167718  -1.5842402158  -1.1327876571   1.3923758878  -2.7598395713  -0.5184130036  -2.0726696047   0.8899118901  -2.7258110942  -2.7405224056  -2.1346036129  -0.1972350357   0.8244092246  -1.1822236705  -1.3916108082  -0.6469858971  -3.0151256684  -0.2387930047   0.7270224252   1.9658610841  -0.6478392719  -0.6779154864  -0.6310571171   0.1891472831   0.0235798577   0.2066220358   0.0890703835  -1.4162763963  -1.5890206331  -0.2577987794  -2.0077605680  -1.1828888110  -1.5821083562  -0.8615297776  -0.2375793131  -2.1604875992   0.4588165028  -1.7022603190  -0.2367360407  -1.2858610967  -2.2235184775  -1.0788630858   0.3290125988  -0.1667587908  -1.3275980705  -2.2723884796  -2.1557527152   0.0606563031  -0.6167995157   1.3244594760  -0.2074195766   0.6828305086  -0.0976768829  -0.9174059377   0.0142538134   0.5345231597  -0.7256055912   0.0209845914   0.0918527028  -2.8440487033  -1.2131072248  -0.4844274647  -0.4653432510  -1.7564068837  -1.1562531275  -1.2354147554  -1.6963211855   0.4810075184   1.8799016288  -1.5129316483  -1.1540217477  -2.3155057197  -2.1586721584  -0.9157487304   0.3444843338  -2.2048728400   1.8999512246   0.3357062715  -0.5762595122  -0.4542520364  -2.3252732580  -1.6006977428  -0.5514154349 
		 30 	 -0.5857783720  -2.6370156734  -2.1336469188  -1.5045309138  -0.2934606235  -1.5929884602  -1.4862464464  -2.5334531339   0.2517764613  -1.5681816450  -0.2323123184  -1.1572757608  -0.7280252653  -1.3969409293  -0.1142208503  -0.3138786292  -0.8691055561  -0.3663331223  -0.8918555410  -0.7768671425  -1.8060863367  -0.7469880438  -1.4079781850  -0.5717332584  -0.7609780949   0.3973689947   0.0394885181  -0.2991244556  -2.2521578013  -0.9912361081  -1.6103183118  -1.1250291574  -0.0543411591  -2.4010553969  -2.0708185732  -1.4441635241  -0.1407784795  -1.6708178386  -0.7001850173  -0.0884294804  -2.3190983623  -3.2731526014  -1.1384185035  -0.2207652784  -0.5971197972  -3.5119889413  -1.1582718308  -1.1394752245  -2.2081316996  -0.7653517398  -0.1851984001  -0.4823215439  -0.9429316182   0.5919374211  -0.4310739568  -1.1067101344  -0.4044734266  -1.0383049093  -1.5256170129  -0.6857738838  -3.3887383537  -1.6181479560  -0.7139954932  -1.6655157857  -1.0381040747  -2.1344757374   0.3952421957  -1.7742427931  -0.4256071147   1.0925310748  -0.7557316324   1.6779752704   0.4497589161  -1.5107473378  -1.0778365495  -0.1210276808  -0.0048665050  -2.6713924115  -0.0129880093  -0.8205757847  -0.2157614924  -1.4661422493  -1.0146098391  -1.8274137206  -1.4198546018  -0.9729815293  -0.7668462005   1.3104558828  -0.9408584745   0.5354640617  -0.0743781223   0.7357826926  -1.2655536814  -2.2656182223  -2.2929695726  -1.4014434134  -1.1033656158  -1.1510789744   0.7585439026  -3.0870119446 
		 31 	 -1.6840651084   0.4203905414  -0.7365929086  -0.0814565901  -0.6322574454   0.7798046484   0.7692195785  -1.5319779702  -1.0462033134  -1.4759062420   0.4911863439  -0.1521230520  -2.2036560172  -0.6828580062  -0.3002852162  -1.4632165463  -2.0031207697  -0.1351406828  -1.2535230357  -1.9311416909   1.3446960310   0.9400591215  -0.0259886986  -0.8800114840  -0.8038000543  -0.6760476101  -2.0332152144  -0.5177562861  -0.8915137875  -2.5808232033  -1.0743223556  -1.8230394128   0.0728556067  -1.3314180754  -1.0932717666   1.0201083139  -0.3012893515   0.3661514199   1.9857265977  -0.2109781688  -0.5442210329   0.6784143172   0.8215368720  -2.3532371191   0.0286569780  -0.4447458926  -0.4044676959  -1.8710188394  -2.1746537483   0.2173450989  -1.4976295153  -1.9199677718   0.1805900850  -1.2722168705   1.2536303254  -0.8704344824   0.3819349148  -1.1239741428   0.0028995028  -1.3969823055  -0.9726740420  -1.5339122168  -0.2509421118  -0.9708373013  -0.5464771839  -0.1125090896  -0.4164197386  -1.7586550639  -0.4982760808  -1.1212295344  -1.1927843218   0.3849895365   0.2314507666  -1.4791356594  -1.1164582932  -0.7220461882  -0.0756142408   0.1316503385  -0.3605694569  -0.9165916682  -0.2676354701  -0.4984315021   1.8902509118   0.1528101847  -0.3899482029  -1.9723721598   1.1828020723  -1.3167968767  -0.7659865417  -2.9860341414   0.4810075184   0.3809922496  -0.0903583482  -1.5123417865  -0.9662850759  -2.3058904386  -2.9470412992  -1.6051932273  -1.5891027476   0.9474384630 
		 32 	  2.1106050635  -0.5026746010  -0.5216025724  -0.2864093391  -1.0967094856  -0.4407237338  -0.2007867959  -0.5606753475  -0.8383169427  -0.7648166242  -1.2585826056  -0.4613239338  -0.0457261105  -0.9931800241  -0.7721130481  -1.6663446620   0.2287541833  -1.2792497996  -1.2292108689  -0.1342214543  -1.0682882412   1.3347746975  -1.0239136340  -1.4184447487  -1.1599065785   0.1556533885  -0.2578030130  -1.6426130970  -0.1284045651  -1.5783919941  -0.4467502044   0.1200285528   1.2743408049  -2.2801092917   0.1299184811  -0.6874139329  -0.4895196443   1.0734109107  -1.8377777337  -0.8885841045  -2.6408529707  -0.4623624844  -0.0828514834  -0.5492655069  -1.9120065940  -0.1683462711  -0.5029509125  -0.3867924060  -0.7878103086   0.3972102858  -1.9493630644  -2.2956365683   0.1539860289   1.0203818337  -2.7801324456  -0.9499857352  -0.7087379971  -1.8929002410   1.0900283349   0.0463392163   0.5797436714  -1.2403535362  -0.9856295432  -0.3589742826  -2.9675050747  -2.6474466061  -0.4570084049   0.2755956457   1.4309795592   0.2185513598  -0.1346027073  -0.7582597912   1.6120071778  -0.6335715874  -0.6513374568  -0.8597381308  -0.3163138858  -2.6179308776   0.2362942489   0.4665053571  -0.7475236895  -1.9581897199  -1.7539522071  -0.5848336707  -0.5258115170  -1.0550353686  -0.2319748753  -1.1406084521  -1.2547704019  -0.3716902739  -0.3650280390   0.1313945978  -1.9906428103  -0.0159274054  -0.1273699447  -2.9665651344  -1.7906760719  -0.8526671979  -1.9538614542  -1.4124419201 
		 33 	 -0.4757178682  -0.4998063092   0.3881635590   0.5161381216  -0.4072702048  -1.8639026200  -1.0588871588   1.3573152361   0.2100174058  -1.3456366134  -1.2317789885   0.0855831017  -1.4357763636  -1.2239645234  -0.6086357823  -0.4680329183  -1.0975382194  -0.6738873821  -0.1361619668   0.2219079242  -0.9630195178  -1.3149244879   0.6824844525  -1.4116862875  -0.5751719815  -0.8437573253  -3.4675307787   1.6735058164   0.4137921862  -1.1039667182  -3.9267654363  -1.8221541167  -0.3894556416   0.4274400239  -0.5485335936  -0.2421328256  -1.1860056441  -2.1533270931  -2.1008489540  -1.6221071609  -1.3844077514   0.2686195058   0.7159840887  -0.1292121122  -2.3275798476  -2.4360690029   1.8625045596  -1.1816565770   1.5203852115  -0.1936356531  -1.4436153288  -1.6153272274   0.4479633037  -1.1041084440   0.2957986519  -1.2400247389  -0.7433142328   0.4743325016  -0.1122983911  -0.1783966801  -0.4428508263  -1.8780896298  -2.1113958834  -0.0554863240  -1.8263907007  -0.5768867678   0.4029155948  -1.9493296413  -0.9364145039  -1.3156750207  -0.3134972837  -0.2057227621  -1.0992447928  -0.5893360099  -0.9600081330  -0.6738162482  -1.5618029148  -1.0221333548  -0.7472054054  -0.9032012402  -0.8474925674   0.1911283063  -0.1303218623  -0.4417675875  -1.0864704411  -1.9431494270  -0.0316059867   0.1898577507  -1.2862127904  -1.7583392247  -0.6096988783  -0.2503487361  -0.5802156981  -1.9772063749   0.9899931881  -0.7437422477   0.4105179598  -0.7817125126  -1.6814647000  -0.6907473067 
		 34 	  0.0472512225  -1.0321325276  -0.0974793546   0.2314480919   0.4180372814   0.0264062721  -1.2351065567  -1.4842404202  -0.5136030726  -1.3188806535   1.6615184575  -1.4470373920  -0.2056526340   0.5266660612   0.3484291655   0.0811629742  -2.5889659092  -0.8021194355  -1.3477630639  -0.4646673328   0.3244291138  -0.7198302982  -1.7902285936   0.8789594502  -0.8269362645   0.6376192974   0.5565014323  -1.1854906688   0.6786842943  -0.5742262076   0.0581590966   0.0687930389   0.2067725364  -0.1813770728  -0.3675445344  -1.7857621715  -0.5224722753   1.1256699511  -1.9003024631  -1.8804741376   1.3911975019  -0.5323817427  -0.4476873369  -1.9630812071   0.7687226428  -0.2617174302  -0.2870657214  -0.3540870482   0.4881720314  -1.8365034794   0.0852916877  -1.8822070424   0.7167921444   0.7769209545   0.0048844779  -3.1239009444  -1.4400030481  -1.0726570556  -0.0218025593   0.9373576995  -1.2148083431  -0.4962386464  -1.3971795454  -1.0807167533   0.3644215105   0.3707572641  -0.7923151280  -0.0926554145  -1.4718350757  -0.1203449084  -2.4823563266   0.0505643641  -2.3607700204   3.2041695575  -0.7930871902   0.2463777845  -0.5476797319  -0.8759462794  -0.3415901098  -2.0414519371  -1.3616667487  -2.3366717182   1.2268065294  -2.0054909606  -2.2005771441  -2.4643789892  -0.1895800541   0.5414698522   0.3890640380  -1.5117277538  -2.3514236391  -1.1720163248  -0.1689453306  -0.5337029942   0.2114164747  -0.1276667715  -1.6504263137  -1.1131469460  -2.0693363905   0.0577946879 
		 35 	  1.1056004852  -1.2194547892  -0.4943800863  -1.0836483061  -1.2543043824  -1.6210588204  -1.7416122202  -2.2297180688  -1.3493431555  -1.2613940058  -1.3740472783  -2.2888758194  -0.8195967938  -1.1878945091  -0.2297013426  -1.3191516864  -1.2634018187  -0.4187236574  -1.4293365953   0.3820240493   0.8772022454  -2.8896246313  -0.9356181275  -0.7110284927  -0.7495457032  -1.3991442417  -0.9666053201   0.3385227119  -1.3573177784  -2.5027189286  -1.1624402122  -0.4587005303  -0.2584932267  -0.2015925373  -2.5306312775   0.1982290253   1.3667820209  -1.0893531739   0.4189928604  -2.0330298182  -0.7823063091   1.1660404154   0.1002120744  -1.8912993288  -1.5024610719  -0.3817556280   0.3748828297   0.8791160426  -1.0995864735  -0.5150449276  -0.1669341689  -0.8875141074  -0.8967998978  -1.0827811877  -0.5159365385  -1.3602017649   0.1737503850   0.3782200592  -0.1580148208  -0.1885797812  -1.3118988237  -1.9624305579  -0.5064820666  -0.3306264118  -2.8811693600  -0.6182440673  -1.1733034884  -1.0086387001  -0.3122482235  -1.0784586543  -1.6159218286  -2.2946579056  -0.7891364287   0.4015998583  -0.7584323788  -0.5808604591  -0.8750054760  -2.6686121631  -0.6210713102  -1.6890530446  -0.0349376778  -1.3427085403  -1.2427382544  -1.7455474676  -1.5361873199   0.5841132274  -1.4517958284   0.9326895320  -1.6816789047  -0.7399411603  -1.6731636265  -1.7761739876   0.7289213038   0.1564354252   0.0581303507   0.2894065291  -0.8339477774  -0.2971073494   0.7690294020  -0.1936205565 
		 36 	  0.9626969109  -0.0736637710  -1.5487960345  -1.0555005840  -1.4719974806  -0.5591937123  -0.8589409731  -0.6540781095  -1.8255661885  -2.9571311813  -2.0036438977  -0.8297691770  -0.8665526169   0.8986536924  -0.9174059377   0.0574939137  -0.3886829148  -2.5409173711   0.5586304065  -2.7632476205  -0.5731755619  -0.1276487814  -0.4156593301  -2.5033017300  -2.6286622818   0.5666420182   0.6754850883   0.0366170978  -1.2098579877  -0.2820631623  -0.6354113153  -0.3409968440  -2.3449744848  -0.6732159263  -0.7372641253  -1.5718132435  -0.9321577228  -2.6348662913  -1.5551905988  -1.7553972948  -1.0337194397  -1.4582814793  -1.3156948039  -2.0867192658  -0.2055569507  -1.0650857507   0.2208671154  -0.5771984729  -1.2774342092  -1.5475332921  -1.2037082272   0.0870912160   0.3611750468   0.9135108904  -1.0065756243  -1.8984067402  -0.8053816371   0.8910433230  -1.7854296602  -2.6835499181  -0.7316653181   0.1634303474  -0.4544692567  -3.2713817967  -0.2715942896  -1.6339774426   0.6911912502   0.3599152838  -1.7509172208  -1.1886726763  -1.3075223147  -0.4246942727  -0.4491274315   1.3179392532  -1.0236704902  -0.1440191214  -3.0156402699  -1.6867426379  -0.2785735376  -0.2814612894  -0.1672431116  -2.2072428859   1.8602318581  -1.0564025136  -0.7545110565  -0.6576001714  -0.1044970887  -0.3560739478  -0.9141219037  -2.2128573500   0.0661707235   0.2922450702   0.2749744091  -0.8709413443  -1.1750784106  -1.5307569138   1.8199002111  -2.7774818219  -1.9828020873  -0.1734586755 
		 37 	 -0.3373172371   0.2560973258  -0.6047973583   1.5285961873  -0.8562803390  -1.1101521838   0.2195512011  -1.1905760500   0.6712722330   0.7677061241  -0.4140616662  -1.6018086135  -2.8541082797  -0.1812398248  -0.3349471564  -1.8768480291  -0.0466249268  -1.6572484141  -0.4957170205   1.8863032658  -0.9272127944  -0.3769019788  -0.3458056842   0.5560720988  -1.4563974951  -2.1822092123   0.6636992541   0.6017021748  -1.1874828899  -0.7544739204  -1.1380179058  -0.8416615042  -2.5601834846  -2.9240700535  -2.0615688992   0.7889698022  -2.4170033313  -3.8159611167  -0.5590751896  -1.8692625530  -1.9225773469  -0.7938542991   0.6417492467  -2.8847270253   0.1284111312  -0.7729180861   1.0290253251  -2.2533665943  -0.2831985925  -1.4584229324  -0.6707791007  -0.4083653371  -0.6229583899   0.3886559872   0.4696616019   0.0254976995   0.6055871419  -2.7628936692  -0.5891526558  -2.2216838820  -0.0954943325  -2.5683891430  -1.8594109775  -1.3020662565  -1.7804317236  -1.1445995745  -0.2394960658   0.5268060172  -0.1888042400  -0.5829097000   0.9643616469  -1.1999640107  -1.2371355120  -0.8216018390   0.4888441041  -1.4348584292  -0.3444445420   0.3259525110  -2.2634593626   1.2398973496  -1.9347088752  -0.8665658419   1.2397250322  -0.2841842628  -1.3329034722  -1.5086394151  -2.1351236515  -1.8681743596   0.4575705600  -1.0753453495  -0.1183441352   0.0710635672  -0.6401869092   0.2185156839  -0.4249659034  -0.6090916579   0.0513816615  -1.2706532215  -0.2789784162  -0.2044365893 
		 38 	 -0.2540130463  -0.6060347064  -0.6585718431  -2.2655815207   0.7011141873   0.7838184302  -0.1467924330  -1.9479108496  -1.8616942136  -1.1786322368  -0.3871232225  -2.0618021271  -0.3358953831  -0.7933773211  -1.3477507450  -1.3761966573  -1.6911600361  -2.1321831617  -0.7826104760  -0.3929363978  -0.9166284504  -0.8818044740  -2.7983212542  -0.7886966175  -1.6245258223  -0.7833643876  -0.8453105965  -1.4220545108   0.4972415764  -0.5620911679  -0.5102327738  -1.3564008320   0.2306504329   0.1650145428  -0.9010423433  -0.5505328019  -0.7033017847   0.0177260553  -1.2760347263  -0.6190560032  -0.8612125912  -2.1749538762  -1.2167770981  -0.3266369836   0.3331260094  -0.9823511871  -1.2856493526  -0.8912124071  -1.7595489298  -1.0976870434   0.0315928791  -2.0474018576  -2.0203891104  -1.2493472096  -1.3993260586   0.0877049116  -0.4582659011  -0.0559246301  -1.4844501541  -0.0646748899   0.1651237476  -1.3729216447  -0.4169763585  -2.0844170134   0.5712992692  -1.5004322218   0.0568798970  -2.3838249743  -0.2063420413  -0.2808633387   0.0845991574  -0.8535683478   0.2405619548  -0.1682775525   0.6099974717  -2.7918107069  -3.1906609614  -1.5599006843  -2.0085628267  -1.4150607108  -1.6638892322  -0.3833060976  -0.7672850784  -1.3531441267  -2.1385231117   0.0649637285  -0.0699275082  -0.7427487640  -0.6332771677  -1.3853414558   0.6878759385  -2.1094049197   0.8564672761  -0.2684224969  -0.5731646807  -0.9160683599  -0.1662621797  -1.3127635477   0.7549539194  -0.0131640587 
		 39 	 -1.2844487540  -0.7188884198   0.0770767115  -1.9099835502  -0.1263135185   0.4274400239  -0.0941726197  -1.4502867126  -0.8182237085  -0.0820774464  -0.6243224421  -1.0379007497  -1.3750139394   0.5592471120  -1.1143563403  -2.8458398060  -1.1524120331  -0.8524970789  -0.9056552053  -3.5155788497  -2.1794608480  -1.0889382317  -0.2723459990  -3.2448156322  -1.0742684783   0.8653030589  -1.5288733177  -1.2431982508  -2.0315355153   0.0718235225  -3.1549137695  -0.7983473021  -0.4036076351   0.1525073598   0.5424534992   0.8614435897  -1.3361951132  -0.4616272382   0.7990938006  -1.5991090236  -0.8760713648   0.0784656526  -1.9039451115  -0.5843967795  -1.4246657882   0.5659247510   0.3621555988  -2.9385430120   0.9551665719  -1.0786513561   0.7556483666   0.0537660314  -0.2048885233  -0.5940145084  -1.0816522559   0.2210450959  -0.9181222519  -0.9594289738  -0.4173537978  -1.8518458322  -0.9960508107  -1.1390844689  -0.4092137314  -0.1663988771  -1.4218533748   0.0224277134  -1.3158234427  -0.0573946221  -2.9134207391  -1.2644641506  -1.3787810683  -0.8374967420  -2.1527468503  -0.7562724311  -0.8494900676  -2.2996566659   0.0947048156  -0.4370270670  -1.6975166369  -1.0084246006  -0.9819658581  -1.0994897676  -0.8979990004   0.4557023961  -1.5230966967  -0.7284689361  -1.4689215435  -2.0361026431  -0.8320128595  -1.4312957053  -0.1036506197  -0.5590969605  -1.5560901268   1.3334501338  -0.4229240886   0.4409831460  -0.2591076130  -3.4111715703   1.1346379647  -0.0357523693 
		 40 	  1.1703428451  -1.0748142023  -1.0526445375  -0.5429256339  -0.6702770597  -0.4168981659  -2.2290920756   0.2531511196  -0.7610672721  -1.7563267092   0.7254690340  -0.4639810081  -1.3700899937  -1.1137354687  -1.0059606464   0.8133711563  -0.0789821384   0.1662416161  -1.8168222785  -1.7607157941  -0.0279820764   0.5546345066  -0.2682342088   1.4436400873  -0.3745767005   0.0532044569  -1.8783351158  -1.1092694530   0.4035411364   1.3818284035   0.0827556604   0.0634344726   0.3407853413  -0.6283533029  -1.8578142407  -0.1227820041  -0.9897882824  -0.3070112492   0.3427655309  -1.0131532815  -1.5968547897  -0.2240709701  -1.3712522548  -0.5080177441  -2.4465127181   0.4532273819   0.9827594547  -1.1635556184  -1.0068632468   1.2783408444  -2.4305370729  -2.6919231292  -1.0541325432  -1.0957001102   0.5429232362  -1.9208657194  -1.0371234358  -0.7212061715  -1.9267821724   1.0931917403  -1.7339177690   0.1911822642  -1.0479388503   0.9366448626  -0.8040515205  -1.4754685466   0.0353823059  -1.4352195204  -1.7355422505  -1.5695080741  -0.9031353617  -1.3906803460  -1.6558041284   1.0105984701  -1.0980748055  -3.5724308121   0.7389249376   1.0955844352   0.4470444485   0.1979071137  -1.3319504454  -0.1723985108   0.8406775131  -1.3707965613  -2.4552257455   0.8483954968   0.3667616632  -1.5508277756  -0.8563345058  -0.3433722249  -0.0575667964  -0.0665838340   0.1592724671  -1.8758649411  -1.6437364942  -2.5196468403  -0.3348400524  -1.5000626791  -2.6459653066  -1.3805352099 
		 41 	 -0.7444520915  -1.8389457091  -0.2094819818   0.4252917174  -0.6930969743  -2.2470229134  -0.6441061499  -1.6169700356  -2.0453169058  -2.3171316778  -1.4403762799  -0.3047350215  -2.3907055520  -1.8068765886  -2.4807993938   0.7715096908  -1.6378106601  -0.2902855525  -0.8177142091  -3.3094706952   0.0924572502  -0.3588229435   0.7373793253  -1.6948994370  -0.9266298406  -1.4277430581  -1.1981855364  -0.0681978947  -1.0435554261   0.1844024617  -1.6731158649  -0.9485496489  -2.6848705106  -0.1552455691  -1.8091379724   0.7135344313  -0.5335170205  -2.4959125429   0.9150349652  -1.6696044992   1.2607844850  -1.6536389859  -0.4759231468  -0.4737431679   0.6156380839   0.1935635123  -2.5111599061  -0.4404836967  -0.3910465029  -0.1016218577  -1.7143941571  -1.1130235814   0.4030306520   0.4995269769  -2.0358434188  -0.3825280161   0.6276796402   1.2084148351  -1.9793547494  -1.2638050604  -0.6175858711  -2.4654607692  -0.6128328379  -0.7787278954  -0.0229892063   0.4100473751  -0.8022321094  -0.0169070877   0.1679188891  -0.4280151400  -0.0206239188  -3.2477863188   0.3788118408  -1.3261483029   0.4757541730  -0.9641672814  -0.9674886786   1.1481501249  -0.5791969051  -2.1566287929  -1.4736103697  -1.1695801757  -0.7523536651  -0.9522695399  -0.5492713307   1.3371193435  -1.3647848246  -1.9278473564  -0.6204703307   0.8813284409   0.8217030530   0.3902048351  -0.4134692445  -0.1099433259  -0.9694605952   0.8913870120  -2.5201600194  -1.9516583447  -0.0006962103  -0.9645520330 
		 42 	 -1.0960855200  -0.3440851603  -2.5430766385  -1.4826824974  -2.2679484518  -2.2231342757  -0.0219538840  -1.7589900186  -0.7537562585   0.0995581630  -0.1267982198  -1.7416915324  -1.9331775650   0.7228392769  -0.4816108200   0.5712713703  -1.0210574363  -1.0722393838  -0.6977870516  -0.6985946728  -0.0533174322  -1.5023579234  -0.7148361996  -0.6475571636  -0.7983982219  -1.1535817663  -0.5547640753  -0.3835511042  -0.8976663719  -0.0895343316  -1.1842920782  -0.5004869902  -0.9099642294  -1.0270946655  -3.1346861471   0.0835236946  -2.1960773487  -2.5327037236   0.3000582710   1.9723847498  -1.5838589232  -0.8359416496  -1.6834931618  -1.9494861050  -1.0642079091  -0.1824072479  -0.4433067702   0.8115923117   1.4318077108  -0.8153803684  -0.1326761469  -1.6814493165  -2.5799001172  -2.4385861884  -2.6936877186  -0.3175508281  -1.8546373987  -2.4072221173   2.0069696061  -1.3196587123  -0.8098655089  -1.5419367244  -0.5473203072  -2.8932557071  -1.0499104125  -1.3429118284  -0.5258208713   1.7063241564   1.2674858313  -0.8199411505  -0.7270701205  -0.9773515630   0.4684666516   1.0441862327  -1.8191417756  -1.5286718803  -0.3295014250  -1.1657806974  -1.7970984456  -2.5492015365  -2.3489388072  -1.2043330358  -2.3898272304  -1.7022909140  -0.3363974599  -2.0982747909  -0.6598270985  -1.5720844264  -0.5348199373  -1.1224734204  -1.1178587255   0.9606523670  -0.1309183070  -2.0356767878   0.4139904349   0.0205447933  -0.3706741334  -2.5488669616  -0.5055414581  -1.6935862339 
		 43 	 -1.7263727359  -1.0557369106  -0.1707958349  -1.9204202387  -1.8521745603   0.1641538468  -0.7448739366  -1.0659615845  -2.1330745277  -1.6073003633  -0.9935871406  -0.0230175543  -0.4221430597  -0.9285229127  -0.4221430597  -1.3399004873  -0.1764102001  -1.2527835988  -2.8676522720  -0.5790193554  -1.7474116845  -0.9658077329  -1.0975062940  -1.4037135161  -1.4644925230  -0.0903885938  -0.0202835889   0.6634560449  -1.2584515739  -1.6664896764  -0.1274673365  -1.4701880168  -0.3589949413  -1.0825767855   1.2960730174  -2.2235708825  -0.4231350650  -1.1017719746  -2.6737853999  -0.9388542438   0.2827106995  -0.2803924622   0.1321192119   1.0627491404  -2.4400577842  -0.9609143342  -0.7777395425  -1.9711612003   0.1317741901  -1.3425685292  -0.5586205725  -0.4543496618  -0.0470823071   1.8021000843  -3.2612967680   1.3301354294  -0.8157169631  -3.4965256811   0.2105924577  -1.8677340117  -0.0154746205  -2.1904982540  -2.2499804995   0.2303974946  -0.7534595794  -0.7288327420   0.1549461330  -1.9291526109  -0.6582684550  -0.9329122692  -0.7469935531  -1.1556621896  -1.8088919650  -1.6161221746  -0.4994883376  -0.4771377018   0.1572521196   0.7225565236  -0.2476434295   0.0475199529  -1.4928409600  -1.6206656821  -2.5243491115  -1.1060922447  -0.9918831794  -0.7743810573  -1.3318198718  -3.1013190760   0.2253936476  -2.7150459287  -0.0032131750  -1.4255358536  -1.6931949086  -0.1368380004  -1.3303903569   0.0476154954  -1.1388213610  -1.6361828381  -0.8488324671  -0.9358142055 
		 44 	 -0.5413292139  -1.8697570821  -0.7864332605  -1.9992885361  -0.5177559353  -0.0558567740  -3.2565540847  -0.0270656026  -1.5824419426  -1.6609305444  -0.3421047032  -2.2866550255  -1.4420415282   0.6923198511  -1.6320010777   0.8002105227  -2.6056864073  -1.0465079756  -1.8154769201  -2.6607221185   0.0296124707  -0.0763178558  -0.4363521991  -0.4500591286  -1.4112370684  -0.3947302567  -0.6100044345  -0.9405218738  -0.5464294159  -1.2364199313  -0.8324016695   0.4826438344  -1.9870808456  -0.9622352652  -0.8989423193   0.2865416921  -0.0192108154  -1.3497947607   0.3648856542  -2.0712080705  -0.2222080155  -1.8724708181   0.1925006071   0.8761746996  -1.4068081114  -1.6302662732  -1.9333763664  -1.6304959492   1.5301487470  -1.7223341312  -0.4648471454   0.1425836194  -0.1444729453  -0.9386097287   0.0734516537  -1.1521271565   0.8057692064  -2.1785316025  -0.2150428352   1.0094468865  -0.0070708424  -1.2106053346  -1.6634218360  -1.5857187603   0.4824733671  -0.8805979472  -1.4532986297  -1.1783092534  -1.9151126363  -0.2027477855  -0.6919786093  -0.7230536408   0.5440295469  -0.4801733456  -0.4206475610  -0.3710486352  -1.7030295052  -1.0969939578  -0.3596510772  -1.4705017256   0.2149464304  -0.4533658986  -1.2214920447  -1.0525690072   0.3232062910  -1.9300680300   1.1655535024  -0.2855632066  -0.1217593183  -0.6236472321  -0.4028626577   0.2938220069   1.5527918085   0.6926675238  -2.3279936768  -1.6390814862  -0.6728321108  -1.0154107640  -1.3824490406  -0.3862126050 
		 45 	  0.2539235295  -2.2668736427  -0.6412137104   0.4913479286  -0.4801871153   0.3387981400  -2.5141016981   0.8949765521  -1.9582248964   0.2922705181  -0.9530182919  -1.7038284759  -0.5276354441  -0.6432871696  -1.6381854467  -2.4479234789  -0.2144499498  -0.3815534561  -1.8027300484   0.0228559760   0.6740582896   0.0563860142  -0.0200710340  -0.1358293629  -0.7525524176   0.1725287894  -1.2517817950  -0.4632936998  -0.6866589719  -1.2022379519  -2.1148577417  -1.5643526876  -1.5951826555  -3.9675904007   0.0790493563  -2.2016540290  -0.1096480315  -0.6083659926   0.3331171722  -0.2598620658  -1.1064691639   0.1725287894   1.0469844735  -1.2586941279  -0.9420641800   0.9193750930  -0.5990257548   0.5865020872   0.8125010180  -0.7872624573  -1.4284625001   0.5918202995  -0.5776802040   0.5982149968  -2.4507838854  -1.7037116671  -2.0855312963  -1.6165363593  -0.7387548814  -1.7989706819  -1.4500488879   0.5978359370  -2.0572795808  -1.1329798623  -0.9247719984  -1.0465588612  -1.3790358078  -1.3296887454  -0.2001365405  -0.7128696515  -2.8650669258   0.0137024297   0.7470546002   0.7425340168  -0.1045634623   0.4884993685   0.9688346749  -1.6508266698  -0.9022915271   0.9365231885  -0.7324946067   0.0867781130  -1.0046811615  -0.5333986135  -0.1817225300  -0.4731289015  -0.1402812476  -0.7519079384  -0.0790396133  -0.8561784890  -0.7734805904   0.2803675799  -0.9617035129  -0.8390735754  -1.2195574683  -1.5524729182  -1.0412975584  -1.3279648156  -0.0011849584  -3.7396270008 
		 46 	 -0.5134104529  -0.4870471912  -1.5990575940   0.6644770471  -1.5471271433   0.8275048580   1.2070894963  -2.0548604809  -1.6320868625  -1.1676329322  -0.1111469009   0.4499355406  -1.0357020075  -0.1957136480  -1.0891329210  -1.8009363992  -0.2100194078  -1.6479114442   0.0229139928  -0.5021215684  -1.0058446941  -0.0402753102  -2.5601767020  -0.0604530679  -1.2906541706  -2.9195260879   0.7337655981   3.0839777786  -2.1595859007  -1.1777653659   0.8159762667  -2.5499765230  -0.5221388858  -1.0083771295  -1.1662296353  -2.3275354248   0.8902827811   0.2910667846  -2.1296014636  -0.6172825127  -0.3280571959  -0.1494044156  -0.4008078360  -1.6381210855   0.0595722351  -0.3221810241  -0.1465069206  -3.2376985908  -2.5114094481  -2.6739559432   0.8041475353   1.0236000252  -1.0248004928  -0.2153045247  -1.7475361921  -0.5970089919   0.9923151945  -0.6396860626  -0.1342678721  -0.1016185911  -0.2702571757  -0.3365722117   1.1080658273  -1.8899070633   0.0812065855  -0.2395767727  -1.2155235842  -0.4946086363  -0.7962597351  -1.8430706308   1.7740763836  -0.8055015270  -0.6486435381  -2.5774538737  -0.3418419143   0.6610547953  -0.6154198943  -2.3926167424  -1.8998983471  -1.3074524312   0.4942595838  -1.6533355917  -1.3927343787   0.7209168750  -0.4771377018  -3.2254089204   1.4510330657   0.3140351375  -1.2681958932   0.1563593731  -2.6705719276  -0.0983930861  -1.2523057797  -0.7771139826   0.6242666120  -1.4258225286  -0.7399707796   0.0389786144   0.1291165804  -2.3919752896 
		 47 	 -0.0027741038  -0.5326487198  -0.2817207843   0.0953233630  -2.1885900354  -0.5443750591   0.0119650063   1.4733853976  -0.1774291140  -1.2090478462   1.1589899199  -0.8715361174  -0.7454904884   0.3344042995  -2.9824527456  -1.7177786627  -1.5880781874   1.1848907563   0.9371952007  -0.1731805911  -1.3328516431  -0.1748402449  -0.5257391739  -0.5447476631  -2.1795180454  -0.2972939308  -0.9137560434  -0.7961617005  -0.8904273059  -0.9073978868  -0.3011085797  -1.5307804736   0.3986406467  -1.0719953772  -0.7747504803   0.4701971841  -0.3744258386   0.7293435983  -1.7653362907   1.1757435613  -1.1883156831  -1.1423952619   0.0018254144  -1.5261432954   0.2937965964  -1.6908452473  -0.4408466969  -0.5306831115  -2.2797586745  -0.4072828097   0.4359346479  -1.2812338510  -1.7683119687  -0.6134231728   0.5892200273  -0.0334860543  -2.6581942334  -2.1269984397   0.3755206615   0.4611339288  -0.4940134721   0.4496812199   0.0681691397  -1.0420164618  -0.5681662887  -1.3061495437   0.2315309469  -1.4379213844  -2.4863387657  -0.3107515601  -0.0963077364  -0.7966658779  -0.0707954840  -1.5204866736   0.1455464410  -0.6797689900   0.1745681825  -1.4773627660  -0.9893451436   0.5410221964  -1.8872913816  -0.3739182575  -2.3963727131  -3.2167703261  -0.1006527071   0.3639455740  -0.5958481076  -1.5644805820  -2.2319599270  -0.3730659972   0.4173915119  -2.1890757492  -2.4994347708  -0.7840366109  -0.8090117183  -0.4561076061  -1.8588173800  -0.5765156668   0.6616774130   1.5408879101 
		 48 	 -0.8300319803  -2.1849850735  -1.4277862539  -1.2218921486  -0.4044137664  -0.6935719187  -0.0524016456  -0.9047144847  -0.9032432579  -0.2702021227  -1.7362589572  -1.7460599261  -1.1985058516  -0.0706217440  -0.4557628191  -1.3552978866  -0.1247281276  -1.1327907113  -1.4783096550  -0.7208031126  -2.3461064759  -0.1799617502   0.2261825121  -0.3904885032   0.6578975281  -1.4592189626   0.6265403656   0.5278825881  -0.9032681210  -0.1676384960  -0.7365238705  -1.0864158619  -1.5306423441  -0.8353263431   0.4814099383   0.7675678022  -0.5349977732  -1.6442699667  -0.9639602598  -0.0950348545  -0.1792297405   1.0768984938   0.7374988661   1.3256008906  -0.5584755933  -0.5079550448  -0.4259850638   0.2818013111   0.8961895368  -0.0019484324  -1.6277497894  -0.1932023888  -0.0077409383  -1.4150910480  -0.0490970485  -0.1415701757  -2.5767629817  -1.7851672179  -0.2976357421  -1.2004724367  -1.1516454197  -1.7604796238   1.8420993471  -0.8767764279  -1.1596406594  -1.7585176736  -1.5989590477   1.0559083479  -0.1299766679  -2.6006039462  -3.2032587895  -0.9218111586  -2.0153353544   1.8337071254  -0.5776338699  -2.8408169224  -0.1144818962  -0.3650060005  -0.9426762570   0.8124465980   2.6935475162   0.0067126926  -0.3896416858   0.9971895212  -0.6737344660  -0.8620960667  -0.8204105975  -0.0147143118  -1.1016058405  -2.1466729353  -0.7018341726  -1.1137523934  -1.6281761298  -0.9830623233  -0.0200265480  -2.5008463558   0.0129400698   0.1479488112  -0.4579641918   0.0488560370 
		 49 	 -2.2575282578   1.9104366330  -0.1102982301  -1.5711062640  -0.2781348649  -0.6407126029  -1.1704591543   0.5435781474  -0.9482183105   0.5540138648  -1.2128554030  -0.5482536222  -1.0872969914  -1.4649654306  -0.1037453545   0.4397634653  -0.8223611452   0.0370589881   0.2739507132  -1.0920314949  -0.9711352856  -0.2566522296  -1.5119123117  -2.4510455205   0.3239969782  -0.3860500808  -2.3015102094  -1.2770016110  -1.4827233698  -1.9758698219  -0.8191035553  -0.7839350720  -0.8439731587   0.3166705785  -2.0427798243  -1.2586811747  -0.8514780864  -1.7802332986  -0.9917166647  -1.0655288039  -1.1445682675  -0.9289053533   1.9202491411   1.5166629175  -1.0418522352  -1.6405724099  -1.1132263522  -0.9646308326   1.5644228314  -1.6578083649  -0.9758325364  -0.9756577911  -1.8790924951  -1.3040183058   0.6441342811  -0.8827601049  -0.2559502881  -1.8222704575  -0.7244580339  -1.0224589402  -1.1660781518  -0.5283814926  -0.1360804236   0.3573868979  -1.2086994760  -1.3129484116  -1.7385405252   0.8753297735  -1.9601573348  -0.5870257601  -1.5503105764  -1.4581199708  -1.5627820801  -0.4222078791   1.0828464785  -0.5331242220  -1.4457139475  -1.5348858073   0.5258219909  -1.5331627478  -1.7610108510  -2.9692387853  -1.3333941252  -1.7629369738  -2.8157711468   0.4687274036  -0.5964804576   0.2390645141  -2.5519589143   0.8188944414  -1.2652119006  -2.2739909455  -1.7421147986   0.1661744733  -0.0431919132  -0.0604715529   0.2633782939  -0.5995065380  -1.8556594829   1.4756576182 
		 50 	 -2.3018783241  -1.1000427526  -1.1169177268   0.1455464410  -0.3721248234   1.1419190614  -1.7679754320  -0.4255906513  -1.9456942763  -1.9866537368  -2.7834946466  -0.3471081459  -0.5677190286   0.0212696138  -2.7463499518  -1.2513814124  -0.4532300884  -0.3246166757   1.3476364243  -0.3614339472  -0.2356203700   0.2915270146  -1.9419760370   0.0914771225  -0.3170056149  -0.9242604473  -1.3898078036   1.2247722064   0.0648570499  -0.4871310581   0.0623691989  -1.9535338768  -0.2411469792  -1.2436487455   1.1348272485   0.3620019662  -2.0714320774  -0.0474649372   0.2382242298  -1.1843615048   1.0335669749   1.5693146229  -2.6723638998  -0.3009574094  -2.2097938837  -1.4921532662  -0.8129955026  -0.9473469532  -3.0255947182  -0.6025166510  -0.4571342841  -2.9280688163  -0.9843852921  -0.6821000922  -0.5158864632   0.2402187368  -1.1715382221  -1.0302422746  -2.1399366760  -0.2013090525  -0.7806231504  -1.7716025914  -0.7120999609  -0.6066151486   0.6812771170   0.4490945460  -1.7792621814   0.4172797229  -2.3568185176  -1.6362910499  -0.9488680973  -0.1478455930  -1.1551165612  -0.1845360256   0.7364726673  -0.0809503537  -1.0429952944  -1.3145655690  -0.6315455542   0.0274287801  -1.9025695707  -1.2719140591  -2.0293610653  -1.3436720073   0.2198211598   0.1880159837  -0.5271047403  -0.8350259135  -0.9165158292  -0.9420641800  -0.9459876826  -1.5140284710  -0.3580654867  -0.1068694163  -1.1917819061   1.5699602273  -1.0141124609   0.7451260591  -2.4141627122   1.5911945887 ];
	Beta_AT=Beta_AT[:,2:end];
	 Beta_FEE=[	 1 	-15.1319150318 -65.7805854953 -29.7476702598 -61.4246922144 -55.7034606148 -28.4545234512 -52.7627981800  -9.8509889766 -26.7510112183 -47.5963891951 -32.7200656103  -8.9314835444 -27.5450622165 -18.1111231099 -54.2907466215 -49.9201272949 -36.3068128271 -33.7490645512 -25.6416405585 -40.4229874915 -31.8524018856 -23.5577044503 -14.6729436159 -30.9071354936 -17.4197973390 -33.6237930149 -76.9702802497 -23.3579600715 -13.9461169849  -8.8455588574 -35.5951369168 -26.6047714294 -18.1649186625 -33.6255606486 -13.7953684340 -62.2552892689 -52.1580378423 -30.3676045683 -39.4158591556 -24.2673378518 -16.6082839538 -27.3171480735 -11.7363812591 -49.6065500379  -3.6121300484 -54.8464596804 -16.9254680950 -29.6184154337 -33.0572496389 -35.6465094247 -32.8701946873 -26.8027622265 -14.0704785186 -35.5652709073 -11.0158070092 -17.4593136189 -51.6382595229 -39.1590662893 -18.9512520378 -28.7036260169  -0.3010969441 -23.8256390588 -13.0039966410  -3.8071573508   0.4320537416 -32.2920369211 -51.0451538626 -34.3263020196 -23.1865149184 -21.4401387426 -17.8017474144 -30.1866756227 -39.7471855498 -46.6996138849 -37.3323243241 -41.8845918406 -24.6299897533 -22.2058462967 -19.8048329274 -34.3850402112 -22.2218606187 -16.2421322093 -33.5787609641 -10.1688606805 -26.9306428324 -52.3775709058 -22.2097245112 -15.5104735233 -19.9843702119 -29.9207812068 -28.9741765854 -62.5226116950 -32.8466858891 -28.7434331025 -37.5334108272 -26.8972012068 -24.2769621500 -51.7676620294 -29.2473523878 -18.4193400061 
		 2 	-21.8550286976 -30.0893141771 -34.6831941804 -46.7690794828 -65.3636689753 -63.6004587376 -26.7795812533 -18.7216275235 -21.4369707725  -8.1857490524 -25.7426905040 -26.2574965177 -38.1886241603 -44.8976343945 -49.1866797721 -27.9845586813 -42.4175351105 -49.4819322171 -24.5605416265 -38.7477796859 -44.8833291638 -23.3128607589 -42.6877765959 -19.4519894825 -12.9291231245 -10.2043729689 -48.7746775786 -24.9629167873 -33.2143455931 -64.9416480225 -31.2221947509 -44.0666255336 -46.3632790675 -23.1981664434 -34.2372866521  -4.2107130668 -35.5595535338 -32.2406989687 -48.4509525853 -31.6712217001 -19.0103527395 -50.9822241146 -36.8493558576 -24.1212571591 -51.7769342706 -27.8696703621 -45.4487423757 -17.7143887632 -43.8595565242 -45.0187141361 -43.0495718121 -47.9886675021 -25.0985055932 -71.8607782150 -23.5157559784 -19.0781988122 -28.7286366132 -12.8575893861 -59.6645442379 -14.5907021641 -41.2799267591  -1.8287753858 -29.4531080773 -64.4496141481 -15.4921822388 -45.4553112217  -3.9499111008 -65.9890814604 -39.0851304662 -45.0686340586 -36.1488203242 -38.1071425981 -19.8343789975 -20.0102074477 -24.5004808306 -32.4422509911 -17.5892535741 -26.0508651099 -33.7597246753 -56.3014089640 -29.2388122889 -18.3205540136  -2.4593823696 -28.7868320436 -12.0820568575 -23.7991790833 -75.1257251400 -44.5006282308 -29.7098067511 -26.1661744577 -19.1270375988 -39.6166843496  -3.2957445369 -46.0814487162 -34.3421765456  -9.4782858992 -24.6318228824 -36.1253612241 -41.7331383552 -29.9351521979 
		 3 	-20.2921815052 -18.9802189944 -21.1628874256 -26.1864568217 -44.6121694783 -26.4059830749 -39.0719972801 -13.0049092893 -30.1866157966 -46.5533886608 -40.6875456723 -10.8464815011 -45.1499267263 -45.5985319458 -45.0224226039 -34.8761143364 -41.6931853205 -57.5374653278 -20.9458093720  -6.6189553461 -30.0410004095 -32.0651191730 -21.7024567290 -43.8026215960 -28.9691233395 -11.0374085841 -42.2278928298 -40.9715273504 -21.2131367141 -51.6617705765 -21.2316340541 -29.0677756653 -56.2510575939 -32.8841870817 -44.1037192182 -12.3902310218 -51.4548759008 -40.1604593739 -41.5001162421 -18.9537229227 -62.8097765416 -33.8542880340 -34.4696763561 -20.6989802670 -19.2652627220 -56.0174518797 -35.9172677049 -35.3603754561 -51.3604351413 -44.7176140704 -46.3512560472  -4.4550637550 -52.7075261972 -37.4372909343 -18.9528328564 -37.2275313397 -53.7198898851  -7.2729962184 -43.1591294875  -4.4032849495 -31.5522300627 -34.2684146991   6.4446991760 -33.3392735473 -43.0189709832 -35.3147877102 -28.7222109979 -33.2830447456 -19.2564473218 -29.3883117606 -35.1783714304 -10.7241270646 -38.0949174663 -10.9870145852 -33.7372745357 -11.6812367792 -59.6438393400 -42.3481063790 -11.8702170220 -14.5106359048 -22.0974149186 -64.5875445983 -45.7075511497 -47.3527317065 -29.0061700369 -20.0655229708 -33.2546394898 -34.0417876861  -7.0203717549 -28.9986051523 -49.4424054893 -19.9747328069   1.7464138706 -26.9751382051 -28.6533684419 -23.2285617268 -26.1332348877 -29.6775402633 -13.8344524250 -15.9538817055 
		 4 	-27.6910611113 -28.2388950470 -37.4497144847 -11.5297456929 -30.0833429503 -44.7967369749 -31.0788479876 -39.3408844197 -24.1548936041 -34.8371019380 -23.5406995684 -27.4707614251 -43.2494997600 -44.8195943166 -13.8116666304 -21.9510282964 -20.4597379763 -24.1229875880 -36.8432511173 -47.7626373786 -42.2977457119 -35.3520921104 -37.8339977862 -53.1916515255 -35.9026849234 -37.5198463042 -32.6671461517 -59.1996047737 -45.3746256605 -40.0311960615 -42.5331439434 -43.9154688555 -22.9436713549 -47.1909199448 -31.6969624471 -23.8491744469 -31.7082198345 -48.8983740613 -20.2947520446 -46.4634358848 -45.1876558400 -37.0012710299 -60.0987727985 -62.7323384816 -18.7326625574 -32.2729100440 -71.5144117008 -13.4911669466 -29.6834505035 -22.3352657082 -15.0543231082 -22.3471405342 -14.5156354442 -28.1292557792 -33.8543987481  -9.0703625757 -40.0744067373 -29.4375704903  -7.4555597090 -16.1043465875 -39.2286903911 -45.4239639772 -27.7598965825 -36.0619642690 -32.4031485852 -18.2756079158  -4.4732673215 -52.5656698456 -44.2219443438 -48.0713151232 -46.9520524633 -12.9626214611 -70.2677942884 -41.3929404722 -43.2269524484 -37.4931161275 -42.0605418929 -45.7100070180 -21.3706965762 -39.9651042108 -45.0031020463  -9.1402698671 -62.0849699082 -31.3267603660 -38.2760562684  -1.1832623440 -34.7967754627 -36.5040737058 -17.8954298066 -21.3787552520 -38.3760012594  -0.3124630484 -17.6945484649 -13.2346108542 -38.3732949562 -34.2246337312  -6.5150916735 -30.5339845984 -34.0728366389 -25.1793413428 
		 5 	 -9.6705791481 -26.9432748661 -39.1034704790 -30.5785746630 -25.3987926591 -28.1891945055 -29.9006216743 -30.7491746631 -23.8388693232 -32.3975544674 -39.6765599607 -28.7676982980 -30.9379550451 -22.0969886614 -36.1484881958  -9.4126791282 -20.5052932805 -29.4518660610 -55.2714985027 -37.0243622520 -18.7849770590 -22.6195284293 -35.7509846322 -19.6096612364 -43.1171446936 -50.8530771428 -30.8404842113 -36.6067222950 -36.4067718853 -22.1122529008 -43.5705831713 -33.3980220451 -28.1830411106 -11.9382610500 -41.5421699334 -30.6409879140 -39.6840453304 -33.9332665010 -18.1096398500 -73.1509186900 -19.0037620674 -15.1153074003 -49.0508334498 -30.0398744435  -5.1555821936 -35.4394717200 -11.1238702385 -26.1393262531 -62.2398672411   1.0588163472 -15.6518595762 -29.7351022152 -16.8044925407 -25.7854564791 -26.7587827281 -45.5538976933 -29.7537754485 -42.0250014084 -35.6724179221 -57.6009614782 -32.5130478448 -23.7913627080 -25.0197640596 -38.4689975493 -58.4506841612 -10.4442756438 -16.7913213655  -9.1280679679 -34.3831851388 -32.8517328139 -10.9000471962 -24.1515809512 -39.9271876556 -61.1212782067 -14.3189424129 -34.0306633277 -22.5414674493 -39.7126900395  -1.2835036946 -30.4353156016 -18.5467924561 -29.9689922620 -73.3822817403 -35.4111564133 -20.2889827237 -28.3802793405 -39.4611177127 -32.1485789303 -58.6874479397 -22.8250285161 -36.3897796350 -62.8554526576 -32.5195487473 -30.4502791886 -69.1545717305 -37.3486925960  -9.1649062962 -33.8095103310 -38.4039091506 -40.7820054560 
		 6 	-45.2308178558  -9.5874506210 -32.4969822299 -18.8733415897 -19.5094724415 -21.9935804807 -40.9640446347 -33.1880306301 -25.3386378888 -38.0279249543 -43.0700232592 -33.1670946912 -61.8943433993  -2.9188295932 -20.3608291080 -44.8214921296 -23.0021842940 -74.7310316529 -51.8782835634   3.2367520746  -9.0055472138 -37.0384749286 -16.3948956965 -25.9889654071 -35.4623185293 -35.3375206242 -47.1544613043 -33.2393722725 -29.0821315928 -24.0930670336 -49.1523122682 -27.4049198795 -55.9000342001 -10.6051727855 -43.3188219936 -41.1290301305 -20.8854271569  -3.1782452184 -49.3236333065 -33.3817327227 -40.6406012418 -33.8901345427 -29.6488585455 -35.0819503968 -34.9151838176 -33.8109061084 -27.3929018894 -41.1627982808 -52.7489061766 -57.0397983004 -32.0560715224 -30.1258720295 -33.3554655684 -27.0532211123 -16.7155225051 -19.7832531877 -20.0953243732 -26.6734864023 -26.3749851646 -34.2071173756 -14.0968002687 -33.1479363293 -45.2135803651 -59.8277159770 -51.7442447530 -25.5098091172  -5.3639888338 -48.4697546942 -52.2543840248 -25.6592094874 -38.5436297096  -6.9640929854 -15.9878964630 -23.4555906515 -29.4444358744 -19.8659244588 -51.1582007533  -7.3091783338 -20.7141917064 -30.6722534059 -50.2610860322   1.8238408478 -52.6961907838 -58.1312539692 -24.7767593589 -28.2078682980 -26.0131917959 -36.2776778034 -37.4572860548 -53.4646968582 -29.9155293397 -29.0981110878 -27.6779093503 -45.0139031206 -27.7065177838 -19.8661618869  11.9549532275 -46.9772325711 -14.7352275117 -13.3315233627 
		 7 	-33.3448193958 -20.8787561536 -28.4944207222 -43.0548366950 -29.5543922115 -13.7369074104 -47.3514807531 -34.0090327690 -57.6054978567 -29.8136833572 -19.8696251773 -43.3446394072 -20.0601357340 -26.5996654171 -12.9534083415 -19.7940552530 -29.5273629568 -27.1728201711 -42.4721481998 -45.0846660441 -40.6615932210 -33.7063529765 -29.3941835290 -43.2571377214 -35.9910489248 -13.8159694866  -9.8631064593 -39.7818818052 -33.9344729059  -8.6089095826 -37.7111195452 -22.8377639564 -29.9750981419 -45.7025065849 -20.2447709091 -43.7624479220 -16.5867039243 -15.2533880633 -32.8115580664 -23.8207465304  -1.0873076341 -37.9284750089 -51.0519958866 -22.0818871781 -38.9526821754 -16.7731278917 -23.1383961256 -16.7770408376 -37.4914582201 -33.7241990457   0.6418448061 -47.4697494124 -39.8033452616 -61.4539430311 -37.1304884000  25.4539962824 -59.2880469684 -35.1592113342 -10.8484112372 -39.1967700328 -30.9507163317 -50.5321420600 -61.5148264253 -33.1229792638 -27.2194713494 -14.2347390720 -29.9564281079 -50.2298835382 -46.8730004619 -41.5527149383 -45.2663634701 -47.3849151466 -21.2092147097 -30.2835521455 -24.7525650447 -34.3928512834 -47.5951065474 -28.3131491058 -33.1434173318 -31.2219704544 -36.9478885435 -25.4018109356 -43.3993911006 -12.5877686021 -59.8506128930 -33.1953358807 -29.8172549698 -14.8671384252 -24.4468089236 -22.5384882910 -47.7287603283 -42.6484744294 -43.0049909393 -32.8798077756 -37.0408380066 -37.9887249134 -47.1963140255 -32.9802424031 -24.1875293835 -62.1712538358 
		 8 	-26.3082390238 -34.6960823584 -22.4490413930 -11.8664140646 -54.6590414738 -38.5210439830 -61.9321650171 -11.3378104156 -44.4668308547 -38.0785720208 -18.9058235787 -31.7515961976 -17.7372618313 -20.6214986357  -2.7376103239 -46.2515171407 -14.4170733269 -27.6568952990 -22.7222631790  -8.8070579715 -65.4057019137  -5.7157991416 -58.7884284629  -7.7245937402 -52.1261971413 -59.0383219688 -19.9630218481 -27.5535171158 -34.4511087809 -27.5000016747 -28.1605900672 -49.2532881277 -40.1828444837 -35.0642114531 -29.4487443507 -11.8786936146 -49.5622672173 -28.6450598062 -31.2905207837 -51.8146071601 -11.9577747978 -43.8340320984 -24.6046228637 -21.0117020453 -47.7036773361 -25.6687373674 -34.9166239702 -17.8959838774 -30.3060759544 -35.1004689513 -39.3683605335 -42.3984080993 -24.5251093027 -12.6847634399 -45.2993446991 -32.8355483849 -56.0815288903 -17.7001520430 -17.2954574592 -54.0818222156 -12.3232685947 -44.4419132204 -30.1208907749 -48.2098153285 -61.0463072460 -25.4769079422 -29.7691432417 -37.8869646907 -43.1873159081 -14.7693401907 -27.3057013954 -39.4328343186 -47.3261558561 -26.3277364791 -23.9388936108   1.0503972901 -36.6111563261 -17.7310532236 -35.6453523885 -38.1382293429 -50.9837792335 -19.6292645062 -70.1094424981 -30.7911053150 -50.7449396131   0.3029648564   4.2849924095 -53.1940286670 -60.2307217916 -44.2733571579 -57.3312592652 -28.5180381119 -32.3121451472 -31.4659911232 -47.5963401594 -33.6185822496 -76.7114068811 -50.9829309746 -59.9386300879 -48.7794365354 
		 9 	-17.9597300208 -38.8835598414 -19.2501363353 -20.1939839072 -29.5514097114 -43.5327015321 -47.3690187473 -38.9725137301 -59.0159788129 -27.9392835633 -35.5776585954 -58.3164761957 -48.7116413435 -41.4314488457 -25.3502201739 -34.6528129681 -49.0202250265 -44.7948555458 -59.1727577468 -25.7405173431 -30.7960543491 -27.1335368987 -31.7839571119 -20.4438731441 -61.4201765960 -25.1587980324 -30.8275761430 -30.4510135523 -32.3697965992 -33.3027782727 -30.2934125277 -26.4964044400 -17.2764984957 -37.6221497070 -26.5444037211 -20.8518450578 -44.0175858828 -58.8795954947 -36.3896710231 -39.0859667274 -21.4117832183 -19.6935960738 -37.0666956925 -30.4697010919 -54.8057502624 -25.9818766907 -11.8054318786 -41.0870344392 -51.6860415857 -52.9103772192  -8.8786070545 -16.6824702622 -45.7345889386 -33.5462706636 -31.3082287463 -33.4392955471 -29.7476304469 -31.9626225519 -32.2667744079 -23.6784062652 -59.1559369348 -19.2236815186 -37.7221979653 -28.3131491058 -32.4851272560 -20.6757555126 -31.7708598435 -23.7076684564 -13.3166769847 -35.5193176065 -62.9579519002 -27.0693931835 -41.5322766638 -51.5046495061 -14.7998107786 -34.4132687433 -24.8652943927 -17.9006907560 -41.2824614860 -23.4748743999 -37.2898003817 -40.7150409076 -10.4297182049 -22.0050618904 -64.4277092564 -23.2466451458 -55.2416461218 -41.4840287221 -12.5814592358 -35.4112216736 -35.1852122300 -27.8612686258 -56.3930126339 -41.1888609013  -8.3127218954 -42.0777178813 -30.6313196381 -48.0943095018 -27.9633790688 -53.7039426364 
		 10 	-22.8619556400 -18.2135301395 -33.7501644124 -27.7431844555 -40.3662419116 -30.6167107877 -33.5722174587 -35.3649142389 -39.2244673066 -44.9489752697 -34.0989597228 -11.5454781037 -51.0208821743 -61.9290598290 -39.5718441787 -37.1816115435  -7.7231884277 -11.1413189923 -24.5726286648 -29.1556836357 -44.8076131236  -8.2754208909 -47.7232039302 -15.3968330346 -16.6707446826 -49.5348337635 -33.8603531277 -31.2952957705 -23.2967566190 -20.0296861981 -18.3857326875 -27.4624010615 -49.2457635882 -17.1318103589 -32.7737326361  -9.0688846483 -53.9428836769 -50.4261652336 -37.4449474098 -32.7521951138 -32.6323850581 -31.6859867684 -35.1290025280 -46.4680214733 -36.7260496932 -29.6534823565 -16.6260052570 -17.5853910093 -16.0937397280 -49.5832663331 -26.9575957976 -40.5384125330 -50.0547557180 -42.2811296862 -45.1963097800 -32.4776320681 -12.7296394671 -18.7406369837 -24.4068424575 -23.2901639098 -34.0846484752 -35.3731243396 -11.9430406160 -38.0337536736 -30.4618708979 -41.8392191577 -34.2430427575 -20.9111199038 -27.0333933634 -37.5612506265 -34.5735837050 -23.1291369832 -15.5596641622 -35.9081776828 -18.3494007903 -41.9868700217 -44.7729019047 -22.1868497899 -33.3342291952 -32.6752060099 -26.0418067441 -37.9337201978 -28.1894205270 -14.1625237945 -21.1598509015 -16.3950527564 -29.8131366197 -40.5493495771 -43.0194188729 -13.4956949334 -68.2237862997 -51.0076380022 -35.9058042435 -51.9101496249 -37.2157163177  -3.4564223621 -37.6101051278 -47.4883232453 -29.7562811256 -31.3328866355 
		 11 	-46.4925911439 -47.5935075454 -42.7642213604 -48.6319870662 -10.9348967718 -49.8238606067 -36.6094034159 -49.3749507822 -71.8144624312 -34.9464670358  -5.1876233011 -38.5915513082 -35.5253253283 -20.8543370117 -33.7731389735 -24.7480284839 -47.1446018311 -34.3423416573 -42.6132174686 -45.9294419159 -54.7488180819 -33.0756112995 -24.5237591500 -54.3728177221 -16.7639554211 -15.2672147394 -45.7490830034 -33.6712377355 -29.2875188586 -22.9938756514 -39.6831908970 -18.7578326885 -64.8849408595 -27.6765430030 -34.6796797506 -14.0405114654 -41.1759001188 -16.9168869543 -28.5012285625 -30.8331011876 -32.7967062431 -28.9836748743 -24.6969693925  -2.4538245046 -46.7887084415 -47.5389209359 -38.5879088264 -50.3185306733  -4.1538389652 -30.7872180936 -18.5002542376 -13.9899539675 -18.7587105095 -35.8767538281 -45.5420998117 -42.6154104369 -36.9282346674 -36.9566866764 -23.4675769201 -47.6182765939 -32.9111298475 -35.8824508882 -39.5318140527 -35.1371669183 -11.0598030514 -32.6253053918 -38.5288332941 -28.2995509169 -10.3476031724 -39.0435354321 -37.1711279942 -28.6169915950 -14.2081214932 -42.0897682749 -50.1076221768 -23.1381542027 -55.3860731964 -45.8273704387 -52.7449451750 -30.3833359347 -52.6661624983 -25.2804865078 -33.7714653343 -12.2083600332 -19.8015672905 -27.7963762623 -25.2363892376 -40.8866331249 -44.2071448310 -28.0662589636   1.2773403884 -48.1340788590 -33.2807552043 -25.7803865863 -42.7819563313 -45.6307082607 -24.9020627816 -58.7830335031 -22.7514702305 -10.7037871415 
		 12 	-23.0829467065 -47.1864212347 -56.8143200950 -36.1104345235 -32.0675491140 -31.6020131189 -40.1170603219  -3.3265002303 -52.5092520202 -55.5554445083 -32.0638036248 -57.5755499148 -30.0412869224 -36.0252776231 -25.7395188563 -14.2926980426 -29.7486015410 -48.7687377785 -22.7820114647 -28.3390434513  -9.5825587035 -17.1862642324 -29.1929734450  -8.6660088009 -68.0509113469 -24.9246220650 -12.3122759158 -28.7926137005 -27.8926475408 -36.7622364303 -44.0829820372 -12.7383166135 -12.8200157218 -27.2110095449 -31.2560895109 -12.8494240595 -44.5435557687 -11.8843613208 -48.1457173582 -23.3743434741 -55.6012072705 -20.7686213007 -31.9580738010 -64.6674714747 -56.3080546160 -12.5435652430 -30.6305377115 -36.1079401787 -43.2527736458 -52.6355805307 -12.6925399461 -11.2407725446 -12.8977566387 -37.8918487369 -54.9138567624 -56.9785805216 -17.8678356686 -39.5526925198 -35.9778529441 -49.0213394981 -48.5083185713  -8.7118383353 -38.8871674777 -40.6328998743 -12.6635864557 -16.1024894485 -29.4096964963 -21.0214628193 -43.2515675578 -19.2664550669 -45.6258803735 -30.8523216019 -60.5263173018 -45.4750719249 -32.6898504574 -27.1025173409 -64.9974345702 -22.7602169217 -23.3137260940 -49.9890372492 -37.1570286738 -42.3841672881 -44.4649007496 -26.5344916763 -32.7410325068 -36.4911728714 -23.7412184098 -51.5771464764 -37.5477108389 -13.6339514308 -27.1563013019 -11.5086942889 -40.2406136086  -3.9594114869 -32.8434963360 -49.5446405727 -28.0450765164 -19.6302891736 -43.4091800655 -25.5816364801 
		 13 	-26.5868665035   0.7793617470 -32.0695320595 -29.1901886935 -31.8768338377 -43.8911047115 -39.1691975469 -42.8167009319 -16.4321814843 -51.8381192741 -40.9182248413 -29.7966634145 -21.7309665118 -37.8251176477 -54.3663901627 -44.3402736920 -45.3709535282 -20.7810108945 -42.2683995411 -39.0744108200  -4.5288797621 -39.9651596897 -23.4422096080 -33.8806702400 -30.1938839118 -58.0753220835 -16.6516755334 -51.9986215476 -25.3874023228 -62.2745632790 -29.4129532191 -19.6004741183 -43.0042664930 -34.9790948301 -54.4154556245 -10.6959433049 -29.5751325375 -11.3625049508 -43.3933972385 -40.0149850121 -15.3841211464 -28.0806156697 -27.7245263177 -40.9409981604 -63.1156087492 -26.9941557535 -40.3511951420 -33.4404042105 -51.8808891541 -28.9718409363 -33.2665605951 -34.6985406129 -16.5459236053 -25.2465157382 -20.4434126758 -15.0126465190 -42.9829278219 -37.6297249828 -40.2297025458 -35.0393275567 -14.3903331731 -50.3483779654 -42.3344573222 -44.0053033424 -41.3880275774 -24.8237941237 -28.8575197749 -28.9634285908 -21.8779772921 -26.4997458656 -39.6303922265 -46.6026978979 -50.4055683890 -39.8763624882 -23.6944193720 -20.6831438288 -33.4321654117 -43.2507304316 -39.2235600287 -31.7726910949 -33.1173563039 -34.9570569017 -39.5390259274 -37.1182625530 -20.9267179318 -40.1805201995 -57.5448395630 -48.7236054965 -33.5256246101 -42.0127659975 -19.6325061085 -26.3717389207 -30.0179749366 -24.6370526039   0.6418448061 -38.7458453475 -45.3890128573   5.5626722923 -46.6750812228 -51.8535979409 
		 14 	 17.6359632244 -19.1550911091 -27.7186923556 -33.4826339428 -37.0667169117 -25.4741217781 -46.0000074857 -30.6827542225 -23.3111222714 -36.3142097282 -26.5319755996 -37.3478915412 -44.9267254869 -24.7955269894 -34.0113845154 -11.6181108891 -31.3994020131 -48.9276350515 -19.1017188983 -21.5584625653 -29.3570769188 -31.8402574365 -38.3204910877 -16.4243976529  -2.1253091646 -50.0824846171 -38.9457882822 -32.2755183299  -9.4788343404 -21.9104077664 -40.5710493470 -28.0933433713 -14.6598862439 -21.8015862880  -9.0892847552 -32.9418384916 -45.4725379392 -27.5451303288 -25.0403906655 -37.2483889781 -48.8798810219 -37.6797351001 -43.3939176451 -55.2069749430 -39.8394728220 -32.3575462889 -34.6414263381 -19.2015568283 -18.5597101336 -26.7465807602 -50.4310062518 -21.9928979316 -40.2709315701 -33.5199192466  13.4022364174 -28.0596640503 -48.7698430903 -48.3975438952 -25.0525667640   0.4503846824 -26.9539092250 -60.9221962388 -39.2446891533 -22.6562348891 -56.3261184826 -52.3956934000 -46.8468361606 -32.4730300083 -65.4040281546 -29.8731321589 -22.6011638546 -22.7041013813  -1.0172000311 -22.6529565737 -24.9237535914  -2.1457816274 -21.4447297254 -50.2493867304 -10.0884550544 -28.4472150205 -25.5849502758 -47.4237092842 -18.4466595517 -21.6631091690 -48.9549356102 -21.3201988869 -17.6184396937 -45.7082101469  -7.2293333798 -16.5603442291 -15.1479797936 -47.5994535317 -38.4438643981 -32.7588589197 -34.3707764106 -24.9293585662 -21.8203566098 -32.4344237641 -39.7782837687 -33.6324993708 
		 15 	-42.7628877842  -5.1223203387 -40.2156828078 -46.6883468237 -26.9612257657 -18.1881234580 -51.8484638285 -25.1336723968 -57.1220174958 -39.2159085448 -34.3545141922 -30.9830299179 -29.7565901312 -38.8221993732  -8.4429760470 -16.0705118371 -23.3293726393 -26.9919335957 -25.5033127174 -40.1149416457 -26.3284521463 -48.1007839042 -36.7146936986 -25.4369236967 -43.6124130154 -43.8948192050 -27.3080395971 -36.3016862107 -39.3671426322 -15.9194703013 -11.0247241731 -46.2823655391 -25.3742759677 -20.3182610500 -22.2974583364 -32.3291741736 -31.1121161596 -57.0858567308 -27.8040254039 -35.6989473161 -45.2099512166 -33.2264970554 -65.3911490148 -29.6602922288 -20.7630157902 -24.2049472560 -46.9314820771 -51.0735633081 -32.9019238615 -64.9485712504 -31.4498632847 -18.4085997298 -19.6572063255 -56.1686113826 -37.2708297059 -25.3676418547 -36.5698553454 -41.6936651694 -24.7546737094 -29.1442741619 -34.9032868353 -14.8392504029 -38.8540002470 -42.5047472798 -24.3457163236 -25.3881511159 -33.3489632150 -27.3015726266 -36.3836882598 -62.6007082071 -48.4378505101 -39.2821973545 -34.7717792895  -8.2352208209 -20.6442792368 -36.3781417919 -25.1124081789 -60.4515728529 -32.5572114845   5.1088250127 -27.9467003112 -31.9268198328 -26.8212974287 -57.3540319416 -39.2791093215 -10.2698928869 -56.8174354519 -50.3210250413 -18.4538488655 -21.9775391035 -34.3411543042 -24.7892200607 -39.0636528306 -41.4582339738 -64.9707959056 -59.8985016621 -41.8983126789 -17.5700543669 -40.1450256775 -15.7880359290 
		 16 	-25.4103377693 -35.6696683503 -49.7553488084 -10.8334820332 -43.6479828903 -35.7685511118 -14.5704864729 -18.3481659410 -54.0868261689 -45.4111604095 -20.2231477616 -22.2469539463 -48.8760603077 -16.7348457614 -40.7375044123 -64.6389068042 -43.1124930787 -48.0489460855 -34.8977710661 -17.0536610017 -39.9001022319  -5.9598277833 -48.7314732331 -21.4222699215  -1.9717091813 -35.9978270556 -19.5911481825 -29.4175558113 -42.6515833016   3.9102947530 -87.8548493596 -24.9983489698 -28.2917114314 -14.2575014684 -24.5541418881  -8.0874600804 -52.6546475447 -36.3798553579  -8.0811335654 -27.3108998365 -19.7599775900 -45.7048000139 -40.8111774132  -7.4739527030 -30.6778548101  -9.1943631101 -28.9491504952 -33.8564012043 -23.0658850431 -29.8764693979 -53.5093154987 -59.5372602593 -28.3507489246 -46.0876641980 -23.5940049411 -28.9915137982 -17.1912898587 -48.5183591684 -24.7547373987 -27.4483032379   0.8327665710 -35.4272091000 -46.2280948993 -34.8055279688 -50.2715042089 -24.3744678315 -23.8246140969 -37.2571303213 -33.1738656789 -43.2571377214 -34.7415760610 -48.4820061405 -63.3550743839 -14.1304538808 -60.8799827411 -28.6180871521 -36.5459285659  -5.8374847571 -50.6633760283 -16.9910713646 -28.4233641615 -36.9276268966 -51.5222882031 -18.4270858589 -35.7582612996 -34.5039861751 -12.2123075702  -7.9020058891 -30.9679577539 -33.7894919200 -47.9581118786 -29.6757967544 -20.5017887961 -25.4075167885 -39.0986732873 -31.9093655777 -34.5321203627 -34.1262093205 -12.9520094537 -17.5849612762 
		 17 	-41.5939875424 -28.8058448092 -37.1031825595 -37.5497818766 -13.3817520030 -37.4115280646 -24.7663457326 -54.3143594661   2.4563338084  -9.1980304956 -40.7553249776 -26.7884920294 -27.1374904728 -16.9373582565 -19.4576886157 -40.4198486091 -17.2886419145 -18.0081652374 -17.5964303659 -23.8603043697 -43.2619118756 -23.3190266243 -47.6233328478 -51.7345447533 -30.5803440566 -25.0586598511 -34.1036520174 -36.9600307709 -31.8072155851 -38.0382428425  10.2755002957 -47.1958891417 -27.7038583850 -36.4407695086 -48.5731188852 -11.6056321524 -40.7442941865 -53.3272933282 -29.6806205445 -35.6743087036 -35.0991519563 -59.8490318587 -36.0889641965 -17.7621137707 -31.0668851528 -26.3939872296 -41.8906332507 -31.1908681089 -42.5943442636 -32.5061133687 -39.4551819038 -18.8642464726 -43.9636712456 -23.0059782756   5.2864507457 -24.5847842513 -48.3113517190 -19.2494148476 -31.3991856184 -30.6262290084 -21.3329794497 -57.1215949765 -52.8523024848 -24.5165606368 -23.3637618819 -42.5604050582 -44.1267674628 -33.2701757925 -41.3820949839 -38.1099367721 -56.5618479277 -20.1183157017 -39.5044655866 -39.2651900337 -55.7687728258 -27.5238699017 -31.6071381672 -50.6530680883 -56.9352684945 -14.2487511137 -52.5840364838 -29.1186271698 -11.8617994691 -26.4116973339 -37.9747506773 -44.5113468496 -64.8381802646 -27.6530969094 -46.9298109601 -49.4112886092 -17.9664720448 -29.8703640607 -14.2034907377 -30.5735353231 -14.8967092641 -16.1957484420 -15.2536591032 -32.4727393864 -55.7705339092 -34.9856776779 
		 18 	-41.7888805835 -41.2997855435 -49.2212656617 -32.3499920505 -59.4410703623 -47.0986677961 -22.5468627244   0.2822096743 -40.8609060123 -34.4164857725 -42.5297831877 -46.1751398460 -58.3068609852 -32.9947249184 -10.6572105266 -19.9050569814 -50.0669103535 -44.6554078801 -43.1395591642 -35.4900900943 -17.1060256971 -52.3785722908 -18.9346231385 -25.0179913234  -2.6137343324 -54.6840213597 -33.5129931164 -35.5781892365  -1.6380589135 -19.3755729395   5.5068246600 -36.9604726800 -64.1351892947 -42.5980004693 -24.7614362958 -42.8887456005 -38.8694287709 -33.2880054939 -44.7738712196  -9.6455335100 -33.7499329467 -33.1808239986 -43.3525238596 -16.7415043993  -0.3165838819 -40.8022982908 -44.7167416352 -28.1924349608 -25.3655256269 -48.0110301676 -26.0427402407 -35.6987943075 -41.4291349581 -44.6968991430 -50.1169218569 -33.1891018545 -47.8928016390 -35.8603102965 -34.1588205103 -12.7811068582 -38.4795402132 -25.6230567826 -53.1000513988 -26.3034672084 -34.4858735660 -15.5166351840 -23.0068487161 -54.0217630052 -33.6276285456 -40.7909946890 -46.9352129368 -39.8353022563 -26.8618028174 -43.0434601982 -60.3216546695 -13.8851031579  -7.5913346926 -36.1899942322 -10.1341430971 -41.5822682055 -13.6891493023 -46.3981622953 -41.4365295853 -21.9434476884 -46.4639632651 -29.2697556515  -5.5558694977 -20.2173455587 -51.0915004726 -30.7268408164 -39.0404316205 -22.3226274183  -9.9745480446 -41.3344785660 -53.7575422548 -27.6122108457 -37.7286104220 -42.1643194352 -30.9846294623 -29.1231040721 
		 19 	-29.4380733658 -50.1506161068 -14.5798964431  -7.0028354892 -24.2651658049 -29.3869564503 -45.8135742015 -47.3016224360 -16.4031863283 -45.7494505418 -48.2558272749 -41.7461633494 -32.1585674000 -19.6408287480 -22.9451785007 -39.9443658772  -5.1051660036 -37.6510543852 -26.5430121042 -48.6327716449 -16.8113113724 -23.7708690538 -54.8269898251 -48.5420395849  -1.5506667876 -27.4309213230 -26.9934602841 -37.0083105897 -56.1887204269 -17.3665956583 -54.3442086558 -17.4127389219 -39.6163981621 -53.4392364810 -49.5509414836 -37.8944066278 -27.6433139540 -57.3063686161 -29.6377987462 -21.9528209112 -28.8712186075 -18.4325161093 -36.4885912866 -39.0860170282 -37.7563529618 -51.1520919142  -8.9266026065 -37.9674046074 -38.6713527743  -4.0919511784 -37.9719206307 -37.6409614614 -23.3342631855 -32.2727608670 -51.8268035083 -13.1337388297 -58.4252784844 -59.0014339000 -56.7973663843 -32.3801530619 -23.6434082145 -43.6278572989 -49.1906434648 -34.6556160409 -38.3584545674 -58.3268664935 -24.5464121843 -13.2361357899 -28.5563371011 -48.5652767039 -26.6652113819 -48.8601466604 -22.5075392677 -37.7432773579 -45.7365030365 -22.2952540344 -46.5294958869 -19.8059138725 -17.1055612732 -33.4425225851 -31.0758553686 -33.4255493497 -22.9802134739 -21.3650861857 -15.3091850937 -22.3505159401 -17.6117545670 -29.4039712715 -36.8932286297 -25.6226881077  15.2912245915 -21.5341007523 -49.0084780486 -42.9932036065 -34.8617140429 -33.0641178517 -44.2080924173 -34.1864973763 -20.8635174290   0.0450557207 
		 20 	-15.8516474343 -49.8237432166 -20.4248642567 -21.9142826671 -63.2392677137 -16.7665225748 -28.9697388213 -30.4653985665 -56.5695551967 -23.9371311461 -26.2693715173 -41.1038810909 -28.5056722365 -12.1013838271 -17.4658842153  -6.4613743605 -29.3340649331 -19.9649483816 -52.8460124224 -39.8327658271 -21.7006318603 -38.4874366862 -33.2494685091 -24.8470993811 -17.5012215219 -38.5791162459 -19.1125526829 -38.8886466584 -30.7578834382 -42.1984856774 -43.4347114533 -21.5883676394 -25.6852504651 -29.7301791409 -30.7995030979 -30.6753099995 -46.2762822124 -32.7338123090 -28.3002390977 -19.0297659977 -27.3824683423 -17.2880808858  -8.2633216921 -29.9140686032 -41.0776208179 -32.6175905784 -40.7509912056 -37.0944801236 -24.6628381517 -32.5829324458 -17.4417755683  -5.2809991584 -30.2281891737 -36.8306334293 -34.5573610812 -27.1097238698 -40.0274308969 -23.6153769293   5.7481357398 -49.0747839712   2.0970884708 -29.5350960288 -41.2975529839 -31.2248476553 -36.6315878949 -43.3493057335 -43.2225943717 -24.7822539749 -20.9934704618 -23.4561739486  -6.1331905105 -21.9853066215 -67.1921555566 -25.3833379535 -24.9080395818 -32.1459912755 -29.6353665032 -28.6680069535 -20.3819585235 -12.4009971924 -37.4224099187 -46.3102561758 -46.4205042384 -28.6955754784 -24.9660045114  -9.4095083788 -32.7920025955 -46.0339309747 -15.4169117240 -33.5050526565 -32.4228405578 -22.6434559971 -25.5304959160 -46.6998112644 -39.5948124298 -41.3955385638 -29.0487458924 -40.6816126801 -27.2327129491 -29.0896825464 
		 21 	-39.7506328959 -38.5970312874 -45.9818595599 -47.1528647350 -34.1373707241 -47.5167867961 -40.2986771967 -27.2125389260 -47.8839425018 -38.6318262110 -25.3462013676 -52.6593902618 -52.6296935885   3.4548924208 -11.9081859881 -48.3387329294 -15.1321806421 -32.7675173404 -67.8977520439 -56.4509273249 -33.0796646070 -45.9299001754 -38.4215226275 -44.2299208269 -50.0159265822 -31.2255595380 -35.0039600046 -24.9526896795 -43.4138938656 -40.5702207392 -55.3452310902 -34.0702973595 -13.6652703342  -3.8205251760  -5.6067663968 -51.1190431899 -51.3950206558 -32.8129704511 -16.0605902931 -43.7914356744   4.6777213676 -24.1686965068 -31.1091262578 -42.8326374112 -55.2621840159 -43.4739093402 -37.6494142204 -23.8858359673 -15.2036513731 -28.5586697622 -20.0762437595 -39.2910093471 -19.4424458780 -31.4340462102 -23.4695984321 -32.4478134873 -45.7147032011 -38.5714240533  -7.2114482381 -29.4653513269 -56.5493501215 -48.3022408205 -26.5135066420 -34.7888453113 -50.0926618887  -3.7408274324 -44.5496146592 -26.7289881210 -16.8799543880 -42.2608385449 -33.9263204142 -42.3801832180 -20.3874247704 -34.2529538765 -33.4867666088 -41.3389685719 -43.5207791875 -43.7174685084 -19.6885594316 -31.2864051214 -34.6924039627 -41.4925762350 -48.4136862336 -42.5730538917 -39.7077353057 -29.1661670552 -51.4441246356 -26.8161158950 -39.7433253321 -29.6217804912 -35.4619918134 -23.8491744469 -30.0329557872 -60.1200020232 -28.7537578106 -44.1513946060 -26.0607430567 -57.4656455034 -34.5358089918  -9.4121005329 
		 22 	-37.3460134320 -71.8090344566 -56.9942331582 -10.5868718062 -19.7805560727 -23.7871097904 -18.0782577804 -43.5991919948 -35.9528211815 -22.7357612355 -39.9776249041 -27.5519842161 -22.7895569503 -41.9498157862 -43.1565920485 -33.4084663518 -42.4876280191 -17.5528095638 -31.1557160966 -46.0785662763 -15.1589256866 -16.0350119500 -44.0417919210 -35.1473199272 -45.1509132169 -11.5354007941 -64.7153423347 -25.5113007549 -20.8528887030 -55.1670115794 -22.0622434260 -35.4479350551 -27.7640421242 -21.5414762501 -33.0437125145 -30.2686432491 -17.6123525375 -26.7872074246 -52.2037108143 -40.0257692500 -39.4754549653 -29.0157904168 -54.6510982406 -17.3421856465 -46.6200407090 -67.3949584517 -30.0919215067 -37.0979113037   1.6057081149 -19.1356076178 -55.5799925131 -22.2348289551 -28.3056452526 -19.6812302798 -27.5450904315 -48.8869565357 -21.4919663227 -29.1329476903 -34.4901811899 -70.9740955403 -26.1224347625 -29.7545289294 -22.5548765428 -45.9658053751 -44.1935770712 -14.6342281301 -33.0247618088 -37.5829266462 -26.8224049183 -23.6572040057 -53.3999110048 -53.0665869244 -23.5919246042 -23.5909571012 -34.2818657660 -24.1885205926 -39.3173109088 -22.3820163171 -25.5258694619 -68.9533545953 -44.2330379916  -6.2816260074 -26.6228123395 -22.7432725987 -38.7894634513 -31.3556621284 -13.2754701363 -40.4831229248 -18.2503954449 -42.9748961360 -26.0053574594  -7.6211541890 -20.9583218317 -58.9165512161 -24.4076698255 -41.5455271349 -46.2228131936 -43.1557883615 -31.4346458076 -52.3403463668 
		 23 	-40.5637565394 -43.4049241548 -49.1317760224 -26.4574720605 -16.2749690405 -36.2995141182 -29.6058933643 -26.8925303321 -30.9598511505  -2.2476295904 -39.6197293404 -50.5529799554 -31.8454128807 -37.5543851522 -10.4024158620 -28.2383270176 -52.9126499914 -30.5707249072 -13.9443979564 -38.4814160882 -19.1911762181 -32.4062030666 -52.3969799117 -39.9193016285 -22.2292334523 -35.1645874432 -65.2610698687   2.4453585478 -32.1660731572 -36.4767050637 -34.7241435302 -21.7021229778 -48.3035450152 -28.2869543171 -19.8204079442 -52.4286103575 -18.5504608017 -34.9105214672 -41.0980165145 -18.1641393592 -40.0820560296  -6.8669592908 -38.0182035877 -45.9374897694 -29.0876821356 -34.0704054724 -20.8543370117 -25.3377667970 -23.1749965757 -25.8778322871 -25.8978563813 -32.0011111181 -38.2083032367 -28.0172105160 -27.5568411514 -40.1009486558 -27.6117583176 -36.7242155846 -29.5373100413 -38.6454662082 -20.8933541893 -48.5836174777 -30.3517739859 -16.6391383268 -11.9563981471 -22.7296966759  -1.2566268985 -19.5168034555 -36.4203120638 -27.8842960011 -60.5153648469 -56.6140065734 -43.6087161257 -56.5314664650 -43.3688927964 -39.8053197588 -16.0604624679 -28.8257244990 -45.0419926629 -44.8086768958 -59.3372505991 -40.9654671314 -10.7885057032 -33.4893846513 -36.9595995074  -0.1440036148 -66.9740007226 -25.3835323627 -30.0432503539 -21.0173878502 -37.5508371518 -41.7843433225 -27.2909717061 -25.8639036619 -31.8401014954 -25.9095710701 -27.4217070953 -59.7494638618 -14.8056904441 -25.9089806387 
		 24 	-31.3467383467 -28.0874305140 -42.0006878676 -36.2686573147 -43.3876493279 -41.3139697167 -20.8510574022 -33.7548054616 -10.2808012654 -27.3205975169  -4.2375481113  -8.9730996645 -27.4973409101 -29.7165257198 -14.5932558479  -9.2615712772 -26.9375906112 -22.2332679328 -29.2827308371 -38.8837659440  -9.2610769320 -56.1422664672 -15.5191055263 -36.8047236382 -21.6998109320 -41.3746476430 -42.1205567154 -42.3551721521  -6.9639327309 -38.7503817247  -3.8307755638 -55.6808845395 -36.6409483243 -54.4680407797 -43.8159356684 -25.6052350427 -26.5789478928 -43.8105737852 -43.8251365814 -24.4132428631 -23.8819716483 -33.5993459518 -35.6092147607 -62.2625984827 -36.6904213817 -38.1069201689 -54.0397298545 -38.9336142795 -44.7051029244  -7.9424062459  -8.5494173574 -62.9099937299 -37.8502366021 -10.3785528065 -27.9368901236 -39.1750717204 -14.4295253526 -33.2192018941 -13.7065403919 -27.5206761332  -7.9439736855 -36.2969398114 -33.4885104930 -24.0845019042 -29.1300455959 -35.3209151798 -52.2150257391 -42.4478924310 -37.7624649808 -14.7021245504 -12.1819433195 -33.6758826244 -26.8935190832 -35.0670619298 -37.5280690720 -45.4403246002 -23.0163022275 -37.8171872597 -41.8621198765 -18.5976789200 -20.8063972463 -36.0726092163 -10.3340347170 -39.3570929761 -22.8787318103 -56.6127810606 -39.2610567030 -46.7825120670 -44.4831446725 -25.9267617829  -8.1334632562 -56.3488003477 -31.9807920720  -9.7533255283 -49.4820838493 -60.4645237731 -25.8606341588 -23.8861549497 -46.2626498225 -43.0986151691 
		 25 	-29.0468406673 -22.8822268247 -16.3241867287 -77.7855021696 -43.4198547504 -35.1210988336 -53.9366418633  -4.8332361529 -26.8012369631   9.2268494247 -29.2407491113 -33.2852211337 -56.1386621557 -35.2512536883 -48.3019013916 -42.5786723430 -22.4600579536 -33.7499329467 -50.1844524370 -24.4438104422 -39.0611037627 -37.8816455872 -16.5015166490 -55.6225986788 -10.4442807304 -36.7992899853 -31.0356349297 -22.9489540635 -17.1995865776 -29.2497912646 -24.4058605195 -34.1812608819 -46.3208524652 -40.9248989084 -69.1133825965 -38.5109702962 -42.5823417888 -30.0625015930 -45.0828379751 -20.6661596862 -37.0460101021 -26.2484684262 -26.9185772303 -43.4839208041 -13.5145302543 -36.2517950233 -33.6183229507 -25.5880236377 -32.2335116587 -39.2593020309 -31.9834708185 -48.5208213534 -53.1421977122 -24.2966971551 -31.0700526950 -29.6232526300 -27.8477656901 -48.6433749907 -48.3909881212 -55.5846253133 -31.3843237413 -26.7295190924 -47.8262385214 -22.4683823053 -15.7582027111 -38.3606899151 -13.7759408959 -32.5243604330 -31.1368153758 -40.2835252149 -32.1261046130 -54.6498306736 -16.1108807303 -36.7627908654 -59.3889184719 -26.8307601801 -31.5313059635  -6.6267092903 -47.1711348406 -63.3172366024  -8.6064387274 -14.9969803396 -35.9726787133 -36.8102852033 -37.6082915182 -37.6734231556 -33.0585404446 -34.5521456026 -52.5278235884 -34.2166900616 -26.2724605543 -52.8225103723  -7.3317466741 -26.5756861109 -30.0169299845 -14.1652536669 -18.5479113008 -65.7198886489 -21.3491074684 -46.9837283379 
		 26 	-50.3108183557 -51.4346976207 -23.3501257442 -16.9123615036 -34.5345137465 -37.8220579232 -48.8419800076 -28.4416096447 -21.8081947236 -33.9843141437 -40.8644690905 -30.7343212735 -26.9197515438  15.7520943181 -44.4819457755 -51.7745379775 -39.9604352993 -34.1435191269 -42.0660970166 -38.4262560443 -23.1740524382   0.1992826561 -50.3577480331 -24.3527961882 -50.6242439133 -34.3831350813 -33.4814219791 -39.4945819126 -22.1709426678 -42.1566198110 -47.0361587081 -45.2688561927  -9.9256500340 -48.7872225844 -40.5284259164 -28.8935502134 -39.4828700325 -40.1444416856 -30.1917105644 -36.9448237372 -33.5337996185 -42.3554351598 -35.3673276814 -57.3698893773 -32.5696639038 -38.9455985607 -50.6959968472 -41.8501590979 -42.2685283241 -29.3067957694 -34.6884089713 -25.6729824039 -31.1679597871 -41.6822216381 -40.2289684939 -29.8806160157 -24.0120315591 -10.7861506110 -25.3033126849 -31.7925217394 -20.6623268488 -31.2983881239   2.3075946655  -6.8864959928 -19.9566062621 -31.8724191544 -27.3008002637  -7.5562889256 -16.8065096198 -45.0970619951 -22.8275615864 -52.7997570977 -35.1618529808 -12.7659566090 -22.0308818416 -25.5741029617 -32.2589995598 -36.9450959427 -17.8640159058 -21.7271840475 -43.1693440851 -17.2033994109 -11.2347742823 -37.5122602195 -61.3715297232 -15.6053243939 -34.9943948734 -39.2616620977 -42.1489193354 -25.3953215174 -20.2760333931 -41.6678504622 -38.6399423004 -48.9792048051 -19.3233476906 -18.7530569764 -48.9201517969 -27.7600497623 -28.4313905472 -40.3232032392 
		 27 	-39.6082522981 -25.2402656889 -36.5279599649 -38.2657789881 -36.0297197339 -16.0897262306 -33.3954077568 -28.5922995517 -18.9631166097 -28.6002818617 -46.1951091379 -44.3941835373 -39.9893492710 -56.6868415872 -44.2347359382 -31.6616702732 -41.5461713494 -21.4180794383  -9.6592390541 -39.5926110665 -30.7454184304 -31.1035126875 -42.5226519543 -33.1152418612 -25.1610521923 -42.6615749407 -14.3696962678 -28.7715796694 -51.9696347514 -34.2631405550 -15.0785571715 -52.4659968905 -29.4550878570 -54.9510949758 -26.9574441652 -25.0640098158 -23.5715994554 -47.4507136345 -49.2230521781 -48.7095655381 -45.2960420783 -39.1406611050 -30.4736227976 -16.2950037679 -22.2819690815 -21.6575250817 -32.5689939883 -45.6456673410 -37.7783887576 -19.1420711359 -19.9530968485 -20.7059986457 -40.0703690319  -7.8732892624 -18.5802493077 -18.1639992005 -10.9811707732 -43.4134035710 -26.3500840662 -71.6862854721 -12.3754565730 -20.0283876117 -57.0120239854 -44.3663536332 -50.8345023042 -71.3181616452 -44.6874879386 -16.0631760598 -11.4764291940 -42.9854027866 -29.9620481899 -18.2501318315 -29.4183222704 -36.2892228470 -34.8857877066  -8.4144716956 -59.4359677475 -55.0938033013 -29.5728620455 -61.2729283985  -1.4528543176 -42.7301234494 -33.4586892902 -55.1967624596  -3.1053734167 -32.7651512338 -46.4796661511 -30.3245067113 -47.2952582465 -14.1072820069 -48.2306392668 -26.4727138137 -14.2688053277 -28.2367817785 -39.1540007810 -35.5938171081 -67.5785678591 -45.3077512568 -29.1717527824 -41.0099955886 
		 28 	-22.7957724197 -34.4623032057 -35.9716108602 -10.4469102327  -9.8880813953  14.2415146445 -30.0517685619 -51.1550123254 -15.8980243062 -40.6663222099 -40.2813245892 -24.0091248766 -27.1561721097 -38.1626011668 -55.9564573907  -2.6020142988 -48.5988900493 -41.7468325655 -25.0320171625 -32.8648431745 -37.6323516526 -29.5516658954 -53.8866725514 -37.6261294348 -17.9934637960 -43.9280553051 -29.4492807608 -32.8739862868 -27.5310740582 -28.7239986816 -27.6870462880 -35.3698598035 -20.1740916533 -11.8600667409 -28.9245544886 -25.8461709883 -35.5758768810 -49.7210595196 -49.3228469615 -50.8100837056 -36.4028147658 -43.6767589916 -28.3917559531 -34.8429172098 -29.1417401023 -15.7102689932 -46.1391242610 -14.5243951547 -64.9370029138 -16.4018067893 -43.0624101318 -30.3671398379 -54.9359206885 -17.6526877835 -40.5765642578 -43.7784341065 -31.2884255972 -33.0358702302 -46.6563393237 -27.4472521745 -18.5873444802 -38.8348673353 -50.9033445334 -30.0120544736 -34.7293800460 -20.7510793347 -40.2914451555 -31.3475637665 -52.3551074785 -61.9478888985 -44.4382323808 -22.7178973375 -17.6499657882  -8.5175132101 -31.5673983113   3.3437760564 -47.9170894011 -43.2328534509 -30.6615358335 -11.0310603179 -28.9284019329 -53.7807986963 -28.2703585984 -37.2432674453 -42.7827652393 -38.8293798862 -31.5132814079 -39.1696564711 -17.9626405918 -26.5319755996 -47.4003496682 -48.0657112415 -27.9428008223 -27.0791081822 -52.0885668079 -24.0864148324 -25.7579895310 -48.1697882597 -55.0465634412 -37.6887987625 
		 29 	-41.8954364758 -50.3172082434 -48.8968817314 -37.5528971448 -25.4646954655 -31.4532665653 -19.2241821251 -19.0401103924 -37.4136850093 -27.9608525724 -55.0651916348 -25.1592250372 -27.0693931835 -24.0289767429 -48.0978449470 -23.7838339322  -2.3515054698 -15.7499914987 -32.7999488067 -58.5550004528 -43.0118970314 -47.2704435518 -44.5501530002 -40.0794231895 -12.4755079335 -16.2263150020 -47.5628399004 -26.2859840251 -45.4218477445 -55.0347749658 -23.6720332665  -9.2113668411 -12.2790422640 -45.8350992082 -41.9622108941 -17.2015975879 -24.1717550719 -32.6688488171 -31.7370533051 -34.1909493906 -65.6196594395 -48.3001613423 -32.8609328045 -35.8806138923 -23.2153870766 -24.2685247274 -24.1034086303 -26.9462067984 -50.6531713357 -33.2485936716 -52.8547449245 -43.1587982890 -34.1353876743 -27.6171167699  -8.0946307136 -32.9498263510 -39.5692059448 -61.7445944452 -37.9043487248 -13.9083828937 -52.1208865693 -38.4213932002 -33.0186159053 -15.4667607741   7.0413228840 -28.0393933233 -11.2096543485 -34.1822733995 -16.0421093689 -31.6806401344 -31.0557155579 -42.2956476761 -46.2225158624 -42.4984230862 -31.3429543275 -37.5894468618 -37.1990959810 -26.0237668650  -6.1850104204 -62.8576172667 -43.7516726516 -42.2761203549 -10.1866357005 -43.6553988827 -25.4337472119 -53.0758002287 -44.1794994732  -1.7922918165 -34.6906601584 -45.6866840886 -42.1135296273 -30.1707500792 -42.8787589451 -35.7527321217 -39.1922550305   2.5847071735 -38.2755269542 -31.3587960067 -44.5456894732 -42.3984080993 
		 30 	-31.9778436866 -51.3419843221 -58.9996284326 -69.7155970770 -21.3842191249 -48.4528232191 -56.8816096140 -38.2413563407  -4.5579140179 -40.7990162814 -27.3186528665 -44.5325454964 -43.2528019968 -41.3223638300 -34.3966167666 -29.6707381160 -54.0217347740 -25.6488957030 -11.9347086462 -19.7241603761 -18.2822018383 -34.4981979846 -24.3751061731 -37.4618745158 -31.4948131740 -15.5182327493 -40.5651688636 -35.5799173010 -27.6040165540  -4.5176885013 -47.8346941401 -51.9919352747 -69.1077877118 -28.7572094866 -38.9492414335 -22.6759171405 -27.0971114524 -27.8797620454 -48.1034038693  -8.7359416565 -51.7092387923 -42.1209534958 -22.7214295548 -46.0015463246 -39.9317759369 -26.2072015354 -40.5881695504 -45.0224089759 -27.6121659064 -35.1142971473 -43.5365859665 -25.9856079695 -22.1455774709 -20.9360496858 -31.2864429346 -34.9605150333  -8.8337942502 -28.9365068409 -40.9375972876 -21.4220407525 -35.6565355698 -25.2987582127 -38.7764540630 -26.5221226268 -58.0062520103 -19.8469444697 -31.3007270969 -25.5767711702 -33.0516128352 -36.6180475378  -6.4866860494 -20.5409704824 -15.1723635302 -24.8008402188 -52.2726614727 -34.2883613603 -52.5841310747 -23.5451677970 -54.7241382023 -60.9608094818  15.7466796023 -18.6798219314 -43.7305954969 -20.9035325495 -22.1825945361 -45.4553112217 -35.7936210147 -15.8573856045 -66.0859216511 -19.0163203447 -31.4377290641 -25.5291867523 -30.6968407371 -34.5305385284 -24.9591365902 -41.3782298563 -27.7647418050 -27.3638180785 -51.8164535726 -45.2678992797 
		 31 	-34.3620528908 -45.9012225196 -39.7754233642 -44.2620566441 -47.9835926934 -42.6204405545 -36.1093554630 -28.7410445299 -28.1152749095 -27.2359072328 -19.7204823530 -23.0707968921 -32.8602907526 -41.7387657944 -42.8293207697 -35.7491966721 -39.6082645686 -32.4429985979 -50.4341382563 -37.3192827141 -53.6962927435 -48.0641305427 -35.6419317833 -36.4609187963 -24.1547281074  -4.9289623377 -25.9819307642 -29.9703438333 -21.1361234522 -24.5916699740 -48.2567983015 -21.7201182963 -16.9842681644 -38.9480863897 -22.4227571174 -15.8280060240 -44.6611398721 -38.0293593191 -32.6388474345 -33.2695024871 -30.7549334850 -34.7519173139 -45.5608006260  -8.9599165602 -23.7490614881 -25.9393391325 -22.0670004137  -9.0315584817   3.6556390830 -25.0887761629 -32.7249755831 -39.3488370913 -31.3309496086 -31.5966395109 -36.4688521556 -32.7917790237 -39.4051876731 -32.2856306226 -46.3621955047 -24.0640836182 -19.3911151322 -27.1710620838 -20.1564290978 -10.9374931666 -46.3600063150 -45.8284522904 -36.7444273347 -41.4380844022 -39.6239053853 -33.8341054109 -37.5703804050 -35.9753646579 -23.9034243811 -36.9059207447 -25.9376444446 -54.1788373143 -44.0917161549 -14.2581455783 -30.2634298115 -35.2169295040 -33.3710690039 -30.4103861982 -35.2480181227 -26.5796623221 -36.9892411142 -26.3412544754 -62.0845494211   1.3233726777 -34.8349690793 -36.4348782244 -35.7918708833 -47.6493198039 -36.0690737220 -58.4346271873 -32.3948452162  -6.1295839604 -25.1056902947 -34.1896778498 -19.4560599079 -61.0710456703 
		 32 	-29.8942072494 -19.8826462471 -48.2023885191 -17.3447504374 -28.2626380454 -14.1133045970 -40.1362775247 -53.6159260391 -18.9588159811 -16.9070333799 -32.4169139493 -22.9277281600 -16.6677803748 -49.6248527097 -42.7106299416 -49.5187535391 -36.1929387251 -39.8116855442 -33.3465893551 -50.1520099950 -47.7347673650 -12.1825871763 -23.4569559973 -31.8403301458 -58.8143172826 -40.7509968473 -15.7701217084 -33.3768176899 -36.3456686158 -53.9447466627  -9.6749486599 -44.2051946086 -33.9324173097 -37.0675232182 -22.8300553520 -18.6835678062 -33.8542880340 -24.4852190884 -30.0378812861 -26.1981345732 -37.6149724594 -30.8505101336 -42.1245266629 -34.2006297923 -49.5142930456 -23.2119979248 -39.8201740686 -36.9177677440 -32.5861856670 -32.9551281112 -42.2169051224 -38.3039618273 -36.3890064347 -41.0653788605 -50.1536167842 -43.3874468547 -19.4226739333 -44.0428348276 -28.8767994210 -46.6545062411 -51.0729655268 -33.4325634926 -29.2551621670 -31.2140560615 -35.7781014421 -38.1600924785 -39.7261248843 -42.3426258969 -27.1365279994 -29.5723201237 -52.4456757858 -39.9741714796 -35.4547138957  -2.3826710259 -30.4619179285 -44.8157121650 -44.1575902800 -24.5323511122 -23.9917644118 -39.5020599545 -31.4233939545 -20.3184049063 -10.1975358234 -52.5078118707 -14.1761120914 -23.1695881017 -34.8732304256 -30.0544669021 -25.9198268574 -34.3458865156 -23.0636704016 -62.4035713125 -32.8855830391 -25.6553463002 -33.3873869742 -18.0599229145 -43.1613263879 -37.8910471066 -21.9048354397 -15.2789478168 
		 33 	-46.4458661993 -35.5787858717 -35.8220630330  -9.6149024153 -29.0384361041 -15.8411430659 -75.2833907456 -44.2851632488 -20.7539929786 -16.8800517614 -35.3684847702 -32.7392301770 -37.6830214544 -43.5174611431 -15.2497072483 -13.3003331520 -39.9677649685 -33.7154180454 -23.4332464847 -23.3653094960 -28.7459047423  -8.6523702365 -32.8766687556 -47.8242793004 -41.6213250324 -31.1272388758 -31.1957160661 -34.6379816414 -33.7323783748 -38.5881551918 -35.3984108171 -19.2741106235 -39.3204984962 -42.7962518632  -2.2569600686 -37.4152944756 -57.0583231590 -40.7660497913 -26.1449350530 -21.9284835647 -26.3557108741 -47.1502693201 -44.3488563781 -36.3392978553 -13.1535119269 -20.7067852715 -33.3342358643 -44.3268644047 -31.3037348564 -29.7340313751 -23.7437611776 -50.4308958411  -8.1323921736 -15.4461465032 -36.8153938038 -37.2625494479 -32.3326286975 -18.4032109625 -35.8761406188 -19.6605748819 -50.0192520416 -14.7872946312 -27.5721739733 -49.4046307172 -43.4742474365 -34.9663093755 -35.7654241591 -40.5289970907 -36.3654541626  -9.2339259923 -37.7066997131 -67.6140521933 -49.2966447810 -46.9135834702 -28.6925451033 -38.1878862044 -25.4902112819 -42.5859024427 -38.4782616787 -35.9217908403 -30.0489183447 -57.4979086673  -4.3742590204 -43.5031748580 -27.4602748202 -47.1902644819 -30.8677891440 -29.8033344417 -34.4313877548 -16.7881343093 -30.2389665755 -35.1290025280 -44.6120394310  -2.8540048112 -25.8739926948   4.4630082235 -68.1726831276 -44.1683664793 -39.3566950273 -39.3740703590 
		 34 	-10.8483140042 -47.6503433570 -42.0208813520 -25.5967785862 -61.6838753999 -67.7623264200 -29.6375370565 -23.8419417501 -26.2701741926 -23.1261096903 -44.7993990343 -29.6034749365 -47.8614935266 -27.5265475045 -13.4625477008 -58.5593842936  -0.3039491396 -29.5011270131 -32.2789323842 -13.8315986401 -32.5555473967  -2.2659842026 -37.4424226284 -45.0697438198 -32.0228819582 -62.4920407964 -26.5397634198 -42.8889972186 -20.5285529908 -21.7539968454 -41.2306194322  -9.8909615947 -22.3910425338 -64.9091894761 -27.9581939531 -74.5362911369 -33.5450205468 -67.1921555566 -45.1785127109 -38.0459363921 -55.6263852043 -30.6544614468 -22.5154243808 -48.6199459821 -14.4855761038 -25.4360801885 -47.4116796214 -22.1178373894 -51.6678406046  -3.5933231954 -26.1652644918 -36.1994401493  11.5497124830  -8.4291921801 -38.0429575067 -47.8314055534 -30.4742938389 -37.7991184054 -15.6530386416 -30.7016560551  -6.7909051076 -23.6782615395 -27.8814297981 -34.0780327179 -42.3520284030 -32.9351237933 -44.4609969483 -49.2651339470 -32.3761965227 -27.3187477299 -50.1111849038 -15.7736034510 -49.3305132840 -15.7480811748 -30.9270398122 -35.6545572415 -55.8186387801  -9.6092635254 -34.5783409713 -38.8331884939 -45.4453416286 -25.5304718010 -43.7121929843 -44.4418317120 -24.8118038376 -41.8095159929 -17.8551836483 -27.0379929611 -23.3995704987 -40.8213694551 -26.6690858200 -44.3406946132 -52.3385950990 -34.4151505912 -34.1254010370 -46.6751871887 -31.7728533764 -20.6205312352 -46.5461459305 -39.4117640764 
		 35 	-30.4645223428 -19.2152429986 -34.3541483417 -35.6844502122 -39.7925166952 -63.8282113579 -47.6142756974 -15.3820957705 -20.9339371942  -9.7181407092 -24.8406476352 -27.6132010776 -55.3609187453 -64.4833314664 -46.2345756495 -22.7375482035 -30.0330822551 -37.3880304356 -47.2267754155 -15.1823466650 -40.4708074627 -62.2106286942 -33.8775371545 -32.6975904800 -19.3706597396 -28.0511874595 -29.7537754485 -19.1794645697 -30.2267327596 -13.5288168409 -39.6435221271 -28.5309216463 -59.9096871251 -21.0002234954 -20.4124901139 -41.5538615616 -35.8819553699 -25.2530565847 -31.3553288586 -45.7697080833 -37.5780406884 -26.1453679105 -10.7879271471 -29.7178213969 -26.1774358123 -39.8275732335 -39.4854739551 -33.0770637540 -32.1557757723 -45.3910357652 -53.5514936362 -65.5319573259 -22.3039880557 -31.4382755912 -26.6303939691 -23.5210306240 -31.3631749206 -63.5722298007 -53.0627635183 -16.4332111845 -56.1555058112  -4.2323865386  -9.2113668411 -24.4056636633 -36.8640510565 -25.4975288985 -53.8554465026 -47.3721822061 -46.3734141524 -51.2865481851 -19.5586337930 -42.6532634976  -4.6770630669 -42.9644668508 -17.5396973241 -23.7821065850 -23.4130845232 -39.7429715181 -37.3846118588 -29.4101881399 -35.0484593370 -20.2371561343 -10.9867571633 -31.8091983761 -51.8834769630  -9.5334774482 -15.6368112689 -37.3948300295 -16.7788841943 -26.5012081337 -40.4430256098 -55.8283902104  -6.6060192464 -73.0333723743 -39.0651187588 -35.5383916960 -21.9998718409 -34.8585701179 -28.1806906380 -28.7060433742 
		 36 	-32.9224484671 -72.3676331860 -31.0885607952 -14.6486999079 -11.3787386961 -36.6673531648 -33.1646338704 -25.8474977261 -54.2477429442  -9.3378978149 -37.6207362832   1.4176305038 -27.6315773707 -13.7202585018 -28.7343266196 -17.2516997759 -46.3674361557 -20.6258171301 -43.8816678699 -56.0518088330 -38.1454191729 -10.3040388934 -14.7258836622 -24.1420547736 -29.7127878543 -27.5267173635 -21.2993733699 -16.3047978230 -64.3805942786 -32.5033962107 -51.8837177236 -30.3211437123 -53.6639719655 -17.2871215960 -56.8476793950 -23.0985010439 -53.1439642005 -36.0041363564 -29.7331681221 -15.6368700181 -33.1327464989 -45.7933098965 -38.5304284135 -51.1376539554 -31.2912756429 -27.1725631540 -34.3442204532 -45.1210261587 -28.0452877829 -44.5622499068 -42.3185838189 -32.3845386037 -26.2102686399 -61.2319134137 -40.0474693494 -49.4990335131 -26.3903499296 -15.6895337590 -41.0375110815 -25.3910538953 -38.8154611263 -43.6319661252 -31.8074582757 -37.9651723354 -38.0179362075 -56.0014490498  -6.4607184713 -22.7615437611 -28.8405027074 -29.8609235297 -32.8463948935 -44.3492862709 -51.5678898243 -50.2135796737 -30.5775651802 -56.3484122906 -28.5288743809 -33.0975561878  -6.4916120774 -26.9628760503 -35.3155684195 -37.7644553034 -32.2274773973 -21.8509610415 -11.7171458164 -31.3503940258  -7.8514412771 -31.7443757687 -15.9935706106 -44.3185035121  -0.6419434019 -18.5314236089 -26.4823488466 -56.9397839433 -33.8948697053 -40.8667188158 -33.2043230254 -25.4097524076 -61.1252092179 -13.7030096586 
		 37 	-45.2711935991 -22.8716840244 -20.3762855794 -23.8975687313 -15.9645928392 -29.0933360535 -28.2307737248 -40.8973688815 -45.3795156347 -42.4059083782 -30.7096995039 -45.4746186833 -36.4239161932 -60.8528186300 -26.2800743122 -44.4454893127 -51.6820717522 -21.6429490154 -43.3813276479 -14.5550993727 -41.6179034251 -33.5068791426 -42.3031137052 -30.2473787424 -31.1158153197 -11.8395225124 -34.8604555147 -26.1300277413 -55.2313528776 -21.6614316294 -28.4112400016 -50.0968681894  -7.3812231950 -44.3139498482 -35.8835534021 -22.3608365983 -29.3230846134 -27.5119777904 -15.6222137759 -34.6195485310 -38.1959351526 -46.3123872622 -39.3551650795  -8.4164170623 -10.6367093487   2.6184377607 -24.6516934440 -26.5513625968 -26.6384932794 -59.5882247468 -52.1980271315 -33.6873629839 -30.7297303243 -39.9353368086 -45.7315729032 -25.9343598272 -22.8354415094 -29.2204460368 -41.1758910272 -31.0112607656   4.6420532900  -0.3992749866 -41.6714688638 -23.5053270368 -40.7302671946 -25.7588587917 -24.9824690460 -52.0582970011 -16.2010194196 -51.8433283747  -3.8393296703 -42.7909312249 -48.6612654705 -29.4351161360 -60.2350773149 -31.3812381027 -56.1450391552 -21.5021084568 -21.3295450813 -24.8680204602 -48.2807492345 -41.5775278005 -22.9106233229 -29.3871434636 -31.3194538964 -35.2758491163 -39.9338861849 -36.8567713919 -26.1538202110 -29.0750781338 -20.5764474113 -33.5817458331 -49.0212812130 -13.5771872043  -8.5103322283 -69.8190835002 -38.1690600198 -38.8684912718  -3.6486847609 -30.8609814761 
		 38 	-20.6558226736 -46.2221717543 -51.6638202293 -35.5446739984 -29.2933341239 -16.0304461435 -31.7257388632 -31.5280616056 -44.8413974643 -21.7572688157 -18.9909928731 -34.3344734957 -27.7848094075 -34.4940049679 -55.8835941766 -14.1757859981 -28.9919460500 -37.3226681862 -25.4777985302 -40.7942749689 -26.1665112222 -38.1743417608 -21.8354586801 -36.6241098639 -29.4773858120 -16.8943114021 -39.8857639841 -29.0845357801 -30.6190193558 -34.1700581947  -1.0930417971 -40.2648071782 -54.9314770169 -53.0816096955 -24.7661382351 -53.2228856254 -19.9545002248 -17.3200736755 -24.8572424752 -29.5172991179 -20.3429349310 -49.3438221681 -10.3712085236 -24.8019251325 -39.9096476314 -16.6896822360 -10.2903004311 -40.0547473656 -13.1670028330 -45.8819324853 -50.0159265822 -24.3549190375 -41.9676932447 -39.3498119127 -21.5485589473 -27.9354486559 -22.0456115219 -11.6234005787  -6.5400256792 -40.2335015254 -47.3655959835 -46.3304139372 -30.2381201185 -23.4820738417 -29.7622113233 -44.7721168388 -34.9166239702 -54.1487792370 -40.2054218088 -35.8041277925 -27.7053452244 -28.1705745720 -39.7050490260 -22.9230386716 -14.9684236413 -28.6349595552  -8.0811335654 -45.5487129166 -39.5273678656 -27.1856283017 -35.0556630917  -8.5500345541  10.9214606865 -29.4860721671 -30.5550993920 -15.9345702614 -52.5625731641 -38.4257829109 -28.3437462440 -56.6022672255 -44.6517291931 -46.6781493077 -46.1207324924 -23.9205230881 -53.6262202001 -26.9837467387 -45.5711500448 -56.5700767486 -17.6138772810 -32.1007973291 
		 39 	-26.7384514507 -23.0431849800 -23.4547048729 -36.8676047769 -27.0091797967 -34.9620466826 -31.1493683560 -22.8251764339 -40.8610660520 -13.8791730834 -41.5249202255 -25.9746302387 -30.1433906430 -59.7283053667 -30.9727919643 -54.5484167553  -7.7544679586 -35.5496469723 -33.6833978794 -53.5595643914 -38.4064105953 -43.5131292824 -20.0723906155 -45.3115218092 -18.7798909231 -25.8250504874 -30.2271174691 -23.4337828833  -5.9852721247 -36.6822580176   8.2535443448 -54.7196886352 -12.1498572751 -12.3491460834 -30.8774543255  -2.8674509438 -45.1242134181 -43.4197001687 -19.8253571655 -30.0628614180 -20.3726845476 -22.4833997437 -23.3230916637 -52.4777757917 -67.6266613314 -32.1154798333 -14.2126925798 -42.5156581122 -42.0595948106 -34.0153860955 -34.6380636837 -38.4768853307 -55.1938143948 -25.9029667875 -37.4289364206 -28.2270803623 -18.2777044825 -47.3713070192 -17.9391169906 -15.7075626178 -29.0644981787 -18.3852266232 -57.3173289061 -66.5767250849 -10.0674130874 -28.7535083162 -30.4283729217 -25.2798101035 -27.2593686703 -37.4340666407 -37.9949385165 -34.7804981261 -27.4329530660 -39.9472465821 -29.3927649862  -9.8191942297 -39.2031027939 -35.0538535044 -27.6981253980 -10.4771916140   2.2372622436 -12.7850963115 -19.2236815186 -11.4996145319 -35.5289215221 -10.6633435564 -44.2812859144 -23.3685215519 -23.8420615650 -48.7714511990 -33.3450797040 -31.0375124719 -32.7062535153 -15.2159493486 -40.2223625863 -53.1306428094 -36.3587880671  -7.2600269741 -21.4716834206 -39.1278515518 
		 40 	-22.8619556400 -56.0466775601 -36.8098646894 -27.9664222840 -57.3632751685 -49.7655256600 -47.3189843734 -36.5783166407 -42.6748440987 -50.6754346602 -27.6895495261 -38.4924777354 -43.5087615458 -24.3414193052 -19.7290073432 -38.0683872987 -21.4358123528 -71.9433553418 -25.3030606627 -51.1931062008 -12.5998173695 -16.7564276796 -33.0254883554 -22.3096510291 -39.2535360615 -47.5208184124  -1.9674356857 -37.8683408145 -23.3885800057 -38.4983109831 -41.0827142625 -36.1906943399 -11.6186223375 -35.3233368203 -40.8602095143 -54.0297689224 -35.9450919638 -43.7255952998 -28.4421738719 -22.0666819506 -32.1479996775 -27.1106931510 -12.8821819185 -56.8242481220 -42.6789302716 -40.8786924590  -5.4644648447 -25.8264303695 -52.2966881828 -28.9810580402 -38.4489864679 -30.5969084270 -31.6755393232  -5.3608231767 -61.5242052691 -20.9214889611 -34.6639443101 -79.1193434996 -21.3529419302 -38.5569028399 -15.1269925800 -21.2345705016 -18.6092325401 -30.4168472170 -36.2750985050 -19.5041356620 -59.5095606865 -27.1695972200 -51.3216498539 -54.9334891251 -24.5136800574 -29.8138367657 -45.0580307537 -63.4662509935 -45.0585767855 -16.2662631227 -22.9755025141 -21.2269235521  -5.0881872211 -53.6047763334 -46.7616470931 -29.1204731961 -20.6574871271 -28.4633957615 -35.0126745826 -59.8269383642 -26.1578843175 -28.8553410739 -14.7673932107 -50.9951078744 -18.8991261874 -24.4645681771 -47.9179595903 -56.5052747639 -26.6005540600 -52.8876904980 -52.7484256279 -64.0760081687 -44.9796319913 -32.8973559751 
		 41 	-36.6856668297   5.8226078670 -45.5947744946 -39.2089664528 -42.2001622300 -28.9331573277 -26.0923171764 -21.4123125154 -16.0800359294 -45.5823234351 -41.7222857557 -20.2478759426 -23.5360149170 -58.6332631405 -33.9738862614 -55.3828159918 -15.0084863272 -47.5329685216 -28.9462596160 -38.1150462662 -49.5619466339 -52.2281757858 -33.0966015759 -25.9599683233 -35.5772397848 -38.7059684184   0.4235839231 -30.5283230175 -22.2765632235 -28.0856947243  -2.6957673914 -39.5860735728 -66.8654754658 -24.5060320698 -21.8462456300 -58.7252693221 -36.8029855983 -52.9420821469 -47.3849715061 -48.3736040402 -45.6054934924 -54.6878065022 -27.2754527201 -25.4485324325 -36.7864310648 -29.6279624645 -26.9277869955  -3.0870415182 -40.8098920973 -36.5325852069 -21.9804597505 -41.2899228903 -22.3669094052 -40.3116484046 -38.5171926263 -34.7800142273 -54.3797513906 -28.2414658584 -70.9336678681 -45.2793689274 -34.6820647043 -38.0014608284 -18.8169232973 -44.8542044706  -9.2448534344  -7.8077261580 -25.9347936525 -18.0077474653 -31.8576057347 -49.0182268379 -37.9366146523 -34.7062727870 -23.2863325424 -48.9555732298 -17.4799242146  -9.1154674574 -31.6207432534 -55.0562278224  -1.9717091813 -44.4783144056 -20.4099126570  -3.5803660335 -27.7561060928 -33.6552257739 -32.3785880679  -8.4482808198 -17.2305879903 -54.7646443952 -34.9118680758 -53.8020979899  -6.6681033285 -40.6328998743 -38.8252787042 -26.9944284942 -28.7064188373 -19.3542987291 -28.3532348617 -31.8780195206  -1.0538485217 -56.6306581271 
		 42 	-23.7293031885 -32.9421084077  -3.2639445224 -39.8593949638 -43.4602717072 -35.5275390623 -43.4383503702 -19.3791524651 -40.1985472795 -25.0778864404 -27.1863062890 -15.8980243062 -15.0898499220 -36.7795652533 -22.5346165515 -32.8315968184 -33.3538730287 -29.9772543957 -36.8548887982 -45.1442356577  -2.1512630335 -35.5498294838  -9.9422643891 -35.9218456361 -55.5216097271 -21.0704076896 -64.9469515722 -24.4900561457 -11.9289603949 -35.1457992653 -27.8463738678 -27.5033110318 -40.5275889932 -50.9793217169 -59.5793646191 -46.4334347922 -21.7135746039 -37.7742862200  -3.2494609666 -36.1728745926 -40.2752458790 -12.2028399947 -39.8564125297  -9.8025329909 -51.1165927910 -32.6118965240 -65.5998272405 -50.7614442474 -29.2066416451 -35.2744102772 -20.3061097538 -40.8301893170 -31.8233441761 -29.6076858830 -45.9376976750 -32.9207359711 -18.6660732242 -36.0695898969 -33.5516217980 -34.5570921614 -20.6336388542 -43.1880076063 -35.0104669353 -11.3012018609 -20.6855870706 -30.5343086716  -2.0864823853 -43.8307574623 -52.7815615051 -18.6304322677 -23.5153353514 -23.0513555682 -43.2801829503 -21.9404465130 -29.2830760043 -36.8351642320 -34.1787501734 -50.7105286795 -14.2497621497 -47.3853611695  -9.8307319790 -52.9059527273 -52.3283504067 -31.6172076844 -62.8523749938 -42.4289865790 -26.3725622680 -14.6194182637 -39.1678531096 -41.3500392412 -54.3078275428 -31.3394934581 -35.5213188793 -18.1469863719 -28.3226093067 -29.8997291167 -26.0822299817 -24.9967952104 -34.0530208989 -21.3256678517 
		 43 	-32.8774659694 -42.3209683570 -39.0157998335 -31.6001580393 -40.8668399783 -24.5706781522 -40.7075251150 -36.8737548159 -25.7684860928 -18.5307071898 -36.5562109099 -47.2721289346 -39.1735626405 -34.1671568845 -26.4951028477 -37.6092932979 -25.8753043782 -45.5866755965 -22.9886564326 -11.0125840962 -30.9033371938 -26.4623667523 -11.9759826746 -32.5988021825 -20.1920128994 -41.4168591089 -25.6936901464 -21.8164785914 -22.4427620277 -38.9324047766 -38.9463846533 -34.3072189723 -35.9471633259 -28.6711264298 -49.8810648318 -38.3274314498 -31.4770164075 -32.4539432152 -32.0860418361 -24.6168294386 -34.5341969093 -32.0142696256 -48.1691968003  -9.1942825672 -30.5628602662   5.0085968954 -30.4223925441 -35.1517036396 -61.3856625754 -49.0313882690 -22.0197083918 -27.5704395845 -14.1887754668 -25.0809131106 -12.8595665113 -49.1478467739 -21.5373933096 -27.6991254678 -28.7827676436 -51.6993337468 -28.9494105221 -44.8445616102 -22.7930152876 -29.3834569341 -33.2159234664 -42.5145310090 -36.4430516913 -35.9242010895 -36.6846058686 -38.5390343831 -44.1272036150 -18.4521910377 -29.4005024827 -18.2707277370 -28.2791216475 -33.6241612360 -20.8201088847 -32.4399341322 -41.9512798213 -54.6092686226 -51.3314647315 -22.5023068621 -40.4525144217 -16.4252095107 -27.1073276977 -48.4135430935 -54.3330296067 -29.1652118635 -52.1562133792 -34.4986227290 -27.7380044501 -33.6536389033 -48.7144572870 -38.5813825376 -17.7850512018 -51.3679391048 -55.0909562203 -17.5705851086 -29.7785831126 -22.7734763344 
		 44 	-66.9474609063 -15.7811355083 -26.8296841551 -32.6752926598 -39.0390046089 -35.7305221191  -2.9278186506 -45.7051159060 -52.9616585846 -43.2471409419 -36.7306776202 -51.3113474408 -56.1390200643 -31.5832890399 -24.8418072603 -46.9740970673 -22.9539589059 -14.7122373877 -26.7194931720 -31.1730394608 -29.7218942068 -41.3601624150 -46.8820369630 -73.4263910176 -20.1340720960 -30.2452567944 -19.9825940090 -29.6430365957 -16.9684483266 -26.7330603343 -45.8690538777 -24.9547195097 -29.4246559975 -16.9185007602 -42.9706288560 -37.0634509455 -52.3535620946  -6.3114096636  -6.6693424935 -16.5105326660 -12.5839094910   7.7242276445 -22.9988863911 -40.0793234539 -33.4203746956 -47.1259423650 -39.0620872989 -32.3061702427  -7.9005453697 -28.1342327035 -14.9239816917 -38.7896678623 -11.0670673118 -24.2105461286 -23.0130998052  -5.6461610451 -28.5366645913 -13.9438255315 -31.6247956772 -11.8065878478 -35.7278653590 -25.3849322071 -41.3274007157 -19.0980574813 -44.5399352227 -28.0197897242 -50.4471121334  -0.8328832692 -24.6189935487 -24.1490433174 -28.9546786546 -23.6469803901 -32.3557323861 -45.9012225196 -30.8390743385 -52.7480044223 -12.1914019618 -37.0774944919 -29.3877887098 -45.7926110593 -31.8241792026 -38.8563481443 -18.2909669957 -42.0858236867 -32.2921878804 -27.4544544786  -6.3478158011 -24.9403307042 -33.2909583480   4.9426638462 -10.9620558556 -17.1840476425 -32.5974209169 -27.8682608398 -37.8901042443 -30.2487211989 -38.4253198189 -28.9145245813 -42.9654054435 -39.3847655864 
		 45 	-30.6430409680 -45.0117776950 -21.9007042662 -43.4501578712 -46.0708001225  -9.2013312607 -18.9887816620 -56.0683084030 -24.5965151152 -55.3224382805 -24.5258189611 -41.0901175455 -24.5088392010 -26.7072882236 -31.8334627452 -51.8297500384  -6.5090503043 -25.0578355291 -27.3770306958 -33.8344601504 -22.1503434783 -34.5662768571 -36.6000696367 -36.6834587915 -35.0738360192 -48.9900834471 -34.1639235788 -38.0437438695 -26.6699355309 -20.3179207071 -46.4473324930 -16.8200001004 -26.7159080749 -35.3233368203 -35.4523556269   0.1563692438 -40.8773951862 -25.9993932213 -48.7191987608 -18.1136675644 -32.7640704027 -46.8949529505 -31.3630391128 -54.6144849386 -31.0008229431 -21.5629525220 -26.7106897526 -15.8635568868 -33.4893077853 -37.0637880901 -41.8841457531 -32.2934303102 -24.4661315921 -35.8704418344 -17.4265337846 -14.5316761622 -43.1022504287 -24.3555630044 -53.5161379430 -37.9944437754  -4.3362244590 -27.3328473047 -50.5368323895 -29.6621835096 -58.8647781309 -45.2020522540 -41.4504216725 -38.8179107872 -31.9718535545 -27.3904043443 -21.4282945311 -18.5611288663 -33.7633959075 -51.9028827883 -37.9491169964 -29.7316592569 -34.2136216696 -35.1268935349 -49.3794918255 -32.0713856017 -41.6891357486 -37.7632326506 -37.0036885510 -29.2843201133 -53.5091225645 -43.9184087025 -40.2201149849 -27.6667819793 -29.8981435909 -49.3716261666 -48.5043219847 -31.1799697532 -29.9074537666 -29.7163330027 -29.6600662448 -47.8049205550 -36.0876050613 -40.6203809008 -30.4233351241 -50.3227847593 
		 46 	-28.9964797925 -27.1341505649 -27.5190652956 -37.8924050379 -21.3365870272 -35.6210184390 -31.0721501014 -43.3108308427  -5.3140017317  -3.9120779291 -39.7223626685 -42.0437711613 -14.2994902689 -14.4051803601 -18.2740991661 -31.9657772667 -39.4359154595 -34.5187972619 -24.1056822977 -29.6916429288 -15.7599125287 -23.6893581389 -37.8040000239 -29.8938910453 -23.7574775448  -3.0601779559 -42.1927735068 -30.5503297276 -22.1933906851 -31.4523114234 -19.7434781278 -31.4086008346 -42.9780557085 -33.1591869781 -37.3133006043 -33.5676278002 -46.6288195421 -43.1737679496 -39.0435261708 -31.6830520612 -17.7572357653 -28.5843687003 -37.6329823691 -22.3980050935 -17.2811707836 -42.5782374205 -20.7249417405 -45.9648089327 -38.3260790507 -16.2307955795 -18.1693486225  -9.7450322297 -10.1121221606 -20.1010607385 -33.4597615913 -45.9744261031 -40.1346944575 -28.2936764658 -48.5876094816 -55.4879303940 -17.8703169778 -38.0487184205 -35.3540473999 -24.8075227481 -26.8710606732 -31.4820323587 -25.1225384065 -55.2765643623 -49.0554458979 -28.7593779238 -54.3783629185 -24.9527740694 -12.2220101325 -37.6321926950 -37.5788365740  -7.0822401118 -30.8700168174 -53.9197031476 -45.9630577721 -41.0503399980 -18.5016724732 -10.7308683316 -26.6657690593 -35.1477114285 -32.5444252515 -45.9184988161 -24.0329757126 -14.3221854011 -54.6869922044 -40.6442513800 -30.3905346841 -18.8693607237 -10.0789860759 -38.9916481296 -15.1232098092 -33.7209875531 -40.3157517171 -31.7959575820 -39.1465966672 -15.3042468041 
		 47 	-38.9484050371 -36.0384315521 -28.0281535165  -7.8317206817 -34.6590295154 -16.9349913500 -39.2782511198 -29.8753091615 -20.6444428118 -37.2155930704 -25.4461110053 -27.0903858616 -28.8553410739 -27.3571350489 -24.9149291416 -43.1416814903 -43.5130881579 -32.9135563495 -34.4109828509 -34.8515196453 -59.5649200559 -12.1227886373 -54.8419899659 -21.0070810553 -32.4052385968 -17.7745531152 -17.8817248486 -39.6154959639  -6.4003871245 -19.4523262586 -26.4502160898 -18.7846819767 -49.9271492636 -31.9238947283 -39.4743153884 -19.4658376763 -29.9462952560 -44.6184387208 -39.1445760233 -30.8928445144 -52.9770197799 -36.2341453523 -59.0014339000 -42.2798440229  -9.7189593669 -12.6468910937 -30.4436968736 -22.7136196075 -23.8238803958 -36.9730474322 -64.5279048443 -46.3144377854 -34.3883757860 -43.2248345170 -29.2403931159 -54.9310203908 -60.9629787719 -29.9788345545 -37.8495000682 -42.6864490079 -26.7593551418 -18.4271515490 -48.0062247151 -44.3560180430  -3.9746506175 -16.4697260522 -32.9619346531 -56.7677486833 -59.8382177559 -24.3820128538 -33.2184132861 -23.6977731672 -63.3342457811 -54.3216095422 -36.1926654936 -34.5481030576 -24.7102664741 -38.6499637982 -27.9451654758 -29.4631570676 -41.3337601611 -22.6526449198 -19.1693273975 -41.3248126696 -50.8340596328 -16.1736106549 -26.9217821161 -41.6124169608 -32.3713384930 -27.1581517386 -50.7684048827 -42.3096747348 -29.6188727529 -16.8001466099 -41.0312665462 -14.1191723354 -58.7361027201 -49.6903021232 -35.3227370859 -29.4848792918 
		 48 	-18.3496748342 -37.4726669041 -37.8244118108 -55.4385514372 -28.1336509787 -22.6102681398 -55.2109101660 -18.4464899044 -14.2814658193 -33.8118768929 -30.7979094525 -19.8483026205 -27.6765925810 -30.8649812290 -32.0586190295 -19.3220235029 -21.7633094891 -13.8675161980 -49.4399312629 -44.0195599776 -12.5033142273 -39.6856429473 -13.2172129587 -36.0047265654 -18.7707023607 -42.7275205632 -43.0271334939 -30.3819901375 -26.4605234307 -24.5725978265 -24.1056245530 -28.3654160243 -19.5018386454 -26.6253786622 -36.0946077412 -26.3424777276 -48.9047990315 -28.4958985961 -39.7549019843 -18.5715521880 -25.9858704166 -21.3595115623 -41.6449473161 -34.9673239650 -51.6696541290 -18.9646552040 -20.3869654027 -28.7941585417 -34.8281300485 -14.4465209113 -47.3655388835 -47.3319906991 -28.3149657688 -24.3395780990 -42.3547805915   6.6454862407 -55.9842353876 -57.7834092527  -7.1313067898 -28.1139423640 -48.1073093854 -22.5900202559 -37.1629337029 -39.3090421687 -52.5036247687 -62.4456299904 -11.7169662077 -37.8033321239 -39.2842409921 -67.3621752526 -36.6175684674 -55.8765324753 -30.9252391380 -23.3973331998 -34.9582472597 -30.4220154706 -12.4651633095 -42.9079203728 -41.2675946166 -28.0240235278 -37.6400930434 -17.1812019185 -16.0119340305 -30.6879915710 -38.5440393837 -10.0722857795 -47.3312929313 -32.3454176528 -38.6309596132 -54.2321265487 -41.2791249520 -39.5643395061 -23.6095756943  -8.1749838960 -34.9227621096 -19.4367300287 -42.2672619487 -15.2315230922 -27.6078695757 -54.2608187026 
		 49 	-53.4677189126 -37.0975552370  -2.3170427112 -33.9904124514 -46.5346763699 -29.0641634531 -31.3712254051 -44.1928080440 -28.9246191211 -17.8276206297 -11.7811527371 -36.3015698537 -43.8152032629 -51.4010876463 -37.4110091709 -37.2718486436 -19.8815714450  -6.4152601775 -31.5421625909 -30.9343010810 -34.2587888925 -23.0643559308 -26.1690305187 -43.2176470170 -14.4262987046 -30.4145508924 -30.5272080898 -55.9043317279 -36.1102790861 -29.0204017193 -29.0923448235  -6.4316096489 -22.6944132080  -3.4207017357 -32.0246533259 -22.0538495566 -62.2107238016 -26.5177830819 -39.0516174763 -22.5392432767 -48.1389823035 -48.4525224154 -12.0913885764 -32.9627612203 -44.3413165626 -22.3042760606 -23.0726412454  -3.6401335846 -31.9421513335 -34.9528764192 -21.6613406283 -31.3150795186 -42.7832309235 -39.5609423963 -35.3822869460 -45.3098722427 -51.5624526161   8.5152582657 -27.2986417558 -23.0769379834 -26.6169507172 -27.4269145761 -32.6504823311 -33.1479997844 -25.5220027206 -50.2335373206 -47.7538835547 -37.5689798425 -39.1444412342 -41.8822327274 -16.6787232033 -35.5556630332 -30.6368467331 -32.0806149613 -35.2154636363 -45.2738859130 -29.2789816489 -52.3215766136 -49.4812801179 -30.3877605164 -14.5856637942  -4.7217322070 -36.5454705058 -39.2318768255 -26.8817502621 -37.7172528338 -55.6592081021 -42.8975412617  -4.7918909949 -32.3639047257 -30.6140947317 -22.2576100053 -37.7590260128 -17.9481930242 -14.9553674367 -18.1642786409 -20.0761619413 -42.5087846554 -33.9522979414 -33.4953939512 
		 50 	-17.9921261050 -35.9933573316 -33.7332433213 -21.9738416596 -23.9851017478 -35.4226740656 -43.1613263879 -19.4447285046 -23.7266621607 -50.9887092408 -14.0294873831 -37.4509203980 -24.9433833217 -47.4436593317 -14.3817447442 -25.7521522231 -55.4660982303 -25.7421285799 -24.5360586663 -52.1528280124   1.3625388130   0.4402387307 -19.4209666989 -17.9900798921 -30.2784654412 -15.7459954664 -26.5810525500 -35.3747317717 -48.8049326412 -27.6765925810 -34.4449325724 -28.9074178373 -29.0488513785   0.9920014103 -45.9866463825 -49.4234930337 -46.6316722787 -34.4393721728 -33.4883708997 -14.6210630955 -14.3247006131 -52.1708042943 -42.8000302306 -46.3212185426 -42.8974552797 -25.9966203683 -34.4825469545 -26.3786978173 -48.1550964104 -32.1498070615 -29.3303089188 -33.0102909309 -58.0832997150 -44.8893568415 -18.2986221047 -23.4150911252 -24.8621283983 -51.9702773934 -28.3545810588 -21.9456146974 -25.6733731472 -32.0528919424 -20.0138936947 -18.9853116704 -17.3286291158 -37.1043254143 -52.9741202255 -53.5094274770  -9.9930851004 -42.9645749167 -21.2408968959 -44.7405951906 -52.7982285790 -21.3319617662 -55.3079076900 -14.1242820858 -21.9501550013 -34.2606091656 -72.5953403546 -21.4605840932 -29.6454431968 -36.2582632771 -19.4695251598 -36.6702237278 -18.8863062033 -12.3170649541 -35.0583531827 -37.3234170633 -29.8057038073 -47.7870493352  -5.3748265903 -29.6316918327 -33.1723474981 -29.6570548335 -25.9460246649 -38.1078973065 -31.7097627161 -29.1332951690 -56.7949270265 -23.1266396218];
		Beta_FEE=Beta_FEE[:,2:end];
	#Variables choice model
	  AT_FSP = 10;
	  TD_FSP = 10;
	  AT_PSP = 10;
	  TD_PSP = 10;
	  AT_PUP = 5;
	  TD_PUP = 10;
	  Origin =[	1	0
	2	1
	3	1
	4	0
	5	0
	6	1
	7	0
	8	0
	9	1
	10	0
	11	1
	12	1
	13	0
	14	0
	15	0
	16	1
	17	0
	18	0
	19	0
	20	1
	21	0
	22	0
	23	0
	24	0
	25	0
	26	0
	27	0
	28	0
	29	1
	30	1
	31	0
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	1
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Origin=Origin[:,2];
	  Age_veh =[	1	0
	2	0
	3	0
	4	1
	5	0
	6	0
	7	1
	8	0
	9	0
	10	0
	11	0
	12	0
	13	0
	14	1
	15	0
	16	1
	17	1
	18	0
	19	1
	20	1
	21	1
	22	1
	23	0
	24	0
	25	1
	26	0
	27	0
	28	0
	29	0
	30	0
	31	1
	32	0
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	1
	40	0
	41	1
	42	0
	43	1
	44	1
	45	0
	46	1
	47	0
	48	0
	49	0
	50	0];
	Age_veh=Age_veh[:,2];
	  Low_inc =[	1	1
	2	1
	3	1
	4	1
	5	1
	6	0
	7	1
	8	1
	9	1
	10	1
	11	1
	12	0
	13	1
	14	0
	15	1
	16	1
	17	1
	18	1
	19	1
	20	0
	21	0
	22	1
	23	1
	24	1
	25	0
	26	1
	27	1
	28	1
	29	1
	30	1
	31	1
	32	1
	33	0
	34	1
	35	0
	36	1
	37	0
	38	1
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	0
	49	0
	50	0];
	Low_inc=Low_inc[:,2];
	  Res =[	1	1
	2	1
	3	1
	4	0
	5	1
	6	1
	7	0
	8	0
	9	1
	10	1
	11	1
	12	1
	13	0
	14	0
	15	1
	16	1
	17	0
	18	0
	19	0
	20	1
	21	1
	22	0
	23	0
	24	0
	25	0
	26	1
	27	0
	28	0
	29	1
	30	1
	31	1
	32	1
	33	0
	34	0
	35	0
	36	0
	37	0
	38	0
	39	0
	40	1
	41	0
	42	1
	43	1
	44	0
	45	1
	46	1
	47	1
	48	1
	49	1
	50	1];
	Res=Res[:,2];
	###---PREPROCESSING
    q_parameter=ones(NUM_POINTS,N,R);
	#Calculate the part of the utility that does not depend on the endogenous variables (g)
	for n=1:N
	  	for r=1:R
	  		q_parameter[1,n,r] = Beta_AT[n,r] * AT_FSP +  Beta_TD * TD_FSP + Beta_Origin * Origin[n];
	  		q_parameter[2,n,r] = ASC_PSP + Beta_AT[n,r] * AT_PSP +  Beta_TD * TD_PSP;
			q_parameter[3,n,r] = ASC_PUP + Beta_AT[n,r] * AT_PUP +  Beta_TD * TD_PUP + Beta_Age_Veh * Age_veh[n];	
	  	end 
	end
	  
	#Calculate the beta fee based on the draws from the distribution
	Beta_parameter=ones(NUM_POINTS,N,R);
	for n=1:N
	  	for r=1:R
	  		Beta_parameter[1,n,r]=0;
	  		Beta_parameter[2,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PSP * Low_inc[n] + Beta_FEE_RES_PSP * Res[n];
	  		Beta_parameter[3,n,r] = Beta_FEE[n,r] + Beta_FEE_INC_PUP * Low_inc[n] + Beta_FEE_RES_PUP * Res[n];
	 	end 
	end
 
    return Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p;
end