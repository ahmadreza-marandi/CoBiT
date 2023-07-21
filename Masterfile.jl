# Julia 1.8.5
#This document contains codes written for the paper
# A. Marandi, V Lurkin (2020), Static Pricing Problems under Mixed Multinomial Logit Demand
#written by Ahmadreza Marandi
# All rights are reserved.
#
#
# using ForwardDiff
using LinearAlgebra
using JuMP
using Gurobi
using HiGHS
# using CPLEX
using MAT
using CPUTime
# using HCubature
using JLD2
using Distributions
# using Cuba
using FiniteDiff
using NLopt
using MosekTools
# using SCS
using Combinatorics
using Dates
struct MyProblem
	model
	p
	y
end



function meshgrid()
	p=collect(0:0.02:1);
	size_p=size(p,1);
	F=zeros(size_p,size_p,size_p);
	for i=1:size_p
		for j=1:size_p
			for k=1:size_p
			# println(io,"$i==$j")
				p_value=[0;p[i]*(UB_p[2]-LB_p[2])+LB_p[2];p[j]*(UB_p[3]-LB_p[3])+LB_p[3];p[k]*(UB_p[4]-LB_p[4])+LB_p[4]];
				F[i,j,k]=discrete_Mixed_logit_function(p_value);
			end
		end
	end
	file = matopen("C:\\Users\\20176914\\surfdrive\\scientific\\University\\Codes\\Julia\\NonlinearProgramming_Virginie\\intel.mat", "w")
	write(file, "F", F)
	# write(file, "F", F)
	close(file)
end
function eye(n)
	Matrix{Float64}(I,n,n);
end
function discrete_Mixed_logit_function_i(p::Vector,i)
	#component of the discrete mixed logit inside the summations for alternatives and customers
	#p is the cost Vector
	#i is the index of the alternative service
	for i_p=1:NUM_POINTS
		for i_n=1:N
			for r=1:R
				if i == (i_p-1)*N*R+(i_n-1)*R+r
					return (sum(exp((Beta_parameter[j,i_n,r]*p[j]+q_parameter[j,i_n,r])-(Beta_parameter[i_p,i_n,r]*p[i_p]+q_parameter[i_p,i_n,r])) for j=1:NUM_POINTS))
				end
			end 
		end
	end
	# f=sum(sum((w_k[r]*exp(Beta_parameter[i,n,r]*p[i]+q_parameter[i,n,r])/(1+ sum(exp(Beta_parameter[j,n,r]*p[j]+q_parameter[j,n,r]) for j=2:NUM_POINTS))) for r=1:R) for n=1:N) 
	# return f;

end
function discrete_Mixed_logit_function(p::Vector)
	#objective function of the discrete mixed logit 
	#p is the cost Vector
	# return sum((p[i_p]*w_k[r])/(sum(exp((Beta_parameter[j,i_n,r]*p[j]+q_parameter[j,i_n,r])-(Beta_parameter[i_p,i_n,r]*p[i_p]+q_parameter[i_p,i_n,r])) for j=1:NUM_POINTS)) for i_p=1:NUM_POINTS for i_n=1:N for r=1:R );
	return sum((p[i_p]*C[(i_p-1)*N*R+(i_n-1)*R+r,i_p]+d[(i_p-1)*N*R+(i_n-1)*R+r])/(sum(exp((Beta_parameter[j,i_n,r]*p[j]+q_parameter[j,i_n,r])-(Beta_parameter[i_p,i_n,r]*p[i_p]+q_parameter[i_p,i_n,r])) for j=1:NUM_POINTS)) for i_p=1:NUM_POINTS for i_n=1:N for r=1:R );
end
function Distribution_function(p::Vector)
	#objective function related to MNL
	f=sum( (p[i_p]*C[(i_p-1)*N*R+(i_n-1)*R+r,i_p]+d[(i_p-1)*N*R+(i_n-1)*R+r])*( exp(Beta_parameter[i_p,i_n]*p[i_p]+q_parameter[i_p,i_n])/(sum(exp(Beta_parameter[j,i_n]*p[j]+q_parameter[j,i_n]) for j=1:NUM_POINTS)))  for i_n=1:N, i_p=1:NUM_POINTS, r=1:R)		
	return f;
end
function homogenious_MNL(time_limit,hatp_1)
	# Algorithm 2 of the paper Li et al. (2019)
	# Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p=Mixed_Logit_n10_random(100);
	# w=1/R*ones(R,1);
	b=Array{Float64,2}(-Beta_parameter[:,1,:]);
	a=Array{Float64,2}(q_parameter[:,1,:]); 
	n=NUM_POINTS;
	######
	# Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p,w_k=Intel_cooperation();
	w=w_k;
	b=Array{Float64,2}(-Beta_parameter[:,1,:]);
	a=Array{Float64,2}(q_parameter[:,1,:]); 
	n=NUM_POINTS;
	#######
	# Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p,w=random_instance_N_1(10);
	# b=Array{Float64,2}(-Beta_parameter[:,1,:]);
	# a=Array{Float64,2}(q_parameter[:,1,:]); 
	# n=NUM_POINTS;
	#considering one customer behaviour

	A=exp.(a);
	m=R;
	

	# hatp_1=1 ./ b[:,1];
	# hatp_1=rand(NUM_POINTS,1);
	hatp_1[1]=0;
	diff=Inf;
	Tolerr=1e-6;
	time_elapsed=0;
	maxf=Inf;
	while diff>Tolerr
		display(" ===$maxf===$hatp_1 =====time=$time_elapsed")
		CPUtic()
		d=zeros(NUM_POINTS,1);
		q_0=zeros(m,1);
		q=zeros(n,m);
		for k=1:m
			q_0[k]=1/(1 + sum( A[j,k] * exp( -b[j,k] * hatp_1[j] ) for j=2:n));
			for i=2:n
				q[i,k] = q_0[k] * A[i,k] * exp( -b[i,k] * hatp_1[i] );
			end
		end
		q_cum= sum(w[k] * q[:,k] for k=1:m);

		for i=2:n
			d[i]=(1 / ( sum( b[i,k] * w[k] * q[i,k] / q_cum[i] for k=1:m ) ) )+
					sum( (w[k] * b[i,k] * q[i,k] / (sum( w[ℓ] * b[i,ℓ] * q[i,ℓ] for ℓ=1:m ))) * (sum( hatp_1[j] * q[j,k] for j=2:n)) for k=1:m ) - hatp_1[i];
		end
		function π_obj(α)
			p=hatp_1+α*d;
			q_0=zeros(m,1);
			q=zeros(n,m);
			for k=1:m
				q_0[k]=1/(1 + sum( A[j,k] * exp( -b[j,k] * p[j] ) for j=2:n));
				for i=2:n
					q[i,k] = q_0[k] * A[i,k] * exp( -b[i,k] * p[i] );
				end
			end
			return sum( w[k] * sum(p[i] * q[i,k] for i=1:n) for  k=1:m);
		end
		opt = Opt(:GN_DIRECT_L, 1);
		opt.lower_bounds = [0];
		opt.upper_bounds = [1];
		# opt.xtol_rel = 1e-;
		opt.maxtime= time_limit-time_elapsed;
		opt.max_objective = π_obj;
		(maxf,α_p,ret) = NLopt.optimize(opt, [1])
		maxf=π_obj(α_p[1]);
		hatp_1=hatp_1+α_p[1]*d;
		diff=norm(α_p[1]*d);
		time_elapsed = time_elapsed+ CPUtoc();
	end
	return hatp_1,discrete_Mixed_logit_function(hatp_1[:,1])
end
function SCIP_degenerate(data)
	model = Model(SCIP.Optimizer)
	set_optimizer_attribute(model, "limits/time", 600)
	# set_optimizer_attribute(model,"println/verblevel",5)
	Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p=data();
	@variable(model,p[1:NUM_POINTS])
	@constraint(model, p .>= LB_p)
	@constraint(model, p .<= UB_p)
	@variable(model, tau)
	@objective(model, Max,tau)
	@NLconstraint(model,tau <= sum(sum((exp(Beta_parameter[i,n]*p[i]+q_parameter[i,n])/(1+ sum(exp(Beta_parameter[j,n]*p[j]+q_parameter[j,n]) for j=2:NUM_POINTS))) for n=1:N) for i=1:NUM_POINTS) )

	sol_scip=JuMP.optimize!(model)
	total_time_scip = MOI.get(model, MOI.SolveTime())
	status_scip=termination_status(model);
	gap_scip=relative_gap(model);
	return total_time_scip, status_scip, gap_scip, value.(p)
end
function Mixed_Logit_distribution(p::Vector,β,i)
	#component of the continous mixed logit inside the integral
	#p is the cost Vector
	#β is the random variable used in the continous mixed logit
	#i is the index of the alternative service

	Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p=Mixed_Logit_10(β);# given β, reading the data.
	NUM_POINTS=size(p,1);
	μ=[-0.788;-32.3];# mean of the random varialbe
	Σ=[(1.06)^2 -12.8;-12.8 (14.2)^2];# covariance matrix of the random varialbe
	f=zeros(N,1);
	# for each customer, calculate the component function
	for n=1:N
		f[n]=(exp(Beta_parameter[i,n]*p[i]+q_parameter[i,n])/ sum(exp(Beta_parameter[j,n]*p[j]+q_parameter[j,n]) for j=1:NUM_POINTS));
	end
	
	return sum(f[n] for n=1:N)*pdf(MvNormal(μ, Σ), β) #continuous

end
function Mixed_logit_function_i(p::Vector, i)
	#component of the continous mixed logit related to the integral
	#p is the cost Vector
	#i is the index of the alternative service
	a=[-3.6;-68.52];# lowerbounds on the integral
	b=[1.94;3.92];#upperbounds on the integral

	## using Cuba
	f=cuhre((x,val)->val[1]=Mixed_Logit_distribution(p,a+x.*(b-a) ,i),2,1)# calculating the integral
	return 1/(prod(b-a)*sum(f[1]));#sum is to make the value scalar
end
function Mixed_logit_function(p::Vector)
	#calculating the objective function with respect to the continuous mixed logit
	# p is the vector of prices
	return sum(p[i]/Mixed_logit_function_i(p, i) for  i=1:NUM_POINTS)
end
function Mixed_logit_function_i_discreteDis(p::Vector, i)
	#component of the discrete mixed logit related to the integral
	#p is the cost Vector
	#i is the index of the alternative service
	μ=[-0.788;-32.3];# mean of the random varialbe
	Σ=[(1.06)^2 -12.8;-12.8 (14.2)^2];# covariance matrix of the random varialbe
	a=[-3.6;-68.52];# lowerbounds on the integral
	b=[1.94;3.92];#upperbounds on the integral
	R=10;
	# for n=1:N
	# 	Prob_β=h[n].weights;
	# 	R_1=size(Prob_β,1);
	# 	R_2=size(Prob_β,2);
	# 	R=maximum([R,R_1,R_2]);
	# end
	# println(io,N)
	f=zeros(N,R,R);
	for n=1:N
		# h = fit(Histogram, (Beta_AT[n,:],Beta_FEE[n,:]),nbins=R);
		# Prob_β=h[n].weights;
		# R_1=size(Prob_β,1);
		# R_2=size(Prob_β,2);
		# println(io,"$R_1 $R_2")
		for i_1=1:R
			for i_2=1:R
				β=[a[1]+(i_1-1)*((b[1]-a[1])/R); a[2]+(i_2-2)*((b[2]-a[2])/R)];
				Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p=Mixed_Logit_50(β);#
				# for each customer, calculate the component function
				f[n,i_1,i_2]=(prod(b-a)/(R^2))*(exp(Beta_parameter[i,n]*p[i]+q_parameter[i,n])/ sum(exp(Beta_parameter[j,n]*p[j]+q_parameter[j,n]) for j=1:NUM_POINTS))*pdf(MvNormal(μ, Σ), β);
			end
		end
	end
	return 1/sum(f)
end
function Mixed_logit_function_discreteDis(p::Vector)
	#calculating the objective function with respect to the continuous mixed logit
	# p is the vector of prices
	# h=Array{Any,1}(undef,N);
	# for n=1:N
 #       h[n] = fit(Histogram, (Beta_AT[n,:],Beta_FEE[n,:]),nbins=R);
 #    end
	return sum(p[i]/Mixed_logit_function_i_discreteDis(p, i) for  i=1:NUM_POINTS)
end
function Mixed_logit_function_nlopt(p::Vector, grad::Vector)
	# the objective function used in NLopt package for continuous mixed logit
	if length(grad) > 0# gradient 
		g = x -> FiniteDiff.finite_difference_gradient(Mixed_logit_function, x);
		grad=g(p);
	end
	val=sum(p[i]/Mixed_logit_function_i(p, i) for i=1:NUM_POINTS);
	return val
end

function Distribution_function_NLopt(p::Vector, grad::Vector)
	# the objective function used in NLopt package for discrete mixed logit
	if length(grad) > 0 #gradient calculations
		g = x -> FiniteDiff.finite_difference_gradient(Distribution_function, x);
		grad=g(p);
	end
	return sum( p[i]*( exp(Beta_parameter[i,n]*p[i]+q_parameter[i,n])/(sum(exp(Beta_parameter[j,n]*p[j]+q_parameter[j,n]) for j=1:NUM_POINTS)))  for n=1:N, i=2:NUM_POINTS)
end

function Distribution_function_i(p::Vector,i)
	# objective fucntion of the discrete mixed logit
	for i_p=1:NUM_POINTS
		for i_n=1:N
			if i == (i_p-1)*N+i_n
				return sum(exp((Beta_parameter[j,i_n]*p[j]-Beta_parameter[i_p,i_n]*p[i_p])+q_parameter[j,i_n]-q_parameter[i_p,i_n]) for j=1:NUM_POINTS)
			end 
		end
	end
end

function P_inr(p::Vector, i::Int64, n::Int64, r::Int64)
	return ( exp(Beta_parameter[i,n,r]*p[i]+q_parameter[i,n,r])/(sum(exp(Beta_parameter[j,n,r]*p[j]+q_parameter[j,n,r]) for j=1:NUM_POINTS)))
end
function derivative_Mixed_general(p::Vector)
	∇=zeros(NUM_POINTS);
	for k=1:NUM_POINTS
		∇[k]= sum(w_k[r,n]*(-Beta_parameter[k,n,r]) * P_inr(p,k,n,r) *((-1/Beta_parameter[k,n,r]) - p[k] +sum(p[i] * P_inr(p,i,n,r) for i=1:NUM_POINTS))    for n=1:N for r=1:R)

	end
	return ∇
end
function P_in(p::Vector, i::Int64, n::Int64)
	return ( exp(Beta_parameter[i,n]*p[i]+q_parameter[i,n])/(sum(exp(Beta_parameter[j,n]*p[j]+q_parameter[j,n]) for j=1:NUM_POINTS)))
end
function derivative_MNLV1(p::Vector)
	∇=zeros(NUM_POINTS);
	for k=1:NUM_POINTS
		∇[k]= sum( P_in(p,k,n) *(1 -(-Beta_parameter[k,n])* p[k] +(-Beta_parameter[k,n])*sum(p[i] * P_in(p,i,n) for i=1:NUM_POINTS))    for n=1:N)

	end
	return ∇
end
function derivative_MNL(p::Vector)
	∇=zeros(NUM_POINTS);
	for k=1:NUM_POINTS
		∇[k]= sum( P_in(p,k,n) *(C[(k-1)*N+n,k] -(-Beta_parameter[k,n])* (C[(k-1)*N+n,:]'*p  +d[(k-1)*N+n]) +(-Beta_parameter[k,n])*sum(((C[(i-1)*N+n,:]'*p  +d[(i-1)*N+n])) * P_in(p,i,n) for i=1:NUM_POINTS))    for n=1:N)

	end
	return ∇
end
function derivative_MNL_discrete(p::Vector)
	∇=zeros(NUM_POINTS);
	for k=1:NUM_POINTS
		∇[k]= sum( P_inr(p,k,n,r) *(1 -(-Beta_parameter[k,n,r])* p[k] +(-Beta_parameter[k,n,r])*sum(p[i] * P_inr(p,i,n,r) for i=1:NUM_POINTS))    for n=1:N for r=1:R)

	end
	return ∇
end

function solve_local_opt_Btree(type_of_problem,time_limit,location)
	#returning the optimal value by conducting branching
	# type_of_problem is whether the problem is "Mixed" (for continuous mixed logit problem), "MNL" (for multinomial logit problem), or "Discrete" (for discrete mixed logit)
	#time_limit is the time limit in seconds
	# cd(location)
	dt=now();
	today_s=Dates.format(dt, "yyyy-mm-dd HH.MM.SS");
	io = open(location*"\\log"*type_of_problem*" N="*string(N)*" I="*string(NUM_POINTS)*" R="*string(R)*" "*today_s*".txt","w")
		# touch(location*"\\log"*today*".txt")
		if type_of_problem=="MNL" 
			obj_fun=Distribution_function;
			obj_fun_i=Distribution_function_i;
			for i=1:NUM_POINTS #lower bounds are similar to LEmma 1 of the MS paper. 
				LB_p[i] = min( minimum(1 ./ abs.(Beta_parameter[i,:])),LB_p[i]);
			end
		elseif type_of_problem=="Mixed" 
			obj_fun=Mixed_logit_function;
			obj_fun_i=Mixed_logit_function_i;
		elseif type_of_problem == "Discrete"
			obj_fun=discrete_Mixed_logit_function;
			obj_fun_i=discrete_Mixed_logit_function_i;
			for i=1:NUM_POINTS #lower bounds are similar to LEmma 1 of the MS paper. 
				LB_p[i] = min( minimum(1 ./ abs.(Beta_parameter[i,:,:])),LB_p[i]);
			end
			# obj_fun=Mixed_logit_function_discreteDis;
			# obj_fun_i=Mixed_logit_function_i_discreteDis;
		end
		
		iteration_history=[-Inf Inf 0];
		break_points=[Float64.(LB_p)];
		push!(break_points,UB_p);
		Number_of_nodes=[0;0];
		node_counter=0;
		global error_tol=0.01
		global rounding_decimals=4
		Tolerr=error_tol*maximum(UB_p-LB_p);
		global relative_err = error_tol
		best_solution=[];
		bestf=-10^18;
		best_old=-Inf;
		A=[eye(NUM_POINTS);
			-eye(NUM_POINTS)];
		b=[LB_p;
			-UB_p];
		Diff=Inf;	
		UB_parent=[UB_p];
		LB_parent=[LB_p];
		fixed_diff=0;
		CompleteLoop=0;
		depth=1;
		# The algorithm do the split, then find new points to break, then find the cuts and do the splitting of the previous sets.  
		Num_breaks_prev=1;
		A_parent=[A];
		b_parent=[b];
		point_parents=[break_points];
		axis_with_large_size=[argmax(UB_p-LB_p)];
		value_axis_with_large_size=[(UB_p-LB_p)[axis_with_large_size]/2];
		Upperbound_value=Inf;
		time_elapsed=0.0;
		first_iteration=1;

		while abs(bestf-Upperbound_value)/abs(bestf)>relative_err && abs(bestf-Upperbound_value)>Tolerr &&time_elapsed<time_limit
			println(io,"====================New_branch======================")
			display("====================New_branch======================")
			# f=open(location*"\\log"*today*".txt", "w")
			# 	println(io,f," ====================New_branch======================")
			# close(f)
			CPUtic()
			best_old=copy(bestf);
			#split the feasible reagon
			m_break=2;
			A_leaf=Array{Any}(undef,m_break,Num_breaks_prev);
			b_leaf=Array{Any}(undef,m_break,Num_breaks_prev);
			LB_p=Array{Any}(undef,m_break,Num_breaks_prev);
			UB_p=Array{Any}(undef,m_break,Num_breaks_prev);
			points_leaf=Array{Any}(undef,Num_breaks_prev);
			Best_solution_leaf=Array{Any}(undef,m_break,Num_breaks_prev);
			Bestf_leaf=ones(m_break,Num_breaks_prev) * (-Inf);
			Best_Upperbound_value=zeros(m_break,Num_breaks_prev);
			status=ones(m_break,Num_breaks_prev);# to track whether we have feasible region. all the regions are assumed to be infeasible at the beginning
			time_elapsed=time_elapsed+CPUtoq();
			println(io,"======$m_break=====$Num_breaks_prev")
			display("======$m_break=====$Num_breaks_prev")
			for j=1:Num_breaks_prev
				for i=1:m_break
					CPUtic()
					
					e=zeros(1,NUM_POINTS);
					e[axis_with_large_size[j]]=1;
					A_lea=[e;-e];
					b_lea=[value_axis_with_large_size[j];-value_axis_with_large_size[j]];
					
					
					A_leaf[i,j]=[A_parent[j];A_lea[i,:]'];b_leaf[i,j]=[b_parent[j];b_lea[i]];
					Lbound_help=Array{Float64}(LB_parent[j]);
					Ubound_help=Array{Float64}(UB_parent[j]);
					if i==1
						Lbound_help[axis_with_large_size[j]]=value_axis_with_large_size[j][1];
						LB_p[i,j]=copy(Lbound_help);
						UB_p[i,j]=copy(Ubound_help);

					else
						Ubound_help[axis_with_large_size[j]]=value_axis_with_large_size[j][1];
						LB_p[i,j]=copy(Lbound_help);
						UB_p[i,j]=copy(Ubound_help);
					end
					# println(io,LB_p[i,j])
					# println(io,UB_p[i,j])
					node_counter += 1 ;
					# println(io,"i=$i....j=$j...$(status[i,j])")
					# println(io,A_leaf[i,j])
					# println(io,b_leaf[i,j])
					time_elapsed=time_elapsed+CPUtoq();
					if time_elapsed<time_limit
						if ((Num_breaks_prev-j)*(2)+2-i) % 10 ==0
							println(io,"==Number of points $(size(break_points,1))===========Remaining node $((Num_breaks_prev-j)*(2)+2-i)========$(round(bestf,digits= rounding_decimals*2))==========$(round(Upperbound_value,digits= rounding_decimals*2)) ========$best_solution===========$(round((Upperbound_value-bestf)/bestf,digits= rounding_decimals*2)*100)%=============$(round(time_elapsed,digits= rounding_decimals*2))===")
							display("==Number of points $(size(break_points,1))===========Remaining node $((Num_breaks_prev-j)*(2)+2-i)========$(round(bestf,digits= rounding_decimals*2))==========$(round(Upperbound_value,digits= rounding_decimals*2)) =================$(round((Upperbound_value-bestf)/bestf,digits= rounding_decimals*2)*100)%=============$(round(time_elapsed,digits= rounding_decimals*2))===")
						end
						CPUtic()
						
						
						time_elapsed=time_elapsed+CPUtoq();
						if time_elapsed<time_limit
							CPUtic()
							trust_Result=trust_region_iteration(A_leaf[i,j],b_leaf[i,j],Tolerr,obj_fun,type_of_problem, rand(NUM_POINTS,1),time_limit-time_elapsed)
							
							if trust_Result == "infeasible"
									status[i,j]=0
									node_counter -=1;
							elseif norm(UB_p[i,j]-LB_p[i,j],Inf) < Tolerr	
								status[i,j]=0
								node_counter -=1;
								Bestf_leaf[i,j]=trust_Result[1]; #for -0.0 to be positive
								Best_solution_leaf[i,j]=abs.(trust_Result[2]);
								Best_Upperbound_value[i,j] = obj_fun(Best_solution_leaf[i,j]);
							else
								# println(io,"==Number of points $(m_break)===========Remaining node =$(sum(status))=== $(size(break_points,1)-m_break)======$(round(bestf,digits= rounding_decimals*2))==========$(round(Upperbound_value,digits= rounding_decimals*2)) ========$best_solution===========$(round((Upperbound_value-bestf)/bestf,digits= rounding_decimals*2)*100)%=============$(round(time_elapsed,digits= rounding_decimals*2))===")
								Bestf_leaf[i,j]=trust_Result[1]; #for -0.0 to be positive
								Best_solution_leaf[i,j]=abs.(trust_Result[2]);
								LB_τ,UB_τ,LB_p[i,j],UB_p[i,j] = bound_identifier_Btree(C, d,A_leaf[i,j],b_leaf[i,j],LB_p[i,j],UB_p[i,j],obj_fun_i,type_of_problem,break_points,N,Sum_size,time_limit-time_elapsed,Best_solution_leaf[i,j]);
							end

							if time_limit<=time_elapsed
								LB_τ="Algorithm is terminated"
								return bestf,best_solution,Upperbound_value,iteration_history,break_points;
							end
							#finding the upperbounds
							if status[i,j]==1 && LB_τ !="Algorithm is terminated"  
								if LB_τ!="solved" 
									# help=overestimator_general(C,d,A_leaf[i,j],b_leaf[i,j],break_points,obj_fun_i,LB_τ,UB_τ,LB_p[i,j],UB_p[i,j],N,time_limit-time_elapsed,bestf);
									# println(io,help)
									help=overestimator_perspective_convex(C,d,A_leaf[i,j],b_leaf[i,j],Best_solution_leaf[i,j],type_of_problem,LB_τ,UB_τ,LB_p[i,j],UB_p[i,j],N,time_limit-time_elapsed,bestf,Best_solution_leaf[i,j],Tolerr);
									# display(help[2])
								else
									help="solved"
								end
								#updating the information
								if help=="solved"
									Best_Upperbound_value[i,j]=round.(Bestf_leaf[i,j],digits=rounding_decimals*2);
									status[i,j]=0;
								elseif help[1] =="Others"
									Best_Upperbound_value[i,j]=round.(Upperbound_value,digits=rounding_decimals*2);
									status[i,j]=1;
								else
									Best_Upperbound_value[i,j]=help[1];
									if Bestf_leaf[i,j]<obj_fun(help[2]) 
										Bestf_leaf[i,j]=obj_fun(help[2]);
										Best_solution_leaf[i,j]=help[2];
									end
								end
								time_elapsed=time_elapsed+CPUtoq();
								
								if Bestf_leaf[i,j] > bestf
									bestf=copy(Bestf_leaf[i,j]);
									best_solution=copy(Best_solution_leaf[i,j]);
									println(io,"===========*$(bestf)==========$(Upperbound_value) ================$(best_solution[:,1])===========$time_elapsed====")
									display("===========*$(bestf)==========$(Upperbound_value) ================$(best_solution[:,1])===========$time_elapsed====")
								end

								
								#terminate branching if the upperbound is at most the same as the best found solution's objective function
								if Best_Upperbound_value[i,j]<= bestf || abs(bestf-Best_Upperbound_value[i,j])/abs(bestf)<relative_err || abs(bestf-Upperbound_value)<Tolerr
									status[i,j]=0;# we don't further branch this leaf as we can't get a better solution there
								end
								##
								#checking whether the found solution is among the existence points, if not put it in the set 
								CPUtic()
								m_break_help=size(break_points,1);
								point_exist=0;
								for m=1:m_break_help
									if maximum(abs.(break_points[m]-Best_solution_leaf[i,j]))<Tolerr #check whether the solution is the same as before or not. if so mark it
										point_exist=1;
										break;
									end
								end
								# randomly choosing a point if the obtained solution is not new
								if point_exist==0 
									push!(break_points,Best_solution_leaf[i,j]);
								end
								time_elapsed=time_elapsed+CPUtoq();
							end
						end
					end
				end

			end
			first_iteration=0; #after the first iteration, we set this value to 0;
			#checking the tolerance found between the previous best and this iteration
			CPUtic()
			if best_old==-Inf
				Diff=bestf;
			else
				Diff=abs(best_old-bestf);
				if Diff<Tolerr
					fixed_diff=fixed_diff+1;
				end
				if bestf>best_old
					fixed_diff=0;
				end
			end
			# checking whether there is a new point added to the list, otherwise report that the points loop has been occurred
			m_break_help=size(break_points);

			if m_break_help==m_break #check whether a new point has been added. If not the loop will be terminated
				CompleteLoop==1;
			end
			if time_elapsed<time_limit
				# reconstruction of the parents: we go to the next branching
				ind_feasible=findall(status.==1);
				Num_breaks_prev=size(ind_feasible,1);
				# push!(Number_of_nodes,Num_breaks_prev)
				if Num_breaks_prev>0
					A_parent=Array{Any}(undef,Num_breaks_prev,1);
					b_parent=Array{Any}(undef,Num_breaks_prev,1);
					UB_parent=Array{Any}(undef,Num_breaks_prev,1);
					LB_parent=Array{Any}(undef,Num_breaks_prev,1);
					axis_with_large_size=Array{Any}(undef,Num_breaks_prev,1);
					value_axis_with_large_size=Array{Any}(undef,Num_breaks_prev,1);
					for i=1:Num_breaks_prev
							A_parent[i]=A_leaf[ind_feasible[i]];
							b_parent[i]=b_leaf[ind_feasible[i]];
							LB_parent[i]=LB_p[ind_feasible[i]];
							UB_parent[i]=UB_p[ind_feasible[i]];
							# println(io,argmax(UB_p[ind_feasible[i]]-LB_p[ind_feasible[i]]))
							axis_with_large_size[i]=argmax(UB_p[ind_feasible[i]]-LB_p[ind_feasible[i]]);
							value_axis_with_large_size[i]=(UB_p[ind_feasible[i]]-LB_p[ind_feasible[i]])[axis_with_large_size[i]]/2+(LB_p[ind_feasible[i]])[axis_with_large_size[i]];
					end
					depth=depth+1;
					Upperbound_value=min(Upperbound_value,maximum(Best_Upperbound_value))
				else
					Upperbound_value=bestf;
				end

				time_elapsed=time_elapsed+CPUtoq();
			end
			println(io,"===========$(bestf)==========$(Upperbound_value) =======$best_solution============$time_elapsed====")
			display("===========$(bestf)==========$(Upperbound_value) =================$time_elapsed====")
			iteration_history=[iteration_history; bestf Upperbound_value time_elapsed];
			push!(Number_of_nodes,copy(node_counter))
			node_counter=0;
		end
	close(io)
	return bestf,best_solution,Upperbound_value,iteration_history,break_points,Number_of_nodes;

end

function trust_region_iteration(A,b,Tolerr,obj_fun::Function,type_of_problem ,initial_point,time_limit)
	# do the iteration of the trust-region
	# the feasible region of Ap⩾b 
	
		function g(x)#gradient function
			if type_of_problem == "MNL"
				# return FiniteDiff.finite_difference_gradient(obj_fun, x);
				return derivative_MNL(x)
			elseif type_of_problem == "Discrete"
				return derivative_MNL_discrete(x)
			end
		end
		

	

	NUM_POINTS=size(A,2);
	constraint_size=size(A,1);
    f0=0;
	## initialization of the radius
	# Model_r = Model(CPLEX.Optimizer);
	# set_optimizer_attribute(Model_r,"CPX_PARAM_SCRIND",0);
	# set_optimizer_attribute(Model_r,"CPX_PARAM_PERIND",0);
	# set_optimizer_attribute(Model_r,"CPX_PARAM_TILIM",time_limit);
	# Model_r= Model(Gurobi.Optimizer)
	# set_optimizer_attribute(Model_r,"LogToConsole",0)
	# set_optimizer_attribute(Model_r,"CSClientLog",0);
	# set_optimizer_attribute(Model_r,"OutputFlag",0);
	# set_optimizer_attribute(Model_r,"TimeLimit",time_limit);
	Model_r=Model(HiGHS.Optimizer)
	set_optimizer_attribute(Model_r,"log_to_console",false)
	set_optimizer_attribute(Model_r,"time_limit",time_limit);
	@variable(Model_r,p[1:NUM_POINTS]>=0)
	@constraint(Model_r,A*p.>=b)
	@objective(Model_r,Max,0)
	JuMP.optimize!(Model_r)
	if termination_status(Model_r)== MOI.INFEASIBLE
		return "infeasible"
	end
	# U=sum(value.(p));
	initial_p=copy(value.(p));
	# Model_r = Model(CPLEX.Optimizer);
	# set_optimizer_attribute(Model_r,"CPX_PARAM_SCRIND",0);
	# set_optimizer_attribute(Model_r,"CPX_PARAM_PERIND",0);
	# set_optimizer_attribute(Model_r,"CPX_PARAM_TILIM",time_limit);
	# @variable(Model_r,p[1:NUM_POINTS]>=0)
	# @constraint(Model_r,A*p.>=b)
	# @objective(Model_r,Min,sum(p))
	# JuMP.optimize!(Model_r)
	# L=sum(value.(p));
	# initial_p=0.5*(value.(p)+copy(initial_p)); #initial point in the middle of the feasible region. Since it is polytope it lies in the set
	###
	# if minimum(A*initial_point-b)<-0.0000001
	# else
	# 	initial_p=copy(initial_point);
	# end
    f0=obj_fun(initial_p);
    radius=1;
    time_elapsed=0;
  	f1=Inf;
    while abs(f0-f1)>Tolerr && norm(g(initial_p))>0.01*relative_err &&time_limit>time_elapsed &&radius >10^(-4)
		CPUtic()
		# display(A)
		# display(b)
		# display(radius)
		# display(initial_p)
		p1=trust_region(initial_p,obj_fun,type_of_problem,A,b,radius,time_limit-time_elapsed)
		if p1 != "infeasible"
			f1=obj_fun(p1);
		else
			p1=copy(initial_p);
			f0=copy(f1);
		end
		time_elapsed=time_elapsed+CPUtoq();
		f1=obj_fun(p1);
		while f1>f0 && time_limit>time_elapsed

			CPUtic()
			initial_p=copy(p1);
			f0=copy(f1);
			p1=trust_region(initial_p,obj_fun,type_of_problem,A,b,radius,time_limit-time_elapsed);
			if p1 != "infeasible"
				f1=obj_fun(p1);
			else
				p1=copy(initial_p);
				f0=copy(f1);
			end
			radius=1;;
			time_elapsed=time_elapsed+CPUtoq();
		end
		radius=radius/10;
		
    end
    return f0,initial_p;
end


function overestimator_perspective_convex(C, d, A,b,P_0,type_of_problem,LB_τ,UB_τ,LB_p,UB_p,N,time_limit,lowerbound,point,Tolerr)
	if LB_τ=="solved" #because of the rounding errors we get infeasible. So, the area is small
		return "solved"
	end
	if sum((UB_τ-LB_τ))<error_tol
		return "solved"
	end
	if sum((UB_p-LB_p))<error_tol
		return "solved"
	end

	
	A=round.(A,digits=rounding_decimals*2);
	b=round.(b,digits=rounding_decimals*2);
	LB_τ=round.(LB_τ,digits=rounding_decimals*2);
	UB_τ=round.(UB_τ,digits= rounding_decimals*2);
	LB_p=round.(LB_p,digits= rounding_decimals*2);
	UB_p=round.(UB_p,digits= rounding_decimals*2);
	Sum_size = size(C,1);

	##CPLEX
	# LP = Model(CPLEX.Optimizer);
	# set_optimizer_attribute(LP,"CPX_PARAM_SCRIND",0);
	# set_optimizer_attribute(LP,"CPX_PARAM_PREDUAL",1);
	# set_optimizer_attribute(LP,"CPX_PARAM_TILIM",time_limit);
	# Ipopt
	# LP = Model(Ipopt.Optimizer);
	# set_optimizer_attribute(LP, "print_level",0)# for ipopt
	# set_optimizer_attribute(LP, "max_cpu_time", time_limit)
	# set_optimizer_attribute(LP, "tol", Tolerr)
	#for gurobi
	# set_optimizer_attribute(LP,"DualReductions",0);
	# set_optimizer_attribute(LP,"LogToConsole",0)
	# set_optimizer_attribute(LP,"CSClientLog",0);
	# set_optimizer_attribute(LP,"OutputFlag",0);
	# set_optimizer_attribute(LP,"TimeLimit",time_limit);
	# NLopt
	# LP = Model(NLopt.Optimizer);
	# set_optimizer_attribute(LP, "algorithm", :LN_COBYLA)
	# set_optimizer_attribute(LP, "ftol_rel",0.00001)
	# set_optimizer_attribute(LP, "xtol_rel",0.00001)
	# conic solver
	LP = Model(Mosek.Optimizer);
	set_optimizer_attribute(LP,"MSK_IPAR_LOG",0)
	set_optimizer_attribute(LP,"MSK_IPAR_LOG_INTPNT",0)
	set_optimizer_attribute(LP,"MAX_NUM_WARNINGS",0)
	set_optimizer_attribute(LP,"MSK_DPAR_INTPNT_CO_TOL_REL_GAP", error_tol)

	# LP = Model(SCS.Optimizer);
	# set_optimizer_attribute(LP,"verbose",0)
	# set_optimizer_attribute(LP,"time_limit_secs", time_limit)
	# set_optimizer_attribute(LP,"eps_rel", Tolerr)
	P_size=size(A,2);
	@variable(LP, W[1:Sum_size,1:P_size]>=0);

	@variable(LP, 1>=τ[1:Sum_size]>=0)

	violance=0;
	@variable(LP, p[1:P_size]>=0)

	@constraint(LP, A*p.>=b)


	for i=1:Sum_size
		@constraint(LP, τ[i]<=UB_τ[i]+violance);
		@constraint(LP, τ[i]>=LB_τ[i]-violance);
		for j=1:P_size
			@constraint(LP,W[i,j] >=LB_τ[i]*p[j]+τ[i]*LB_p[j]-LB_τ[i]*LB_p[j])
			@constraint(LP,W[i,j] >=UB_τ[i]*p[j]+τ[i]*UB_p[j]-UB_τ[i]*UB_p[j])
			@constraint(LP,W[i,j] <=UB_τ[i]*p[j]+τ[i]*LB_p[j]-UB_τ[i]*LB_p[j])
			@constraint(LP,W[i,j] <=τ[i]*UB_p[j]+LB_τ[i]*p[j]-LB_τ[i]*UB_p[j])
			@constraint(LP,W[i,j] -UB_τ[i]*p[j] <=violance)
			@constraint(LP,W[i,j] -LB_τ[i]*p[j] >=-violance)
			@constraint(LP, W[i,j] <=τ[i]*UB_p[j]+violance);
			@constraint(LP, W[i,j] >=τ[i]*LB_p[j]-violance);
			@constraint(LP, W[i,j]<=UB_p[j]*UB_τ[i])
			@constraint(LP, W[i,j]>=LB_p[j]*LB_τ[i])
		end
		
		for t=1:size(A,1)
			@constraint(LP, sum(A[t,j]*W[i,j] for j=1:P_size)-b[t]*τ[i]>=LB_τ[i]*(sum(A[t,j]*p[j] for j=1:P_size)-b[t])-violance)
			@constraint(LP, sum(A[t,j]*W[i,j] for j=1:P_size)-b[t]*τ[i]<=UB_τ[i]*(sum(A[t,j]*p[j] for j=1:P_size)-b[t])+violance)
			@constraint(LP, sum(A[t,j]*W[i,j] for j=1:P_size)-sum(A[t,j]*p[j]*LB_τ[i] for j=1:P_size)>=(τ[i]- LB_τ[i])*b[t]-violance)
			@constraint(LP, sum(A[t,j]*W[i,j] for j=1:P_size)-sum(A[t,j]*p[j]*UB_τ[i] for j=1:P_size)<=(τ[i]- UB_τ[i])*b[t]+violance)
		end
		
	end
	############ MNL ###############
	@variable(LP, θ[1:Sum_size,1:P_size]>=0)
	# θ = Variable(Sum_size,P_size);
	# LP.constraints += [θ>=0	];
	for i=1:P_size
		for n=1:N
			# LP.constraints += [sum(exp(Beta_parameter[j,n]*(W[(i-1)*N+n , j] / τ[(i-1)*N+n])+q_parameter[j,n]-(Beta_parameter[i,n]*(W[(i-1)*N+n , i]/ τ[(i-1)*N+n])+q_parameter[i,n])) for j=1:NUM_POINTS)*τ[(i-1)*N+n]<=1];
			# @NLconstraint(LP, sum(exp(Beta_parameter[j,n]*(W[(i-1)*N+n , j] / τ[(i-1)*N+n])+q_parameter[j,n]-(Beta_parameter[i,n]*(W[(i-1)*N+n , i]/ τ[(i-1)*N+n])+q_parameter[i,n])) for j=1:NUM_POINTS)*τ[(i-1)*N+n]<=1)
			
			if type_of_problem == "MNL"
				@constraint(LP, sum( θ[(i-1)*N+n,j] for j=1:NUM_POINTS)<=1)
				for j=1:P_size
					@constraint(LP, [Beta_parameter[j,n]*(W[(i-1)*N+n , j] )-(Beta_parameter[i,n]*(W[(i-1)*N+n , i])) + (q_parameter[j,n]-(q_parameter[i,n]))*τ[(i-1)*N+n],τ[(i-1)*N+n],θ[(i-1)*N+n,j] ] in MOI.ExponentialCone())
				end
			elseif type_of_problem == "Discrete"
				for r=1:R
					@constraint(LP, sum( θ[(i-1)*N*R+(n-1)*R+r,j] for j=1:NUM_POINTS)<=1)
					for j=1:P_size
						@constraint(LP, [Beta_parameter[j,n,r]*(W[(i-1)*N*R+(n-1)*R+r , j] )-(Beta_parameter[i,n,r]*(W[(i-1)*N*R+(n-1)*R+r, i])) + (q_parameter[j,n,r]-(q_parameter[i,n,r]))*τ[(i-1)*N*R+(n-1)*R+r],τ[(i-1)*N*R+(n-1)*R+r],θ[(i-1)*N*R+(n-1)*R+r,j] ] in MOI.ExponentialCone())
					end
				end
			end
		end 
	end

	###################
	@objective(LP, Max, sum(sum(C[i,j]*W[i,j] for j=1:P_size) + d[i]*τ[i] for i=1:Sum_size))
	@constraint(LP, sum(sum(C[i,j]*W[i,j] for j=1:P_size) + d[i]*τ[i] for i=1:Sum_size) >= lowerbound)
	# solve!(LP, Mosek.Optimizer; silent_solver = true)
	sol=JuMP.optimize!(LP);
	# display(termination_status(LP))
	# display(solution_summary(LP))
	# println(io,termination_status(LP))
	# println(io,value.(p))
	# println(io,value.(W))
	if termination_status(LP)== MOI.INFEASIBLE || termination_status(LP)==MOI.INFEASIBLE_OR_UNBOUNDED || termination_status(LP)==MOI.LOCALLY_INFEASIBLE
		return "solved"# this is because the infeasibility comes from rounding errors
	end
	if termination_status(LP)== MOI.NUMERICAL_ERROR
		return "Others"
	end
	if termination_status(LP) == MOI.SLOW_PROGRESS
		summary = solution_summary(LP);
		return summary.dual_objective_value, value.(p)
	end
	return objective_value(LP), value.(p)

end


function trust_region(P_k,obj_fun::Function,type_of_problem,A,b,radius,time_limit)
	# Ap ⫺ b is the feasible set 
	# P_k is the obtained feasible solution in the previous round
	# radius shows the size of the ball
	# file = matopen("C:\\Users\\20176914\\surfdrive\\scientific\\University\\Codes\\Julia\\NonlinearProgramming_Virginie\\debug_trust.mat", "w")
	# write(file, "A", A)
	# write(file, "b", b)
	# write(file, "p_0", P_k)
	# write(file, "r", radius)
	# close(file)
	function g(x)#gradient function
		if type_of_problem == "MNL"
			# return FiniteDiff.finite_difference_gradient(obj_fun, x);
			return derivative_MNL(x)
		elseif type_of_problem == "Discrete"
			return derivative_MNL_discrete(x)
		end
	end
	P_size=size(A,2);
	# Model_k = Model(CPLEX.Optimizer);
	# set_optimizer_attribute(Model_k,"CPX_PARAM_SCRIND",0);
	# # set_optimizer_attribute(Model_k,"CPX_PARAM_OBJLLIM",-Inf)
	# set_optimizer_attribute(Model_k,"CPX_PARAM_TILIM",time_limit);

	# Model_k = Model(Gurobi.Optimizer)
	# set_optimizer_attribute(Model_k,"LogToConsole",0)
	# set_optimizer_attribute(Model_k,"CSClientLog",0);
	# set_optimizer_attribute(Model_k,"OutputFlag",0);
	# set_optimizer_attribute(Model_k,"TimeLimit",time_limit);
	Model_k=Model(HiGHS.Optimizer)
	set_optimizer_attribute(Model_k,"log_to_console",false)
	set_optimizer_attribute(Model_k,"time_limit",time_limit);
	@variable(Model_k,p[1:P_size]>=0)
	@constraint(Model_k,A*p.>=b)
	grad=g(P_k);
	if norm(grad)<=10^10
		@objective(Model_k,Max,sum((grad).*(p)));
	else
		@objective(Model_k,Max,sum(rand(P_size,1).*(p)));
	end
	@variable(Model_k, auxABS[1:P_size]>=0);
	@constraint(Model_k, sum(auxABS)<=radius)
	for i=1:P_size
		@constraint(Model_k, auxABS[i]>=p[i]-P_k[i])
		@constraint(Model_k, -auxABS[i]<=p[i]-P_k[i])
	end
	JuMP.optimize!(Model_k)
	# display(termination_status(Model_k))
	if termination_status(Model_k) == MOI.OPTIMAL
		P_new=value.(p);
	else
		P_new="infeasible"
	end
	return P_new
end


function bound_identifier_Btree(C, d, A,b,LB_p,UB_p,obj_fun_i::Function,type_of_problem,P_0,N,Sum_size,time_limit,P_current)
	#bound identifier for τ and p
	#Ap⫺b
	#obj_fun_i is the objective function
	#N is the number of customer
	# A=round.(A,digits= rounding_decimals*2);
	# b=round.(b,digits= rounding_decimals*2);
	Num_points=size(P_0,1)
	P_size=size(A,2);
	# println(io,typeof(P_size))
	p_init=copy(P_current);
	LB_τ=zeros(Sum_size);
	UB_τ=zeros(Sum_size);
	Err = error_tol;
	for i=1:P_size
		if abs(UB_p[i]-LB_p[i]) <= Err
			LB_p[i]=P_current[i];
			UB_p[i]=copy(LB_p[i]);
		end
		
	end
	
	#we know that τ_in = 1/f_in. So, to find the min we can minize 1/f_in, as it is convex, and to find the
	#max we can simply minimize f_in and then use the argmin as the argmax for τ. To solve the convex min we use ipopt
	for i_p=1:NUM_POINTS
		for i_n=1:N
			for r=1:R
				i= (i_p-1)*N*R+(i_n-1)*R+r;
				model=Model(HiGHS.Optimizer)
				set_optimizer_attribute(model,"log_to_console",false)
				set_optimizer_attribute(model,"time_limit",time_limit);
				@variable(model,p[1:NUM_POINTS]>=0)
				@constraint(model, A*p .>=b)
				@constraint(model, p .>=LB_p)
				@constraint(model, p .<=UB_p)

				for j=1:NUM_POINTS
					if type_of_problem == "MNL"
						@objective(model, Max, Beta_parameter[j,i_n]*p[j]+q_parameter[j,i_n] - (Beta_parameter[i_p,i_n]*p[i_p]+q_parameter[i_p,i_n]))
					elseif type_of_problem == "Discrete"
						@objective(model, Max, Beta_parameter[j,i_n,r]*p[j]+q_parameter[j,i_n,r] - (Beta_parameter[i_p,i_n,r]*p[i_p]+q_parameter[i_p,i_n,r]))	
					end
					JuMP.optimize!(model)
					if termination_status(model) !=  MOI.INFEASIBLE
						LB_τ[i] += exp( objective_value(model))
					else
						LB_τ[i] = obj_fun_i(p_init,i)
					end
				end
				LB_τ[i] = 1/LB_τ[i];
			
				τ_up=τ_upperbound_conic(C,d,A,b,LB_p,UB_p,P_size,P_0,i_n,i_p,r,time_limit, p_init,Err,type_of_problem)
				
				if τ_up != "infeasible"
					UB_τ[i]=round(1/τ_up,digits= rounding_decimals*2)
					# println(io,UB_τ[i]-LB_τ[i])
				else
					UB_τ[i]=1/obj_fun_i(p_init,i);
					LB_τ[i]=copy(UB_τ[i]);
				end

				if UB_τ[i] - LB_τ[i] < Err
					UB_τ[i]=1/obj_fun_i(p_init,i);
					LB_τ[i]=copy(UB_τ[i]);
				end
				
			end
		end
	end

	return LB_τ,UB_τ,LB_p,UB_p;
end

function τ_upperbound_conic(C,d,A,b,LB_p,UB_p,P_size,P_0,i_n,i_p,r,time_limit, p_init,Err,type_of_problem)

	###JuMP
	m =Model(Mosek.Optimizer)
	set_optimizer_attribute(m,"MSK_IPAR_LOG",0)
	set_optimizer_attribute(m,"MSK_IPAR_LOG_INTPNT",0)
	set_optimizer_attribute(m,"MAX_NUM_WARNINGS",0)
	set_optimizer_attribute(m,"MSK_DPAR_INTPNT_CO_TOL_REL_GAP", Err)
	# m = Model(SCS.Optimizer);
	# set_optimizer_attribute(m,"verbose",0)
	# set_optimizer_attribute(m,"time_limit_secs", time_limit)
	# set_optimizer_attribute(m,"eps_rel", Err)

	@variable(m,p[1:NUM_POINTS]>=0)
	@variable(m,beta[1:NUM_POINTS]>=0)
	@constraint(m, A*p .>= b)
	@objective(m, Min, sum(beta))
	if type_of_problem == "MNL"
		@expression(m, u[ j=1:NUM_POINTS], Beta_parameter[j,i_n]*p[j]+q_parameter[j,i_n] - (Beta_parameter[i_p,i_n]*p[i_p]+q_parameter[i_p,i_n]))

	elseif type_of_problem == "Discrete"
		@expression(m, u[ j=1:NUM_POINTS], Beta_parameter[j,i_n,r]*p[j]+q_parameter[j,i_n,r] - (Beta_parameter[i_p,i_n,r]*p[i_p]+q_parameter[i_p,i_n,r]))
	end
	for j = 1:NUM_POINTS
		@constraint(m, [u[j],1,beta[j]] in MOI.ExponentialCone())
		# @NLconstraint(m, exp(u[j]) <= beta[j])
	end
	
	sol=JuMP.optimize!(m);
	# println(io,termination_status(m))
	if termination_status(m) != MOI.INFEASIBLE
		return objective_value(m)
	else
		return "infeasible"
	end 
	
	##Convex.jl
	# p=Variable(NUM_POINTS);
	# Constraint=[p >=0 , A*p >= b];
	# problem = minimize(sum(exp(Beta_parameter[j,i_n]*p[j]+q_parameter[j,i_n] - (Beta_parameter[i_p,i_n]*p[i_p]+q_parameter[i_p,i_n])) for j=1:NUM_POINTS), Constraint)
	# solve!(problem, Mosek.Optimizer(); silent_solver = true)

	# if problem.status != :Infeasible
	# 	return problem.optval
	# else
	# 	return "infeasible"
	# end 
end



function τ_upperbound(C,d,A,b,LB_p,UB_p,P_size,P_0,obj_fun_i,i,time_limit, p_init)
	### for only the denominator using our trust_Region
	function upper_help(x)
       return -obj_fun_i(x,i)
    end
    # function g!(G,x)
    # 	m=size(A,1);
    # 	for ind_i=1:m
    # 		G[ind_i]=-A[ind_i,:]'*x+b[ind_i]
    # 	end
    # 	for ind_i=1:P_size
    # 		G[m+ind_i]=x[ind_i]-UB_p[ind_i]
    # 		G[m+P_size+ind_i]=-x[ind_i]+LB_p[ind_i]
    # 	end
    # end
    # res = Optim.optimize(upper_help, g!, p_init, NewtonTrustRegion())
    # minf=res.minimum;
    # p=res.minimizer;

    # Tolerr=0.0001
    minf,p=trust_region_iteration([A; eye(P_size);-eye(P_size)],[b; LB_p;-UB_p],Tolerr,upper_help,p_init,time_limit) 
	# minf,p=local_solver_NLopt([A; eye(P_size);-eye(P_size)],[b; LB_p;-UB_p],time_limit,Tolerr,upper_help,p_init)
	# minf,p=local_solver_IPOPT_MNL_i([A; eye(P_size);-eye(P_size)],[b; LB_p;-UB_p],time_limit,i)
	### for only the denominator using NLopt
	# function bound_f_nlopt(x::Vector, grad::Vector)
	# 	if length(grad) > 0# gradient 
	# 		g = x -> FiniteDiff.finite_difference_gradient(x->(-obj_fun_i(x,i)), x);
	# 		grad=g(p);
	# 	end
	# 	val=-obj_fun_i(x,i);
	# 	return val
	# end
	# function bound_myconstraint(x::Vector, grad::Vector, a, b)
	#     if length(grad) > 0
	#     	grad=a;
	#     end
	#     a'*x-b
	# end
	# # A_nlop=A;
	# # b_nlop=b;
	# # # file = matopen("C:\\Users\\20176914\\surfdrive\\scientific\\University\\Codes\\Julia\\NonlinearProgramming_Virginie\\debug_trust.mat", "w")
	# # # write(file, "A", A)
	# # # write(file, "b", b)
	# # # write(file, "p_0", p_init)
	# # # write(file, "r", r)
	# # # close(file)
	# P_size=size(A,2);
	# # opt = Opt(:LN_COBYLA, P_size);
	# opt = Opt(:LD_MMA, P_size);
	# opt.ftol_rel=0.00001;
	# opt.xtol_rel=0.00001;
	# opt.lower_bounds = LB_p;
	# opt.upper_bounds = UB_p;
	# for i_cons=1:size(A,1)
	# 	inequality_constraint!(opt, (x,g) -> bound_myconstraint(x,g,A[i_cons,:],b[i_cons]), 1e-8)
	# end
	# opt.maxtime= time_limit;
	# opt.max_objective = bound_f_nlopt
	# (minf,p,ret) = NLopt.optimize(opt,p_init)
	# return minf,p

	####################################
	###################################
	################################
	#for nominator and denominator 
	# Num_points=size(P_0,1);
	# LP = Model(CPLEX.Optimizer);
	# set_optimizer_attribute(LP,"CPX_PARAM_SCRIND",0);
	# set_optimizer_attribute(LP,"CPX_PARAM_PERIND",0);
	# set_optimizer_attribute(LP,"CPX_PARAM_TILIM",time_limit);
	# ##gurobi
	# # LP = Model(Gurobi.Optimizer);
	# # set_optimizer_attribute(LP,"CSClientLog",0);
	# # set_optimizer_attribute(LP,"OutputFlag",0);
	# # set_optimizer_attribute(LP,"TimeLimit",time_limit);
	# @variable(LP, W[1:P_size]>=0);
	# @variable(LP, τ>=0);
	# violance=0;
	# @variable(LP, p[1:P_size]>=0);
	# @constraint(LP, A*p.>=b);
	# @constraint(LP, A*W .>= b*τ);

	# for k=1:Num_points
	# 	p_0=P_0[k];
	# 	function f_derivative(p)
	# 		f=x->obj_fun_i(x,(i-1)*N + n);
	# 		g=x->FiniteDiff.finite_difference_gradient(f, x);
	# 		return g(p);
	# 	end
	# 	Matrix_f_Derivative=f_derivative(p_0);
	# 	value_f_i=obj_fun_i(p_0,(i-1)*N + n);
	# 	@constraint(LP, τ*(value_f_i)+ sum((W[j]*Matrix_f_Derivative[j])  for j = 1:P_size ) -sum(((p_0[j]*Matrix_f_Derivative[j])*τ)  for j = 1:P_size )  <= C[(i-1)*N + n,:]'*p+d[(i-1)*N + n])# using the perspective form
	# end

	# @objective(LP,Max,τ)
	# sol=JuMP.optimize!(LP);
	# # if termination_status(LP)== MOI.INFEASIBLE

	# 	# file = matopen("C:\\Users\\20176914\\surfdrive\\scientific\\University\\Codes\\Julia\\NonlinearProgramming_Virginie\\debug_trust.mat", "w")
	# 	# write(file, "A", A)
	# 	# write(file, "b", b)
	# 	# write(file, "P_0", P_0)
	# 	# write(file, "n", n)
	# 	# write(file, "i", i)
	# 	# close(file)
	# # else
	# return value(τ) , value.(p)[:]
end
function τ_lowerbound(C,d,A,b,P_size,P_0,obj_fun_i,i,time_limit)
	Num_points=size(P_0,1);
	# LP = Model(CPLEX.Optimizer);
	# set_optimizer_attribute(LP,"CPX_PARAM_SCRIND",0);
	# set_optimizer_attribute(LP,"CPX_PARAM_PERIND",0);
	# set_optimizer_attribute(LP,"CPX_PARAM_TILIM",time_limit);
	##gurobi
	# LP = Model(Gurobi.Optimizer)
	# set_optimizer_attribute(LP,"LogToConsole",0)
	# set_optimizer_attribute(LP,"CSClientLog",0);
	# set_optimizer_attribute(LP,"OutputFlag",0);
	# set_optimizer_attribute(LP,"TimeLimit",time_limit);
	LP=Model(HiGHS.Optimizer)
	set_optimizer_attribute(LP,"log_to_console",false)
	set_optimizer_attribute(LP,"time_limit",time_limit);
	@variable(LP, W[1:P_size]>=0);
	@variable(LP, τ>=0);
	violance=0;
	@variable(LP, p[1:P_size]>=0);
	@constraint(LP, A*p.>=b);
	@constraint(LP, A*W .>= b*τ);

	for k=1:Num_points
		p_0=P_0[k];
		function f_derivative(p)
			f=x->obj_fun_i(x,i);
			g=x->FiniteDiff.finite_difference_gradient(f, x);
			return g(p);
		end
		Matrix_f_Derivative=f_derivative(p_0);
		value_f_i=obj_fun_i(p_0,i);
		@constraint(LP, τ*(value_f_i)+ sum((W[j]*Matrix_f_Derivative[j])  for j = 1:P_size ) -sum(((p_0[j]*Matrix_f_Derivative[j])*τ)  for j = 1:P_size )  <= 1)# using the perspective form
	end

	@objective(LP,Min,τ)
	sol=JuMP.optimize!(LP);
	# println(io,termination_status(LP))
	# if termination_status(LP)== MOI.INFEASIBLE

		# file = matopen("C:\\Users\\20176914\\surfdrive\\scientific\\University\\Codes\\Julia\\NonlinearProgramming_Virginie\\debug_trust.mat", "w")
		# write(file, "A", A)
		# write(file, "b", b)
		# write(file, "P_0", P_0)
		# write(file, "n", n)
		# write(file, "i", i)
		# close(file)
	# else
	return value(τ) , round.(value.(p)[:],digits= rounding_decimals*2)
end


function solve_nlopt(time_limit)
	#solving the problem using NLopt
	P_size=size(Beta_parameter,1);
	opt = Opt(:GN_DIRECT_L, P_size);
	opt.lower_bounds = LB_p;
	opt.upper_bounds = UB_p;
	opt.maxtime= time_limit;
	opt.max_objective = Mixed_logit_function_nlopt
	t=@timed (minf,p,ret) = NLopt.optimize(opt, [0;0;0])
	return t
end


function local_solver_NLopt(A,b,time_limit,tolerr,obj_fun,p_init)
	function myconstraint(x::Vector, grad::Vector, a, b)
	    if length(grad) > 0
	    	grad=a;
	    end
	    a'*x-b
	end
	function f_NLopt(x::Vector, grad::Vector)
		if length(grad) > 0# gradient 
			g = x -> FiniteDiff.finite_difference_gradient(-obj_fun(x), x);
			grad=g(p);
		end
		val=-obj_fun(x);
		return val
	end
	P_size=size(A,2);
	opt = Opt(:LN_COBYLA, P_size);
	# opt = Opt(:LD_MMA, P_size);
	opt.ftol_rel=tolerr;
	opt.xtol_rel=tolerr;
	# opt.lower_bounds = LB_p;
	# opt.upper_bounds = UB_p;
	for i_cons=1:size(A,1)
		inequality_constraint!(opt, (x,g) -> myconstraint(x,g,A[i_cons,:],b[i_cons]), 1e-8)
	end
	opt.maxtime= time_limit;
	opt.max_objective = f_NLopt
	(minf,p,ret) = NLopt.optimize(opt,p_init)
	return -minf, round.(p,digits= rounding_decimals*2)
end
function local_solver_MadNLP(A,b,time_limit, initial_point)
	NUM_POINTS=size(A,2);
	model = Model(()->MadNLP.Optimizer(print_level=MadNLP.INFO,max_iter=100))
	@variable(model, p[1:NUM_POINTS]>=0, start = initial_point)
	@constraint(model, A*p .>= b)
	@variable(model, obj)
	@objective(model, Max, obj)
	@NLconstraint(model, obj <=sum( p[i]*( exp(Beta_parameter[i,n]*p[i]+q_parameter[i,n])/(sum(exp(Beta_parameter[j,n]*p[j]+q_parameter[j,n]) for j=1:NUM_POINTS)))  for n=1:N, i=2:NUM_POINTS))
	sol=JuMP.optimize!(model)
	return value(obj),value.(p);
end