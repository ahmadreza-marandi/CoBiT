# Static_Pricing_problem
Codes related to the paper: A. Marandi, V Lurkin (2023), An exact algorithm for the static pricing problem under discrete mixed logit demand (https://arxiv.org/abs/2005.07482).

In case of using these codes, you are obliged to properly cite the paper.  

## Table of contents
* [Technologies](#technologies)
* [Data](#data)
* [Functions and their features](#functions-and-their-features)
* [Examples](#examples)

## Technologies
For this project, we use multiple Julia packages. The code has been updated and tested on Julia 1.8.0. The packages and the tested versions are listed below:
* ForwardDiff v0.10.34 (https://github.com/JuliaDiff/ForwardDiff.jl)
* LinearAlgebra (https://github.com/JuliaLang/julia/tree/master/stdlib/LinearAlgebra)
* JuMP v1.6.0 (https://jump.dev/)
* Gurobi v0.11.5 ([https://github.com/jump-dev/Gurobi.jl]): make sure to properly install Gurobi beforehand
* MosekTools 0.13.2 ([https://github.com/jump-dev/MosekTools.jl]): make sure to instal Mosek beforehand
* MAT v0.10.3 (https://github.com/JuliaIO/MAT.jl): used to import and export some parameters
* CPUTime v1.0.0 (https://github.com/schmrlng/CPUTime.jl)
* JLD2 v0.4.29 (https://github.com/JuliaIO/JLD2.jl): used for some importing and exporting
* Distributions v0.25.80 (https://github.com/JuliaStats/Distributions.jl)
* FiniteDiff v2.17.0 (https://github.com/JuliaDiff/FiniteDiff.jl)
* NLopt v0.6.5 (https://github.com/JuliaOpt/NLopt.jl)
* SCIP v0.11.10 (https://github.com/scipopt/SCIP.jl): make sure to instal SCIP beforehand.

To make use of these packages, use ```using <name of the package> ```.
## Data
This package contains many data related to pricing problems related to parking choice. We have stored the generated data used for our numerical experiments in the folder #random-beta. Also, you can find the instaces generated to run the algorithm by Li et al. (2019) in the folder #Intel-instance, and the ones for van de Geer and den Boer (2022) in the folder GB. 

In the Data.jl file, we store all the used data in the paper and in this section we discuss how to use them. All the codes return the following outputs:
* Beta_parameter: a matrix whose rows are related to the parking choices (FSP, PSP, and PUP) and columns are related to customer classes. This matrix represents $`\beta^p_{in}`$.
* q_parameter: a matrix whose rows are related to the parking choices (FSP, PSP, and PUP) and columns are related to customer classes. This matrix represents $`q_{in}`$.
* NUM_POINTS: This is the diemension of the problem. 
* N: This is the number of customer classes.
* UB_p: This is the vector containing the upper bounds on the prices
* LB_p: This is the vector containing the lower bounds on the prices.

All the data are store in a function format. So, in order to generate the parameters, you may need to call the function with suitable arguments. 
* ``` Logit_10() ``` returns the tuple (Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p) for MNL model with N=10;
* ``` Logit_50() ``` returns the tuple (Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p) for MNL model with N=50;
* ``` Mixed_Logit_10(β) ``` returns the tuple (Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p) for Mixed Logit model with N=10, given the value β.
* ``` Mixed_Logit_50(β) ``` returns the tuple (Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p) for Mixed Logit model with N=50, given the value β.
* ``` Mixed_Logit_n10_random(R_AT) ``` returns the tuple (Beta_parameter,q_parameter,NUM_POINTS,N,R,UB_p,LB_p) for discrete Mixed Logit model with R_AT number of points with N=10 customer classes. For this function R_AT can be 10, 100, 1000, 10000, 100000, 1000000. 


## Functions and their features
This package contains many functions. In this section we discuss each functions, what they do, and how they can be used:
* ``` discrete_Mixed_logit_function(p) ``` returns the objective value of the Discrete Mixed Logit model at point p. 
* ``` Distribution_function(p) ``` returns the objective value of the MNL model at point p.
* ``` homogenious_MNL(time_limit,p) ``` returns a 2-tuple. Input: time_limit is the limit on the time to run the code, and p is the initial point.  Output: The first element is the objective value and the second element is the obtained solution. This function is the implementation of Algorithm 2 in [Li et al. (2019)](https://pubsonline.informs.org/doi/abs/10.1287/msom.2017.0675). 
* ``` SCIP_degenerate(function_data) ``` returns a 4-tuple. Input: function_data is the data related to a MNL model (for instance Logit_10 in Data.jl). The first element of the output is the time taken by SCIP, the second element is the termination status of SCIP, the third element is the optimality gap, and the last element is the obtained solution by SCIP.
* ``` Mixed_logit_function(p::Vector) ``` returns the objective value of the pricing problem with a continuous logit model. To use this function, the following considerations are needed:
  * this function uses the function ``` Mixed_Logit_distribution(p,β,i) ```. The data of the logit model should be put here. In the first line of this function, the data is acquared (for instance Mixed_Logit_50(β)). Then, the information about the distribution of β should be given. Currently we use the normal distrbution Normal(μ, Σ) using the code ``` pdf(MvNormal(μ, Σ), β) ```.
  * to take the derivative, we need a bound. These are the vectors ``` a ``` (lower bounds on β) and ``` b ``` (upper bound on β) inside the function ``` Mixed_logit_function_i(p, i) ```.
* ``` Mixed_logit_function_discreteDis(p::Vector) ``` returns the value of the discretize relaxation of the continuous mixed logit at point p. For this fucntion the following considerations are needed:
   * all the information about the data is given in the function ``` Mixed_logit_function_i_discreteDis(p::Vector, i) ```. In this function, we currently use the data for ``` Mixed_logit_50 ``` function where the box determined with  ``` a ``` (lower bounds on β) and ``` b ``` (upper bound on β) are discritized by $`R^2`$ points. You need to specify ``` R ``` in this function (default value is 10);
* ``` solve_nlopt(time_limit) ``` uses solver in NLopt package. The current solver is ``` GN_DIRECT_L ``` and the objective function is defined in ``` Mixed_logit_function_nlopt ``` . It returns the information from ``` @timed ``` next to the value of the objective function and obtained solutions. 
* ``` solve_local_opt_Btree(type_of_problem,time_limit,location) ``` returns a 6-tuple. This fucntion is the implimentation of CoBiT in our paper. Input: type_of_problem is string with values either MNL, Mixed, or Discrete, time_limit is the time limit in seconds for CoBiT, and location is a string with the link to the folder where the output is stored. Output: element 1 is the lower bound on the optimal value, element 2 is the best obtain solution, element 3 is the upper bound on the optimal value, element 4 is the list containing the information extracted from CoBiT after each iteration, element 5 is the list containing the break-points, and element 6 is the list containing the number of nodes in each iterations.
   * MNL: uses the fucntions ``` Distribution_function ``` and ``` Distribution_function_i ```
   * Mixed: uses the functions ``` Mixed_logit_function ``` and ``` Mixed_logit_function_i ```
   * Discrete: uses the fucntions ``` Mixed_logit_function_discreteDis ``` and ``` Mixed_logit_function_i_discreteDis ```.
 We emphasize that in definition of the objective function we have $ \sum_{k}\frac{C_k*p+d^k}{f_{k}(p)} $ (see problem (7) in the paper).
## Examples
* parking case with N=10:
```
loc ="C:/Users" #the link to the folder where the output is saved
Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p=Logit_10();  #defining the parameters
R=1; #no probability
##making the C matrix and the d vector
C=zeros(Sum_size,NUM_POINTS);
d=zeros(Sum_size);
for i_p=1:NUM_POINTS
	for i_n=1:N
			C[(i_p-1)*N+i_n,i_p]=1;
	end
end
sol= solve_local_opt_Btree("MNL",7200,loc)
```
* parking case with N=50:
```
loc ="C:/Users" #the link to the folder where the output is saved
Beta_parameter,q_parameter,NUM_POINTS,N,UB_p,LB_p=Logit_50();
Sum_size=NUM_POINTS*N;
C=zeros(Sum_size,NUM_POINTS);
d=zeros(Sum_size);
for i_p=1:NUM_POINTS
	for i_n=1:N
			C[(i_p-1)*N+i_n,i_p]=1;
	end
end
sol= @timed solve_local_opt_Btree("MNL",7200,loc)
```
* solving one of the instances of Intel case
```
instance =1 # the index of the instance

data_inst=JLD2.jldopen("C:/Users/20176914/surfdrive/scientific/University/Codes/Julia/NonlinearProgramming_Virginie/Intel instance/rand_instance_$instance.jld2", "r")
Beta_parameter_old=data_inst["Beta_parameter"];

q_parameter=data_inst["q_parameter"];
NUM_POINTS=data_inst["NUM_POINTS"];
N=data_inst["N"];
R=data_inst["R"];
UB_p_old=data_inst["UB_p"];
LB_p=data_inst["LB_p"];
w_k_old=data_inst["w_k"];
Beta_parameter=copy(Beta_parameter_old);
UB_p=copy(UB_p_old);
w_k=copy(w_k_old);
Sum_size=NUM_POINTS*N*R;
C=zeros(Sum_size,NUM_POINTS);
d=zeros(Sum_size);
close(data_inst)
for i_p=1:NUM_POINTS
	for i_n=1:N
		for r=1:R
			C[(i_p-1)*N*R+(i_n-1)*R+r,i_p]=w_k_old[r]*UB_p_old[i_p];
			w_k[r]=w_k_old[r]*UB_p_old[i_p];
			Beta_parameter[i_p,i_n,r]=Beta_parameter_old[i_p,i_n,r]*UB_p_old[i_p];
		end		
	end
	if UB_p_old[i_p] != 0 
		UB_p[i_p]=1;
	end
end
sol= @timed solve_local_opt_Btree("Discrete",7200,loc)
```
* solving one of the instances of van de Geer and den Boer (2022)
```
using BSON #to read the file
r=0
n_=10
m_=1
d_=BSON.load("GB\\" * string(n_) * string(m_) * string(r)*".bson")  #this should be the folder where the instances are stored
global NUM_POINTS = d_[:n] + 1; #we need to incorporate the no price option
global N= d_[:m];
lb = copy(d_[:p_lb])
global LB_p = pushfirst!(lb ,0);
ub = copy(d_[:p_ub])
global UB_p = pushfirst!(ub ,0);

global Beta_parameter=zeros(NUM_POINTS,N)
beta = copy(d_[:b])
pushfirst!( beta ,0);
for n=1:N
Beta_parameter[:,n] = -  beta
end

global q_parameter=zeros(NUM_POINTS,N);
for n=1:N
q = copy(d_[:a][n])
q_parameter[:,n] = pushfirst!(q ,0)
end
global w = d_[:w]

global Sum_size=NUM_POINTS*N;
global C=zeros(Sum_size,NUM_POINTS);
global d=zeros(Sum_size);
for i_p=1:NUM_POINTS
for i_n=1:N
	C[(i_p-1)*N+i_n,i_p] = w[i_n];
end
end
loc ="C:/Users/20176914/surfdrive/scientific/University/Codes/Julia/NonlinearProgramming_Virginie/MS_pricing paper/CoBIT"
sol= @timed solve_local_opt_Btree("MNL",7200,loc)
```

