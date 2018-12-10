# SETS

set V;
set W within {V,V};
set WW = {i1 in V, i2 in V: (i1,i2) in W || (i2,i1) in W};


# PARAMETERS

param X {V};
param Y {V};
param phi_w; # maximal temperature
param kappa_w; # coefficient for blending heat from laser (kappa_w) and existing node heat (1-kappa_w), (>=0 <=1)
param kappa_e; # decay of temperature in node per time step (>=0 <=1)

check: kappa_w in [0,1];
check: kappa_e in [0,1];

# GENERATING MISSING DATA

set A = {i1 in V, i2 in V: i1!=i2}; 
param dist{(i,j) in A} = ceil(sqrt((X[i]-X[j])**2+(Y[i]-Y[j])**2)*100)/100;

# NODE DEGREE

set INCIDENCE {v in V} = {(i,j) in W: i==v || j==v};
param degree {v in V} = card(INCIDENCE[v]);

param number_of_odd_nodes = card({v in V: degree[v] mod 2 == 1});
param number_of_steps = number_of_odd_nodes / 2 + card(W) -1;

set P={1..number_of_steps};
set P0= P union {0};

#VARIABLES

var x {(i,j) in WW, p in P} binary;
var y {(i,j) in A, p in P} binary;

var temp {i in V, p in P0} in [0,phi_w];
var maxtemp in [0,phi_w];

#OBJECTIVE

minimize time: sum{(i,j) in A, p in P} dist[i,j]*y[i,j,p]+maxtemp;# sum {i in V, p in P} temp[i,p];

#CONSTRAINTS

subject to start_somewhere:
	sum{(i,j) in WW} x[i,j,1] == 1;

subject to end_somewhere:
	sum{(i,j) in WW} x[i,j,number_of_steps] == 1;

subject to weld {(i,j) in W}:
	sum {p in P} (x[i,j,p]+x[j,i,p]) == 1;

subject to unique {p in P: p in (1,number_of_steps)}: 
	sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p] == 1;
	
subject to path_cont {j in V, p in P: p < number_of_steps}:
	sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p] == sum {(j,k) in WW} x[j,k,p+1] + sum {(j,k) in A} y[j,k,p+1];

subject to no_cons_y {p in P: p in (1,number_of_steps)}:
	sum {(i,j) in A} y[i,j,p] + sum {(i,j) in A} y[i,j,p+1] <=1;

subject to start_temp {i in V}:
	temp[i,0] == kappa_w * phi_w * sum {(i,j) in WW} x[i,j,1];

subject to compute_temp1_lb {j in V, p in P}:
	temp[j,p] >= (1-kappa_w)*kappa_e*temp[j,p-1] + kappa_w * phi_w - phi_w * (1-(sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p]));

subject to compute_temp1_ub {j in V, p in P}:
	temp[j,p] <= (1-kappa_w)*kappa_e*temp[j,p-1] + kappa_w * phi_w + phi_w * (1-(sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p]));

subject to compute_temp2_lb {j in V, p in P}:
	temp[j,p] >= kappa_e*temp[j,p-1] - phi_w * (sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p]);	

subject to compute_temp2_ub {j in V, p in P}:
	temp[j,p] <= kappa_e*temp[j,p-1] + phi_w * (sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p]);	

subject to compute_maxtemp {i in V, p in P0}:
	temp[i,p] <= maxtemp;		 	
	
	
	
	
	