# SETS

set V;
set W within {V,V};
set WW = {i1 in V, i2 in V: (i1,i2) in W || (i2,i1) in W};


# PARAMETERS

param X {V};
param Y {V};
param l {W};
param phi_w; # maximal temperature
param kappa_w; # coefficient for blending heat from laser (kappa_w) and existing node heat (1-kappa_w), (>=0 <=1)
param kappa_e; # decay of temperature in node per time step (>=0 <=1)
param v_w; #speed of welding head when welding

check: v_w >=0;
check: kappa_w in [0,1];
check: kappa_e in [0,1];

# GENERATING MISSING DATA

set A = {i1 in V, i2 in V: i1!=i2}; 
param dist{(i,j) in A} = ceil(sqrt((X[i]-X[j])**2+(Y[i]-Y[j])**2)*100)/100;
param dt_W{(i,j) in W} = ceil((l[i,j]/v_w));
param dt_WW{(i,j) in WW} = if (i,j) in W then dt_W[i,j] else dt_W[j,i];

# NODE DEGREE

set INCIDENCE {v in V} = {(i,j) in W: i==v || j==v};
param degree {v in V} = card(INCIDENCE[v]);

param number_of_odd_nodes = card({v in V: degree[v] mod 2 == 1});
param number_of_steps = sum {(i,j) in W} dt_W[i,j] + number_of_odd_nodes / 2 -1;

#INDEX SETS TIME

set P={1..number_of_steps};
set P0= P union {0};

#VARIABLES

var x {(i,j) in WW, p in P} binary;
var y {(i,j) in A, p in P} binary;
var x_ind {j in V, p in P0} binary;
var x0 {j in V} binary;
var c {j in V, p in P} binary;
var prod_temp_ind {j in V, p in P} >=0;


var temp {i in V, p in P0} in [0,phi_w];
var maxtemp in [0,phi_w];

#OBJECTIVE

minimize time: sum{(i,j) in A, p in P} dist[i,j]*y[i,j,p]+maxtemp;

#CONSTRAINTS

subject to start_somewhere:
	sum{i in V} x0[i] == 1;

subject to first_decision {i in V}:
	sum{(i,j) in WW} x[i,j,1] == x0[i];

subject to end_somewhere:
	sum{(i,j) in WW} x[i,j,number_of_steps] == 1;

subject to limit_y:
	sum {(i,j) in A, p in P} y[i,j,p] == number_of_odd_nodes / 2 - 1;

subject to weld {(i,j) in W}:
	sum {p in P} (x[i,j,p]+x[j,i,p]) == dt_W[i,j];

subject to more_period_work_t1{(i,j) in WW}:
	sum {k in {1..dt_WW[i,j]}} x[i,j,k] >= dt_WW[i,j] * x[i,j,1];
	
subject to more_period_work {(i,j) in WW, p in P: 1<p<=number_of_steps-dt_WW[i,j]+1}:
	sum {k in {p..p+dt_WW[i,j]-1}} x[i,j,k] >= dt_WW[i,j] * (x[i,j,p]-x[i,j,p-1]);

subject to unique {p in P: p in [1,number_of_steps)}: 
	sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p] == 1;
	
subject to path_cont1 {j in V, p in P: p < number_of_steps}:
	sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p] <= sum {(j,k) in WW} x[j,k,p+1] + sum {(j,k) in A} y[j,k,p+1] + 10 * c[j,p];

subject to path_cont2 {(i,j) in WW, p in P: p < number_of_steps}:
	x[i,j,p] + y[i,j,p] <= x[i,j,p+1] + 10 * (1 - c[j,p]); 

subject to no_cons_y {(k,m) in A, p in P: p in (1,number_of_steps)}:
	sum {(i,j) in A} (y[i,j,p] + y[i,j,p+1]) <=1;

subject to limit_c {p in P}:
	sum {j in V} c[j,p] <= 1;
	
subject to init_indicator {j in V}:
	x_ind[j,0] == x0[j];

subject to set_indicator1 {j in V}:
	sum {(i,j) in WW} (x[i,j,1] - x0[i]) <= x_ind[j,1];

# Indikator beim Verlassen der Kante nicht richtig gesetzt!
subject to set_indicator2 {j in V, p in P: p>1}:
	sum {(i,j) in WW} (x[i,j,p] - x[i,j,p-1])  <= x_ind[j,p];

subject to set_indicator3 {j in V, p in P: p>1}:
	sum {(i,j) in A} y[i,j,p] <= x_ind[j,p];

subject to limit_indicator1 {j in V}:
	sum {p in P} x_ind[j,p] <= ceil(degree[j] / 2);
	
subject to limit_indicator2 {p in P}:	
	sum {j in V} x_ind[j,p] <= 1;

subject to prod_constr1 {j in V, p in P}:
	prod_temp_ind[j,p] <= phi_w * x_ind[j,p];

subject to prod_constr2 {j in V, p in P}:
	prod_temp_ind[j,p] <= temp[j,p-1];

subject to prod_constr3 {j in V, p in P}:
	prod_temp_ind[j,p] >= temp[j,p-1]-phi_w * (1 - x_ind[j,p]);	

subject to start_temp {i in V}:
	temp[i,0] == kappa_w * phi_w * x0[i];

subject to compute_temp1_lb {j in V, p in P}:
	temp[j,p] >= (1-kappa_w)*kappa_e* prod_temp_ind[j,p]+ kappa_w * phi_w * x_ind[j,p] + kappa_e * (temp[j,p-1] - prod_temp_ind[j,p]) - phi_w * (1-(sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p]));

subject to compute_temp1_ub {j in V, p in P}:
	temp[j,p] <= (1-kappa_w)*kappa_e* prod_temp_ind[j,p]+ kappa_w * phi_w * x_ind[j,p] + kappa_e * (temp[j,p-1]- prod_temp_ind[j,p]) + phi_w * (1-(sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p]));


subject to compute_temp2_lb {j in V, p in P}:
	temp[j,p] >= kappa_e*temp[j,p-1] - phi_w * (sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p]);	

subject to compute_temp2_ub {j in V, p in P}:
	temp[j,p] <= kappa_e*temp[j,p-1] + phi_w * (sum {(i,j) in WW} x[i,j,p] + sum {(i,j) in A} y[i,j,p]);	

subject to compute_maxtemp {i in V, p in P0}:
	temp[i,p] <= maxtemp;		
	
	
	