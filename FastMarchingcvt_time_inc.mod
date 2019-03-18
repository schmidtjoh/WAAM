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
param v_m; # speed of welding head when moving without welding
param delta_t; #length of one timestep
param a; # thermal diffusity
param b; # radiation factor (eps*sigma*A) 

check: v_w >=0;
check: v_m >=0;
check: kappa_w in [0,1];
check: kappa_e in [0,1];

# NODE DEGREE

set INCIDENCE {v in V} = {(i,j) in W: i==v || j==v};
param degree {v in V} = card(INCIDENCE[v]);

#WORKTIME

set A = {i1 in V, i2 in V: i1!=i2}; 
param dist{(i,j) in A} = ceil(sqrt((X[i]-X[j])**2+(Y[i]-Y[j])**2)*100)/100;
#param dt_A{(i,j) in A} = floor((dist[i,j]/v_m)/delta_t);
param dt_W{(i,j) in W} = ceil((l[i,j]/v_w)/delta_t); #eventuell floor ?
param dt_WW{(i,j) in WW} = if (i,j) in W then dt_W[i,j] else dt_W[j,i];

#NEEDED TIME

param number_of_odd_nodes = card({v in V: degree[v] mod 2 == 1});
param number_of_steps = (number_of_odd_nodes / 2 -1) * 0 + sum {(i,j) in W} dt_W[i,j];

#INDEX SETS TIME

set P={1..number_of_steps};
set P0= {0} union P;

# RADIATION

param intervals;
param bp{i in {1..intervals}, j in {0,1}};
param temp_bp{i in {1..intervals}, j in {0,1}};
param alpha{V,V};
param Pi := 4*atan(1);

# GENERATING MISSING DATA

set WW_Exp = {i in V, pi in P0, j in V, pj in P: (i,j) in WW && pj = pi+dt_WW[i,j]};
set A_Exp = {i in V,j in V, pj in P:(i,j) in A};
 

#VARIABLES

var x {(i,ti,j,tj) in WW_Exp} binary;
var y {(i,j,tj) in A_Exp} binary;


var temp {i in V, p in P0} in [0,phi_w];
var maxtemp in [0,phi_w];

var delta{v in V, p in P0, i in {1..intervals}} in [0,1]; # related to right interval endpoints
var w{v in V, p in P0, i in {1..intervals-1}} binary;

#OBJECTIVE

minimize time: sum{(i,j,tj) in A_Exp} dist[i,j]*y[i,j,tj]+maxtemp;# sum {i in V, p in P} temp[i,p];

#CONSTRAINTS

subject to start_somewhere:
	sum{(i,0,j,tj) in WW_Exp} x[i,0,j,tj] == 1;

subject to end_somewhere:
	sum{(i,ti,j,number_of_steps) in WW_Exp} x[i,ti,j,number_of_steps] == 1;

subject to limit_y:
	sum {(i,j,tj) in A_Exp} y[i,j,tj] == number_of_odd_nodes / 2 - 1;

subject to weld {(i,j) in W}:
	sum {(i,ti,j,tj) in WW_Exp} x[i,ti,j,tj]+sum {(j,tj,i,ti) in WW_Exp} x[j,tj,i,ti] == 1;

subject to path_cont {j in V, p in P: p < number_of_steps}:
	sum {(i,ti,j,p) in WW_Exp} x[i,ti,j,p] + sum {(i,j,p) in A_Exp} y[i,j,p]
	 == sum {(j,p,k,tk) in WW_Exp} x[j,p,k,tk] + sum {(j,k,p) in A_Exp} y[j,k,p];

subject to no_cons_y {p in P: p in (1,number_of_steps)}:
	sum {(i,j,p) in A_Exp} y[i,j,p] <=1;

subject to start_temp {i in V}:
	temp[i,0] == kappa_w * phi_w * sum {(i,0,j,tj) in WW_Exp} x[i,0,j,tj];
	
subject to compute_temp1_lb {i in V, p in P: p<number_of_steps}:#Idee Faktor: kappa_e*a*10**(-6)*dt/(l[j,k]**2)
	temp[i,p] >= (1-kappa_w)*kappa_e*temp[i,p-1] + kappa_w * phi_w 
				- sum {(i,j) in WW} kappa_e*a*exp(-dist[i,j])*(temp[i,p-1]-temp[j,p-1]) # conduction
				- sum {j in V: alpha[i,j] > 0 and i != j} (sum {k in {1..intervals}} delta[i,p-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[i,i]/360)*sin(alpha[i,j]*Pi/180)
				+ sum {j in V: alpha[j,i] > 0 and i != j} (sum {k in {1..intervals}} delta[j,p-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[j,j]/360)*sin(alpha[j,i]*Pi/180)
				- phi_w * (1-(sum {(i,p,j,tj) in WW_Exp} x[i,p,j,tj] + sum {(i,j,p) in A_Exp} y[i,j,p]));	

subject to compute_temp1_ub {i in V, p in P: p<number_of_steps}:
	temp[i,p] <= (1-kappa_w)*kappa_e*temp[i,p-1] + kappa_w * phi_w 
				- sum {(i,j) in WW} kappa_e*a*exp(-dist[i,j])*(temp[i,p-1]-temp[j,p-1])
				- sum {j in V: alpha[i,j] > 0 and i != j} (sum {k in {1..intervals}} delta[i,p-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[i,i]/360)*sin(alpha[i,j]*Pi/180)
				+ sum {j in V: alpha[j,i] > 0 and i != j} (sum {k in {1..intervals}} delta[j,p-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[j,j]/360)*sin(alpha[j,i]*Pi/180)
				+ phi_w * (1-(sum {(i,p,j,tj) in WW_Exp} x[i,p,j,tj] + sum {(i,j,p) in A_Exp} y[i,j,p]));	

subject to compute_temp1_end_lb {j in V}:
	temp[j,number_of_steps] >= (1-kappa_w)*kappa_e*temp[j,number_of_steps-1] + kappa_w * phi_w 
							- sum {(i,j) in WW} kappa_e*a*exp(-dist[i,j])*(temp[i,number_of_steps-1]-temp[j,number_of_steps-1])
							- sum {i in V: alpha[j,i] > 0 and i != j} (sum {k in {1..intervals}} delta[j,number_of_steps-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[j,j]/360)*sin(alpha[j,i]*Pi/180)
							+ sum {i in V: alpha[i,j] > 0 and i != j} (sum {k in {1..intervals}} delta[i,number_of_steps-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[i,i]/360)*sin(alpha[i,j]*Pi/180)
							- phi_w * (1-sum {(i,ti,j,number_of_steps) in WW_Exp} x[i,ti,j,number_of_steps]);	

subject to compute_temp1_end_ub {j in V}:
	temp[j,number_of_steps] <= (1-kappa_w)*kappa_e*temp[j,number_of_steps-1] + kappa_w * phi_w 
							- sum {(i,j) in WW} kappa_e*a*exp(-dist[i,j])*(temp[i,number_of_steps-1]-temp[j,number_of_steps-1])
							- sum {i in V: alpha[j,i] > 0 and i != j} (sum {k in {1..intervals}} delta[j,number_of_steps-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[j,j]/360)*sin(alpha[j,i]*Pi/180)
							+ sum {i in V: alpha[i,j] > 0 and i != j} (sum {k in {1..intervals}} delta[i,number_of_steps-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[i,i]/360)*sin(alpha[i,j]*Pi/180)
							+ phi_w * (1-sum {(i,ti,j,number_of_steps) in WW_Exp} x[i,ti,j,number_of_steps]);	

subject to compute_temp2_lb {i in V, p in P:p<number_of_steps}:
	temp[i,p] >= kappa_e*temp[i,p-1] 
 				- sum {(i,j) in WW} kappa_e*a*exp(-dist[i,j])*(temp[i,p-1]-temp[j,p-1])
				- sum {j in V: alpha[i,j] > 0 and i != j} (sum {k in {1..intervals}} delta[i,p-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[i,i]/360)*sin(alpha[i,j]*Pi/180)
				+ sum {j in V: alpha[j,i] > 0 and i != j} (sum {k in {1..intervals}} delta[j,p-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[j,j]/360)*sin(alpha[j,i]*Pi/180)
				- phi_w * (sum {(i,p,j,tj) in WW_Exp} x[i,p,j,tj] + sum {(i,j,p) in A_Exp} y[i,j,p]);	

subject to compute_temp2_ub {i in V, p in P:p<number_of_steps}:
	temp[i,p] <= kappa_e*temp[i,p-1] 
				- sum {(i,j) in WW} kappa_e*a*exp(-dist[i,j])*(temp[i,p-1]-temp[j,p-1])
				- sum {j in V: alpha[i,j] > 0 and i != j} (sum {k in {1..intervals}} delta[i,p-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[i,i]/360)*sin(alpha[i,j]*Pi/180)
				+ sum {j in V: alpha[j,i] > 0 and i != j} (sum {k in {1..intervals}} delta[j,p-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[j,j]/360)*sin(alpha[j,i]*Pi/180)
				+ phi_w * (sum {(i,p,j,tj) in WW_Exp} x[i,p,j,tj] + sum {(i,j,p) in A_Exp} y[i,j,p]);

subject to compute_temp2_end_lb {j in V}:
	temp[j,number_of_steps] >= kappa_e*temp[j,number_of_steps-1] 
							- sum {(i,j) in WW} kappa_e*a*exp(-dist[i,j])*(temp[i,number_of_steps-1]-temp[j,number_of_steps-1])
							- sum {i in V: alpha[j,i] > 0 and i != j} (sum {k in {1..intervals}} delta[j,number_of_steps-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[j,j]/360)*sin(alpha[j,i]*Pi/180)
							+ sum {i in V: alpha[i,j] > 0 and i != j} (sum {k in {1..intervals}} delta[i,number_of_steps-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[i,i]/360)*sin(alpha[i,j]*Pi/180)
							- phi_w * (sum {(i,ti,j,number_of_steps) in WW_Exp} x[i,ti,j,number_of_steps]);	

subject to compute_temp2_end_ub {j in V}:
	temp[j,number_of_steps] <= kappa_e*temp[j,number_of_steps-1] 
							- sum {(i,j) in WW} kappa_e*a*exp(-dist[i,j])*(temp[i,number_of_steps-1]-temp[j,number_of_steps-1])
							- sum {i in V: alpha[j,i] > 0 and i != j} (sum {k in {1..intervals}} delta[j,number_of_steps-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[j,j]/360)*sin(alpha[j,i]*Pi/180)
							+ sum {i in V: alpha[i,j] > 0 and i != j} (sum {k in {1..intervals}} delta[i,number_of_steps-1,k]*(temp_bp[k,1]-temp_bp[k,0]))*b*(alpha[i,i]/360)*sin(alpha[i,j]*Pi/180)
							+ phi_w * (sum {(i,ti,j,number_of_steps) in WW_Exp} x[i,ti,j,number_of_steps]);

subject to compute_spline1 {i in V, p in P0}:
	temp[i,p] == sum {k in {1..intervals}} delta[i,p,k]*(bp[k,1]-bp[k,0]);
	
subject to comupte_spline2 {i in V, p in P0, k in {1..intervals-1}}:
	w[i,p,k] <= delta[i,p,k];

subject to compute_spline3 {i in V, p in P0, k in {1..intervals-1}}:
	w[i,p,k] >= delta[i,p,k+1];	

subject to compute_maxtemp {i in V, p in P0}:
	temp[i,p] <= maxtemp;
	