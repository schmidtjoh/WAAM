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
param L{V,V}; # Visibility Matrix
param b; # radiation factor (eps*sigma*A) 

check: v_w >=0;
check: v_m >=0;
check: kappa_w in [0,1];
check: kappa_e in [0,1];

# RADIATION

param intervals;
param bp{1..intervals-1};
param slope{1..intervals};
param alpha{V,V};
param Pi := 4*atan(1);

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

# GENERATING MISSING DATA

set WW_Exp = {i in V, pi in P0, j in V, pj in P: (i,j) in WW && pj = pi+dt_WW[i,j]};
set A_Exp = {i in V, pi in P, j in V, pj in P:(i,j) in A && pj = pi};
 

#VARIABLES

var x {(i,ti,j,tj) in WW_Exp} binary;
var y {(i,ti,j,tj) in A_Exp} binary;


var temp {i in V, p in P0} in [0,phi_w];
var maxtemp in [0,phi_w];

#OBJECTIVE

minimize time: sum{(i,ti,j,tj) in A_Exp} dist[i,j]*y[i,ti,j,tj]+maxtemp;# sum {i in V, p in P} temp[i,p];

#CONSTRAINTS

subject to start_somewhere:
	sum{(i,0,j,tj) in WW_Exp} x[i,0,j,tj] == 1;

subject to end_somewhere:
	sum{(i,ti,j,number_of_steps) in WW_Exp} x[i,ti,j,number_of_steps] == 1;

subject to limit_y:
	sum {(i,ti,j,tj) in A_Exp} y[i,ti,j,tj] == number_of_odd_nodes / 2 - 1;

subject to weld {(i,j) in W}:
	sum {(i,ti,j,tj) in WW_Exp} x[i,ti,j,tj]+sum {(j,tj,i,ti) in WW_Exp} x[j,tj,i,ti] == 1;

subject to path_cont {j in V, p in P: p < number_of_steps}:
	sum {(i,ti,j,p) in WW_Exp} x[i,ti,j,p] + sum {(i,p,j,p) in A_Exp} y[i,p,j,p]
	 == sum {(j,p,k,tk) in WW_Exp} x[j,p,k,tk] + sum {(j,p,k,p) in A_Exp} y[j,p,k,p];

subject to no_cons_y {p in P: p in (1,number_of_steps)}:
	sum {(i,p,j,p) in A_Exp} y[i,p,j,p] <=1;

subject to start_temp {i in V}:
	temp[i,0] == kappa_w * phi_w * sum {(i,0,j,tj) in WW_Exp} x[i,0,j,tj];
	
subject to compute_temp1_lb {i in V, p in P: p<number_of_steps}:#Idee Faktor: kappa_e*a*10**(-6)*dt/(l[j,k]**2)
	temp[i,p] >= (1-kappa_w)*kappa_e*temp[i,p-1] + kappa_w * phi_w 
				- sum {(i,p,j,tj) in WW_Exp} kappa_e*a*(1-dist[i,j])*(temp[i,p-1]-temp[j,p-1]) # conduction
				- sum {j in V: L[i,j] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[i,j])*sin(alpha[i,j]*Pi/180)*temp[i,p-1]
				+ sum {j in V: L[i,j] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[j,i])*sin(alpha[j,i]*Pi/180)*temp[j,p-1]
				- phi_w * (1-(sum {(i,p,j,tj) in WW_Exp} x[i,p,j,tj] + sum {(i,p,j,p) in A_Exp} y[i,p,j,p]));	

subject to compute_temp1_ub {i in V, p in P: p<number_of_steps}:
	temp[i,p] <= (1-kappa_w)*kappa_e*temp[i,p-1] + kappa_w * phi_w 
				- sum {(i,p,j,tj) in WW_Exp} kappa_e*a*(1-dist[i,j])*(temp[i,p-1]-temp[j,p-1])
				- sum {j in V: L[i,j] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[i,j])*sin(alpha[i,j]*Pi/180)*temp[i,p-1]
				+ sum {j in V: L[i,j] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[j,i])*sin(alpha[j,i]*Pi/180)*temp[j,p-1]
				+ phi_w * (1-(sum {(i,p,j,tj) in WW_Exp} x[i,p,j,tj] + sum {(i,p,j,p) in A_Exp} y[i,p,j,p]));	

subject to compute_temp1_end_lb {j in V}:
	temp[j,number_of_steps] >= (1-kappa_w)*kappa_e*temp[j,number_of_steps-1] + kappa_w * phi_w 
							- sum {(i,ti,j,number_of_steps) in WW_Exp} kappa_e*a*(1-dist[i,j])*(temp[i,number_of_steps-1]-temp[j,number_of_steps-1])
							- sum {i in V: L[j,i] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[j,i])*sin(alpha[j,i]*Pi/180)*temp[j,number_of_steps-1]
							+ sum {i in V: L[j,i] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[i,j])*sin(alpha[i,j]*Pi/180)*temp[i,number_of_steps-1]
							- phi_w * (1-sum {(i,ti,j,number_of_steps) in WW_Exp} x[i,ti,j,number_of_steps]);	

subject to compute_temp1_end_ub {j in V}:
	temp[j,number_of_steps] <= (1-kappa_w)*kappa_e*temp[j,number_of_steps-1] + kappa_w * phi_w 
							- sum {(i,ti,j,number_of_steps) in WW_Exp} kappa_e*a*(1-dist[i,j])*(temp[i,number_of_steps-1]-temp[j,number_of_steps-1])
							- sum {i in V: L[j,i] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[j,i])*sin(alpha[j,i]*Pi/180)*temp[j,number_of_steps-1]
							+ sum {i in V: L[j,i] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[i,j])*sin(alpha[i,j]*Pi/180)*temp[i,number_of_steps-1]
							+ phi_w * (1-sum {(i,ti,j,number_of_steps) in WW_Exp} x[i,ti,j,number_of_steps]);	

subject to compute_temp2_lb {i in V, p in P:p<number_of_steps}:
	temp[i,p] >= kappa_e*temp[i,p-1] 
				- sum {(i,p,j,tj) in WW_Exp} kappa_e*a*(1-dist[i,j])*(temp[i,p-1]-temp[j,p-1])
				- sum {j in V: L[i,j] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[i,j])*sin(alpha[i,j]*Pi/180)*temp[i,p-1]
				+ sum {j in V: L[i,j] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[j,i])*sin(alpha[j,i]*Pi/180)*temp[j,p-1]
				- phi_w * (sum {(i,p,j,tj) in WW_Exp} x[i,p,j,tj] + sum {(i,p,j,p) in A_Exp} y[i,p,j,p]);	

subject to compute_temp2_ub {i in V, p in P:p<number_of_steps}:
	temp[i,p] <= kappa_e*temp[i,p-1] 
				- sum {(i,p,j,tj) in WW_Exp} kappa_e*a*(1-dist[i,j])*(temp[i,p-1]-temp[j,p-1])
				- sum {j in V: L[i,j] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[i,j])*sin(alpha[i,j]*Pi/180)*temp[i,p-1]
				+ sum {j in V: L[i,j] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[j,i])*sin(alpha[j,i]*Pi/180)*temp[j,p-1]
				+ phi_w * (sum {(i,p,j,tj) in WW_Exp} x[i,p,j,tj] + sum {(i,p,j,p) in A_Exp} y[i,p,j,p]);

subject to compute_temp2_end_lb {j in V}:
	temp[j,number_of_steps] >= kappa_e*temp[j,number_of_steps-1] 
							- sum {(i,ti,j,number_of_steps) in WW_Exp} kappa_e*a*(1-dist[i,j])*(temp[i,number_of_steps-1]-temp[j,number_of_steps-1])
							- sum {i in V: L[j,i] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[j,i])*sin(alpha[j,i]*Pi/180)*temp[j,number_of_steps-1]
							+ sum {i in V: L[j,i] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[i,j])*sin(alpha[i,j]*Pi/180)*temp[i,number_of_steps-1]
							- phi_w * (sum {(i,ti,j,number_of_steps) in WW_Exp} x[i,ti,j,number_of_steps]);	

subject to compute_temp2_end_ub {j in V}:
	temp[j,number_of_steps] <= kappa_e*temp[j,number_of_steps-1] 
							- sum {(i,ti,j,number_of_steps) in WW_Exp} kappa_e*a*(1-dist[i,j])*(temp[i,number_of_steps-1]-temp[j,number_of_steps-1])
							- sum {i in V: L[j,i] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[j,i])*sin(alpha[j,i]*Pi/180)*temp[j,number_of_steps-1]
							+ sum {i in V: L[j,i] > 0.5} << {v in 1..intervals-1} bp[v];{v in 1..intervals} b*slope[v]>> (1-dist[i,j])*sin(alpha[i,j]*Pi/180)*temp[i,number_of_steps-1]
							+ phi_w * (sum {(i,ti,j,number_of_steps) in WW_Exp} x[i,ti,j,number_of_steps]);

subject to compute_maxtemp {i in V, p in P0}:
	temp[i,p] <= maxtemp;		 	
		