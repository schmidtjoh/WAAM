#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 12 09:30:47 2019

@author: jschmidt
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 09:22:22 2019

@author: jschmidt
"""

import numpy as np
import intervals as I
import math

def Intersection(node,a1): # ersetn Schnittpunkt von Richtung a1 und Kanten ermitteln (keine Rückgabe für Nachbarn)
    points = [];    
    for w in W:
        a2 = [V_x[w[1]-1]-V_x[w[0]-1],V_y[w[1]-1]-V_y[w[0]-1]];
        b = np.array([V_x[w[0]-1]-V_x[node-1],V_y[w[0]-1]-V_y[node-1]]);
        A = np.array([[a1[0],-a2[0]],[a1[1],-a2[1]]]);
        if np.linalg.det(A) == 0 and node in w:     # TODO: anderen Knoten der adjazenten Kante ausgeben
            return 0;
        elif np.linalg.det(A) != 0: # singuläre Matrizen von Parallelen ausschließen
            x = np.linalg.solve(A,b);
            if x[0] > 0 and x[1] in I.openclosed(0,1):
                points.append([w[0],w[1],round(x[0],4),round(x[1],4)]);
    temp = [points[p][2] for p in range(0,len(points))];
    if len(temp) != 0:            
        index = temp.index(min(temp));
#        print(points[index])
        if points[index][3] > 0.5:
            return points[index][1];
        else:
            return points[index][0];
    else:
        return 0;

###########################################################################
############################EINGABEPARAMETER###############################
file = '/Users/jschmidt/Documents/Forschung/WAAM/Modell_AMPL/FastMarchingcvt_3gen.dat';
card_V = 8;
card_W = 10;
stepsize = 0.1;
###########################################################################
###########################################################################

V = range(1,card_V+1);
V_x = [];
V_y = [];
W = [];

obj_in = open(file,"r").readlines();
for s in obj_in[1:card_V+1]:
    t = s.split(" ");
    V_x.append(float(t[1].replace('\n','').replace(';','')));
    V_y.append(float(t[2].replace('\n','').replace(';','')));

for s in obj_in[card_V+3:card_V+3+card_W]:
    t = s.split(" ");
    W.append([int(t[0].replace('\n','').replace(';','')),int(t[1].replace('\n','').replace(';',''))]);

alpha = [[0 for v in V] for v in V];            

# Werte für Strahlen berechnen

for v in V:
    for d in np.arange(0,360-stepsize,stepsize):
        direction = np.array([round(math.cos(d/180*math.pi),15),round(math.sin(d/180*math.pi),15)]);
        node2 = Intersection(v,direction)-1;
        if node2 > -1:
            alpha[v-1][node2] += stepsize; 
    

## Formatierung Winkelmatrix als String
#
#z1 = "";    
#z2 = "";
#for v in V:
#    z1 += "   {:d} ".format(v);
#    z2 += "{:d} ".format(v);
#    for k in range(0,len(alpha[0])):
#        z2 += "{:.4f} ".format(alpha[v-1][k]);
#    z2 += "\n ";
#
#
## Ausgabe einer neuen *.dat Datei
#obj_out = open(file.replace('.dat','_rad.dat'),"w");
#obj_out.write("".join(obj_in[l] for l in range(0,card_V+3+card_W+1))+"param intervals := 6;\n"+
#"param bp: 0	 1 :=\n"+
#	"1 	0  200\n"+
#	"2 200  400\n"+
#	"3 400  600\n"+
#	"4 600  800\n"+
#	"5 800 1000\n"+
#	"6 1000 1200;\n\n"+
#"param temp_bp: 0	1 :=\n"+
#	"1 	0  1.6e9\n"+
#	"2 1.6e9  2.56e10\n"+
#	"3 2.56e10  1.296e11\n"+
#	"4 1.296e11  4.096e11\n"+
#	"5 4.096e11 1e12\n"+
#	"6 1e12 2.0736e12;\n\n"+	
#    "param alpha:"+z1+":=\n "+\
#    z2+";\n"+
#    "param phi_w := 1800; # maximal temperature\n\
#param kappa_w := 0.5; # coefficient for blending heat from laser (kappa_w) and existing node heat (1-kappa_w), (>=0 <=1)\n\
#param kappa_e := 0.95; # decay of temperature in node per time step (>=0 <=1)\n\
#param v_w := 0.25; # speed of welding head when welding\n\
#param v_m := 1; #speed of welding head when moving without welding\n\
#param delta_t := 1; # length of one timestep\n\
#param a := 0.1; # thermal diffusity\n\
#param b := 1e-11; # radiation parameter (eps*sigma)");