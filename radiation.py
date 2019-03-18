#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 09:22:22 2019

@author: jschmidt
"""

# ACHTUNG: Funktioniert nur für konvexe Gebiete und Geraden als Kanten

import numpy as np
import intervals as I
import math

def Visible(node,p1,p2): # Test welcher Knoten welchen Kantenmittelpunkt sieht
    points = [];
    a1 = p2-p1;
    for w in W:
        if node in w: # angrenzende Kanten ausschließen
            points.append([-100,-100]);
        else:    
            a2 = [V_x[w[1]-1]-V_x[w[0]-1],V_y[w[1]-1]-V_y[w[0]-1]];
            b = np.array([V_x[w[0]-1]-p1[0],V_y[w[0]-1]-p1[1]]);
            A = np.array([[a1[0],-a2[0]],[a1[1],-a2[1]]]);
            if np.linalg.det(A) == 0: # singuläre Matrizen von Parallelen ausschließen
                points.append([-1000,-1000]);
            else:
                x = np.linalg.solve(A,b);
                points.append([round(x[0],4),round(x[1],4)]);
    flag = 1;
    for p in points:
        if p[0] in I.closed(0,1) and p[1] in I.closed(0,1) and np.linalg.norm(a1*p[0]) < np.linalg.norm(a1):
            flag = 0;
            break;
    return flag;    

def Cosinussatz(pl,pw,pr): # Berechnung eines Winkels mit Kosinussatz (auch über Skalarprodukt möglich)
    sl = np.linalg.norm(pr-pw);
    sw = np.linalg.norm(pl-pr);
    sr = np.linalg.norm(pl-pw);
    numb = math.acos((sl**2+sr**2-sw**2)/(2*sl*sr));
    return numb/math.pi*180;

###########################################################################
############################EINGABEPARAMETER###############################
file = '/Users/jschmidt/Documents/Forschung/WAAM/Modell_AMPL/horseshoe2.dat';
card_V = 14;
card_W = 23;
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

# Berechnung der Kantenmittelpunkte
MP = [[0]*2]*len(W);
for i in range(0,len(W)):
    MP[i] = [0.5*(V_x[W[i][0]-1]+V_x[W[i][1]-1]),0.5*(V_y[W[i][0]-1]+V_y[W[i][1]-1])];

# Berechnung der Sichtbarkeit        
L=[[0 for mp in MP] for v in V];
for v in V:
    for mp in MP: 
        if v in W[MP.index(mp)]:
            L[v-1][MP.index(mp)] = 0;
        else:    
            L[v-1][MP.index(mp)] = Visible(v,np.array([V_x[v-1],V_y[v-1]]),np.array(mp));

# Berechnung der jeweiligen Winkel
alpha = [[0 for v in V] for v in V];            
for v in V:
    pw = np.array([V_x[v-1],V_y[v-1]]);
    temp = np.where(np.array(L[v-1]) > 0);
    ind = temp[0].tolist();
    for i in ind:
        pr = np.array(MP[i]);
        pl1 = np.array([V_x[W[i][0]-1],V_y[W[i][0]-1]]);
        pl2 = np.array([V_x[W[i][1]-1],V_y[W[i][1]-1]]);
        alpha[v-1][W[i][0]-1] += round(Cosinussatz(pl1,pw,pr),4);
        alpha[v-1][W[i][1]-1] += round(Cosinussatz(pl2,pw,pr),4);
    alpha[v-1][v-1] = round(sum(alpha[v-1]),1);

# Formatierung Winkelmatrix als String

z1 = "";    
z2 = "";
for v in V:
    z1 += "   {:d} ".format(v);
    z2 += "{:d} ".format(v);
    for k in range(0,len(alpha[0])):
        z2 += "{:.4f} ".format(alpha[v-1][k]);
    z2 += "\n ";
# Ausgabe einer neuen *.dat Datei

obj_out = open(file.replace('.dat','_rad.dat'),"w");
obj_out.write("".join(obj_in[l] for l in range(0,card_V+3+card_W+1))+"param intervals := 6;\n"+
"param bp: 0	1 :=\n"+
	"1 	0  200\n"+
	"2 200  400\n"+
	"3 400  600\n"+
	"4 600  800\n"+
	"5 800 1000\n"+
	"6 1000 1200;\n\n"+
"param temp_bp: 0	1 :=\n"+
	"1 	0  1.6e9\n"+
	"2 1.6e9  2.56e10\n"+
	"3 2.56e10  1.296e11\n"+
	"4 1.296e11  4.096e11\n"+
	"5 4.096e11 1e12\n"+
	"6 1e12 2.0736e12;\n\n"+	
    "param alpha:"+z1+":=\n "+\
    z2+";\n"+
    "param phi_w := 1800; # maximal temperature\n\
param kappa_w := 0.5; # coefficient for blending heat from laser (kappa_w) and existing node heat (1-kappa_w), (>=0 <=1)\n\
param kappa_e := 0.95; # decay of temperature in node per time step (>=0 <=1)\n\
param v_w := 0.25; # speed of welding head when welding\n\
param v_m := 1; #speed of welding head when moving without welding\n\
param delta_t := 1; # length of one timestep\n\
param a := 0.1; # thermal diffusity\n\
param b := 1e-11; # radiation parameter (eps*sigma)");