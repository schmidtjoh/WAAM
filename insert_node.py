#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 13:24:56 2018

@author: jschmidt
"""
import math;


def insert_node (zeile, dt):
    n = math.floor(float(zeile[2].replace(";",""))/dt);
    if n > 0:
        array = [zeile[0]+' zw_'+zeile[0]+zeile[1]+str(1)];
        nodes.append('zw_'+zeile[0]+zeile[1]+str(1));
        for i in range(2,n+1):
            array.append('zw_'+zeile[0]+zeile[1]+str(i-1)+' zw_'+zeile[0]+zeile[1]+str(i));
            if ('zw_'+zeile[0]+zeile[1]+str(i)) not in nodes:
                nodes.append('zw_'+zeile[0]+zeile[1]+str(i));
        array.append('zw_'+zeile[0]+zeile[1]+str(n)+' '+str(zeile[1]));
        return array;
    else: 
        return [zeile[0]+' '+zeile[1]];


    

file = '/Users/jschmidt/Documents/Forschung/Modell_AMPL/FastMarchingcvt_grid.dat';
dt = 0.2;
edges = [];
nodes = [];
card_V = 14;
card_W = 19;


obj_in = open(file,"r").readlines();
for s in obj_in[card_V+3:card_V+card_W+3]:
    edges.extend(insert_node(s.split(" "),dt));
    
    
obj_out = open('/Users/jschmidt/Documents/Forschung/Modell_AMPL/FastMarchingcvt_extnodes.dat',"w");
obj_out.write("".join(obj_in[:card_V+2])+\
              "set V_added := "+"".join(str(v)+' ' for v in nodes)+';\n\n'+\
              obj_in[card_V+2]+"".join(str(edges[i])+'\n' for i in range(len(edges)-1))+\
              "".join(str(edges[len(edges)-1])+';\n'+"".join(obj_in[card_V+card_W+3:])));    