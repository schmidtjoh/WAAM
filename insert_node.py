#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 13:24:56 2018

@author: jschmidt
"""
import math;
from decimal import Decimal

def generate_pos (data,i,j,count,n):
    start = data[i-1].split(" ");
    end = data[j-1].split(" ");
    x = Decimal(float(start[1])+(float(end[1])-float(start[1]))*count/(n+1));
    y = Decimal(float(start[2].replace(";",""))+(float(end[2].replace(";",""))-float(start[2].replace(";","")))*count/(n+1));
    return str(round(x,4))+' '+str(round(y,4));


def insert_node (zeile, data, dt):
    n = math.floor(float(zeile[2].replace(";",""))/dt);
    if n > 0:
        array = [zeile[0]+' zw_'+zeile[0]+zeile[1]+str(1)];
        added_nodes.append('zw_'+zeile[0]+zeile[1]+str(1)+' '+generate_pos(data,int(zeile[0]),int(zeile[1]),1,n+1));
        for i in range(2,n+1):
            array.append('zw_'+zeile[0]+zeile[1]+str(i-1)+' zw_'+zeile[0]+zeile[1]+str(i));
            if ('zw_'+zeile[0]+zeile[1]+str(i)) not in added_nodes:
                added_nodes.append('zw_'+zeile[0]+zeile[1]+str(i)+' '+generate_pos(data,int(zeile[0]),int(zeile[1]),i,n+1));
        array.append('zw_'+zeile[0]+zeile[1]+str(n)+' '+str(zeile[1]));
        return array;
    else: 
        return [zeile[0]+' '+zeile[1]];


    

file = '/Users/jschmidt/Documents/Forschung/Modell_AMPL/FastMarchingcvt_time.dat';
dt = 1;
edges = [];
added_nodes = [];
card_V = 14;
card_W = 19;


obj_in = open(file,"r").readlines();
for s in obj_in[card_V+3:card_V+card_W+3]:
    edges.extend(insert_node(s.split(" "),obj_in[1:card_V+1],dt));
    
    
obj_out = open('/Users/jschmidt/Documents/Forschung/Modell_AMPL/Test.dat',"w");
if len(added_nodes) > 0:
    obj_out.write("param: V_ext: X Y a :=\n"+"".join(obj_in[l].replace(";","").replace("\n"," 0 \n") for l in range(1,card_V+1))+\
              "".join(str(added_nodes[v])+' 1 \n' for v in range(len(added_nodes)-1))+\
              "".join(str(added_nodes[len(added_nodes)-1]+' 1;\n\n'))+\
              obj_in[card_V+2]+"".join(str(edges[i])+'\n' for i in range(len(edges)-1))+\
              "".join(str(edges[len(edges)-1])+';\n'+"".join(obj_in[card_V+card_W+3:])));
else:
    obj_out.write("param: V_ext: X Y a :=\n"+"".join(obj_in[l].replace(";","").replace("\n"," 0 \n") for l in range(1,card_V))+\
              "".join(obj_in[card_V].replace(";"," 0;"))+"\n"+"set W :=\n"+"".join(str(edges[i])+'\n' for i in range(len(edges)-1))+\
              "".join(str(edges[len(edges)-1])+';\n'+"".join(obj_in[card_V+card_W+3:])));
                  