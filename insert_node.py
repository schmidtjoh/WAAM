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
    if (float(zeile[2].replace(";",""))/dt) % 1 != 0:
        n = math.floor(float(zeile[2].replace(";",""))/dt);
    else:
        n = math.floor(float(zeile[2].replace(";",""))/dt)-1;
    if n > 0:
        array = [zeile[0]+' '+str(card_V+len(added_nodes)+1)];
        added_nodes.append(str(card_V+len(added_nodes)+1)+' '+generate_pos(data,int(zeile[0]),int(zeile[1]),1,n));
        for i in range(2,n+1):
            array.append(str(card_V+len(added_nodes))+' '+str(card_V+len(added_nodes)+1));
            if str(card_V+len(added_nodes)+1) not in added_nodes:
                added_nodes.append(str(card_V+len(added_nodes)+1)+' '+generate_pos(data,int(zeile[0]),int(zeile[1]),i,n));
        array.append(str(card_V+len(added_nodes))+' '+str(zeile[1]));
        return array;
    else: 
        return [zeile[0]+' '+zeile[1]];


    

file = '/Users/jschmidt/Documents/Forschung/WAAM/Modell_AMPL/FastMarchingcvt_3gen.dat';
dt = 0.1;
edges = [];
added_nodes = [];
card_V = 8;
card_W = 10;


obj_in = open(file,"r").readlines();
for s in obj_in[card_V+3:card_V+card_W+3]:
    edges.extend(insert_node(s.split(" "),obj_in[1:card_V+1],dt));
    
    
obj_out = open('/Users/jschmidt/Documents/Forschung/WAAM/Modell_AMPL/FastMarchingcvt_extnodes'+str(dt)+'_test.dat',"w");
if len(added_nodes) > 0:
    obj_out.write("param: V_ext: X Y z :=\n"+"".join(obj_in[l].replace(";","").replace("\n"," 0 \n") for l in range(1,card_V+1))+\
              "".join(str(added_nodes[v])+' 1 \n' for v in range(len(added_nodes)-1))+\
              "".join(str(added_nodes[len(added_nodes)-1]+' 1;\n\n'))+\
              "set W :=\n"+"".join(str(edges[i])+'\n' for i in range(len(edges)-1))+\
              "".join(str(edges[len(edges)-1])+';\n'+"".join(obj_in[card_V+card_W+3:]))+"\n"+\
              "".join("param dt := "+str(dt)+"; #for data-filename only\n")+\
              "param a := 0.1; # thermal diffusity");
else:
    obj_out.write("param: V_ext: X Y z :=\n"+"".join(obj_in[l].replace(";","").replace("\n"," 0 \n") for l in range(1,card_V))+\
              "".join(obj_in[card_V].replace(";"," 0;"))+"\n"+"set W :=\n"+"".join(str(edges[i])+'\n' for i in range(len(edges)-1))+\
              "".join(str(edges[len(edges)-1])+';\n')+"".join(obj_in[card_V+card_W+3:])+"\n"+\
              "".join("param dt := "+str(dt)+"; #for data-filename only\n")+\
              "param a := 0.1; # thermal diffusity");
                  