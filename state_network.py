import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import sys
import os
import math as mt
import itertools

from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from HEADER.file_operation import import_truth_table, import_key, import_adj

fnO=pd.read_csv("network_title.txt",header=None,sep=" ")
#aa=int(sys.argv[1])
#for i in [aa]:
for i in range(0,len(fnO)):
	fnorigin=str(fnO.iloc[i][0])
	adj_header,adj_matrix=import_adj(fnorigin)
	num_nodes=len(adj_header)
	print(i,fnO.iloc[i][0],num_nodes)

	if(num_nodes<20):
		possible_values = [0,1]
		state=list(itertools.product(possible_values, repeat=num_nodes))
		df_state=pd.DataFrame(columns=adj_header)
		
		#print(len(state))
		state0=[None]*num_nodes
		df_state_fi = pd.concat([pd.DataFrame([state0],columns=adj_header)]*int(mt.pow(2,num_nodes)),ignore_index=True)

		for ii in range(len(state)):
			#print(ii)
			df_state.loc[ii]=state[ii]
			
			df_skey=import_key(fnorigin)
			for KEY in df_skey.index:
				non_zero_values = adj_matrix.loc[KEY][adj_matrix.loc[KEY] != 0]
				
				if len(non_zero_values)>0:
					condition=[]
					for key1 in non_zero_values.index:
						condition.append([key1,df_state.loc[ii,key1]])
						#print(key1,df_state.loc[ii,key1])
					condition_str=' & '.join(f"{key2} == {value}" for key2, value in condition)

					TT=import_truth_table(fnorigin,KEY)
					f_i=TT.loc[TT.eval(condition_str),TT.columns[-1]]
					df_state_fi.loc[ii,KEY]=f_i.iloc[0]
					#print(ii,":",condition_str,"==>>",KEY,"==",f_i.iloc[0])
				else:
					#print(ii,KEY,len(non_zero_values))
					df_state_fi.loc[ii,KEY]=df_state.loc[ii,KEY]

		#print(df_state)
		print(df_state_fi)


		"""
		w1_file=open("./STATE_Networks/{:s}.nw".format(fnorigin),"w")
		for ii in range(len(state)):
			row_fi = df_state_fi.loc[ii]
			matching_row_indices = (df_state == row_fi).all(axis=1)

			if matching_row_indices.any():
				#Using the index of the first match found.
				matching_row_index = matching_row_indices.idxmax()  
				w1_file.write("{:d} {:d}\n".format(ii,matching_row_index))
			else:
				w1_file.write("{:d} {:s}\n".format(ii,"error"))
		
		w1_file.close()
		"""




