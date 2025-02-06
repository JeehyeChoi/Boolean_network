import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import math as mt
import itertools
from concurrent.futures import ProcessPoolExecutor
from HEADER.file_operation import import_truth_table, import_key, import_adj
import networkx as nx
import community as community_louvain


def read_network(fnorigin):
    """파일이 존재하면 내용을 읽고, 없으면 프로그램을 종료하는 함수."""
    file_path = f"./STATE_Networks/{fnorigin}.nw"
    if os.path.exists(file_path):
        state_network = []
        with open(file_path, "r") as r_file:
            for line in r_file:
                ii, matching_row_index = map(int, line.strip().split())
                state_network.append((ii, matching_row_index))
        return state_network
    else:
#        print(f"File {file_path} does not exist. Exiting the program.")
        return None
#        sys.exit(1)  # 프로그램 종료

def directed_state_network(num_genes,adj_header,adj_matrix,fnorigin):

		possible_values = [0, 1]
		state = list(itertools.product(possible_values, repeat=num_genes))
		df_state = pd.DataFrame(state, columns=adj_header)

		state0 = [None] * num_genes
		df_state_fi = pd.DataFrame([state0] * int(mt.pow(2, num_genes)), columns=adj_header)

		with ProcessPoolExecutor() as executor:
			futures = [executor.submit(process_state, ii, row, adj_matrix, adj_header, fnorigin) 
								for ii, row in df_state.iterrows()]
	
		for future in futures:
			ii, df_state_fi_row = future.result()
			df_state_fi.iloc[ii] = df_state_fi_row
			
		return df_state,df_state_fi


def process_state(ii, state_row, adj_matrix, adj_header, fnorigin):
  df_state_fi_row = pd.Series(index=adj_header)
  df_skey = import_key(fnorigin)

  for KEY in df_skey.index:
    non_zero_values = adj_matrix.loc[KEY][adj_matrix.loc[KEY] != 0]

    if len(non_zero_values) > 0:
      condition = []
      for key1 in non_zero_values.index:
        condition.append([key1, state_row[key1]])

      condition_str = ' & '.join(f"{key2} == {value}" for key2, value in condition)

      TT = import_truth_table(fnorigin, KEY)
      f_i = TT.loc[TT.eval(condition_str), TT.columns[-1]]
      df_state_fi_row[KEY] = f_i.iloc[0] if len(f_i) > 0 else None
    else:
      df_state_fi_row[KEY] = state_row[KEY]

  return ii, df_state_fi_row


def find_connected_clusters_and_create_map(G):
    # 모든 연결된 클러스터 찾기
    connected_components = list(nx.connected_components(G.to_undirected()))

    # 클러스터 맵 생성
    cluster_map = {}
    for cluster_id, component in enumerate(connected_components):
        for node in component:
            cluster_map[node] = cluster_id  # 노드와 클러스터 ID 매핑

    return cluster_map

def find_attractors_and_loops(G):
    attractors = []
    loops = []

    for node in G.nodes:
        in_degree = G.in_degree(node)
        out_degree = G.out_degree(node) - G.has_edge(node, node)  # Exclude self-loops

        if out_degree == 0 and in_degree > 0:
            attractors.append(node)

    all_loops = list(nx.simple_cycles(G))
    loops = [loop for loop in all_loops if len(loop) > 1]  # Exclude self-loops
    loops.extend([[attractor] for attractor in attractors])  # Add each attractor as a loop

    return attractors, loops

def basin_nodes_by_cluster(cluster_map):
    # 클러스터 ID를 키로 하여 노드를 그룹화
    cluster_dict = {}
    for node, cluster_id in cluster_map.items():
        if cluster_id not in cluster_dict:
            cluster_dict[cluster_id] = []
        cluster_dict[cluster_id].append(node)

    basin=[]
    for cluster_id, nodes in cluster_dict.items():
        basin.append(len(nodes))
    
    return basin / np.sum(basin)

def create_transition_matrix(G,loops,num_genes,cluster_map):
	# 전이 행렬 초기화
	num_states = len(loops)
	transition_matrix = pd.DataFrame(0, index=range(num_states), columns=range(num_states))
	
	for loop in loops:
		for node in loop:
			i_2=bin(node)[2:].zfill(num_genes)
			for ii in range(0,num_genes):
				binary_list=list(str(i_2))
				binary_list[ii]=str((int(binary_list[ii])+1)%2)
				changed_binary_string = ''.join(binary_list)
				changed_node=int(changed_binary_string, 2)
				transition_matrix.loc[cluster_map[node],cluster_map[changed_node]]+=1
#				print(node,i_2,cluster_map[node],"-->",changed_binary_string,changed_node,cluster_map[changed_node])

	transition_matrix_normalized = transition_matrix.div(transition_matrix.sum(axis=1), axis=0)

	return transition_matrix_normalized, loops




def input_check(adj_header,fnorigin):
	num_genes_input=0
	for key in adj_header:
		if os.path.exists("./DATA_truth/{:s}_TT/{:s}_nospace.csv".format(fnorigin,key)):
			num_genes_input+=1

	return num_genes_input

N=0
if __name__ == "__main__":
	fnO = pd.read_csv("numnodes-filenames", header=None, sep=" ")
	for row in fnO.itertuples(index=False):

		fnorigin = str(row[1])
		adj_header, adj_matrix = import_adj(fnorigin)
		num_genes_input=input_check(adj_header,fnorigin)

		state_network=read_network(fnorigin)
		if state_network!=None:
			G=nx.DiGraph(state_network)
			cluster_map = find_connected_clusters_and_create_map(G)

			df = pd.DataFrame(list(cluster_map.items()), columns=['Node', 'Cluster_ID'])
			df.to_csv("STATE_Networks/{:s}.cluster".format(fnorigin), index=False,sep=" ")

			basin=basin_nodes_by_cluster(cluster_map)
			#print(basin)
		
			attractors,loops= find_attractors_and_loops(G)
			#print("Found loops and attractors:", loops)

			transition_matrix, loop_node = create_transition_matrix(G, loops,num_genes_input,cluster_map)
			transition_matrix.T.to_csv("STATE_Networks/{:s}.trans_matrix".format(fnorigin),index=False,header=False,sep=" ")
			#print(transition_matrix.T)
				
			eigenvalues, eigenvectors = np.linalg.eig(transition_matrix)

			
			largest_eigenvalue = np.max(eigenvalues)
			largest_eigenvalue_index = np.where(eigenvalues == largest_eigenvalue)[0][0]
			largest_eigenvector = eigenvectors[:, largest_eigenvalue_index]
			normalized_vector = (largest_eigenvector / np.sum(largest_eigenvector)) 
			
			diagonal_elements = np.diag(transition_matrix.T)
			print("{:s},{:g},{:g}".format(fnorigin,np.dot(diagonal_elements, normalized_vector.real),np.dot(basin, normalized_vector.real)))
			
			print(eigenvectors.shape)
			print(eigenvectors)
		

			if N==0:
				  break
