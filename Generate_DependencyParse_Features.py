# -*- coding: utf-8 -*-
__author__ = 'Doyel'

import argparse;
import sys;
import json;
reload(sys)
sys.setdefaultencoding("utf-8")
import os;
import copy;
import networkx as nx
from networkx.readwrite import json_graph
import string
from pprint import pprint;
from Utilities import read_config_file, silentremove


# python Generate_DependencyParse_Features.py NonOA_graph_dict.json Dependency_Features.json [-c]
# python Generate_DependencyParse_Features.py OA_graph_dict.json Dependency_Features.json


def main_body():
    global parent_location; global txt_location; global json_location
    parent_location, txt_location, json_location = read_config_file("./config.cfg")

    parser = argparse.ArgumentParser(prog='Generate_DependencyParse_Features', usage='Generate_DependencyParse_Features.py <inputFile> <outputFile> [-c]', description='Script to generate dependency features')
    parser.add_argument('inputJSON', help='Input JSON filename that store the graph, tokens, prot and trigger details')
    parser.add_argument('outputJSON', help='Output JSON filename that will store dependency features for each sent')
    parser.add_argument('-c', action='store_true')  # For a clean build of the output json file

    args = parser.parse_args()

    if args.c:
        print "Performing a clean build of the Dependency features json file"
        silentremove(os.path.join(json_location, args.outputJSON))

    if os.path.isfile(os.path.join(json_location, args.outputJSON)):
        print "Found - " + args.outputJSON + "! Extending the JSON file."
        file_handle = open(os.path.join(json_location, args.outputJSON), "rb")
        mydict = json.load(file_handle)
        file_handle.close()
    else:
        print "The dictionary - " + args.outputJSON + " was not found in the JSON directory. Creating a new one!"
        mydict = {}

    if os.path.isfile(os.path.join(json_location, args.inputJSON)):
        print "Found - " + args.inputJSON + "! Everything is perfect in this world!"
        file_handle = open(os.path.join(json_location, args.inputJSON), "rb")
        networkX_data = json.load(file_handle)  # networkX_data is a dictionary
        file_handle.close()
    else:
        print "The dictionary file - " + args.inputJSON + " was not found in the JSON directory. Exiting!"
        exit()

    for pmid in networkX_data:
        if pmid in mydict:
            print "The PMID already exists in the output dict!", pmid
            print "You may want to perform a clean build of the Dependency Features json file next time"
            exit()
        print "PMID: ", pmid
        mydict[pmid] = {}
        for sentid in networkX_data[pmid]:
            if "Proteins" not in networkX_data[pmid][sentid]:
                continue
            #print "Sentence Id: ", sentid
            mydict[pmid][sentid] = {}
            data = json.loads(networkX_data[pmid][sentid]["Graph"])
            G = json_graph.node_link_graph(data)
            adj_dict = get_adjacency_matrix(G)
            for prot in networkX_data[pmid][sentid]["Proteins"]:    # There may be multiple proteins in a sentence
                mydict[pmid][sentid][prot] = {}
                for tok_id_list in networkX_data[pmid][sentid]["Proteins"][prot]:       # networkX_data[pmid][sentid]["Proteins"][prot] is a LIST of lists, eg. [[7, 8], [52, 53]]
                    for source_id in tok_id_list:  # tok_id_list may be of the form [7, 8] if prot is multi-token
                        shortest_paths = nx.single_source_dijkstra_path(G, source_id, cutoff=None, weight='weight')
                        remove_unwanted_paths(shortest_paths, tok_id_list, networkX_data[pmid][sentid]["Tokens"])
                        features = create_features_edgenames(shortest_paths, adj_dict)
                        features.extend(create_lexi_features_edgenames(G, shortest_paths, adj_dict))
                        populate_mydict(mydict[pmid][sentid][prot], features)
            # pprint(mydict[pmid][sentid]);
            # exit()

    print "Saving the dependency parse features JSON file!"
    json_fn = open(os.path.join(json_location, args.outputJSON), 'w')
    json.dump(mydict, json_fn, indent=4, ensure_ascii=False)
    json_fn.close()


def create_features_edgenames(shortest_paths, adj_dict):
    features = []
    for target in shortest_paths:
        feat_strings = ["PROT"]
        curr_path = shortest_paths[target]  # [0, 1, 2, 3, 4] - 0 is the prot (source) id and 4 is the target id
        if curr_path[-1] != target:
            print "Something went wrong!"; exit()
        edge_list = get_edges(curr_path)    # [(0,1), (1,2), (2,3), (3,4)]
        if path_endsin_punct(edge_list[-1], adj_dict):
            continue
        for edge in edge_list:
            if edge[0] in adj_dict and edge[1] in adj_dict[edge[0]]:
                edge_direction = adj_dict[edge[0]][edge[1]]["direction"]
                edge_name = adj_dict[edge[0]][edge[1]]["name"]
                edge_name = edge_name.replace("_", "")
                temp_direction = str(edge[0]) + ">" + str(edge[1])
                if temp_direction == edge_direction:
                    feat_strings.append("_" + edge_name + ">")
                else:
                    feat_strings.append("<" + edge_name + "_")
            else:
                print "Edge not found in adjacency matrix"
                pprint(adj_dict); exit()
        # print "Feature built: ", '_'.join(feat_strings)
        features.append(''.join(feat_strings))
    return features


def create_lexi_features_edgenames(G, shortest_paths, adj_dict):
    features = []
    for target in shortest_paths:
        feat_strings = []
        curr_path = shortest_paths[target]  # [0, 1, 2, 3, 4] - 0 is the prot (source) id and 4 is the target id
        if curr_path[-1] != target:
            print "Something went wrong!"; exit()
        edge_list = get_edges(curr_path)  # [(0,1), (1,2), (2,3), (3,4)]
        if path_endsin_punct(edge_list[-1], adj_dict):
            continue
        target_word = G.node[target]['word']
        for edge in edge_list:
            if edge[0] in adj_dict and edge[1] in adj_dict[edge[0]]:
                edge_direction = adj_dict[edge[0]][edge[1]]["direction"]
                edge_name = adj_dict[edge[0]][edge[1]]["name"]
                edge_name = edge_name.replace("_", "")
                temp_direction = str(edge[0]) + ">" + str(edge[1])
                if temp_direction == edge_direction:
                    feat_strings.append("_" + edge_name + ">")
                else:
                    feat_strings.append("<" + edge_name + "_")
            else:
                print "Edge not found in adjacency matrix"
                pprint(adj_dict); exit()

        curr_feat = ''.join(feat_strings)
        curr_feat = "PROT" + curr_feat + "*" + target_word
        #print "Feature built: ", curr_feat
        features.append(curr_feat)
    return features


def path_endsin_punct(last_edge, adj_dict):
    if adj_dict[last_edge[0]][last_edge[1]]["name"] == "punct":
        return True
    else:
        return False


def get_adjacency_matrix(G):
    adj_dict = dict()
    for n, nbrsdict in G.adjacency_iter():
        check_node_adjdict(G, adj_dict, n)
        adj_dict[n] = {}
        for nbr, eattr in nbrsdict.iteritems():
            check_node_adjdict(G, adj_dict[n], nbr)
            adj_dict[n][nbr] = copy.deepcopy(eattr)
    # print "Adjacency Matrix: "
    # pprint(adj_dict)
    return adj_dict


def check_node_adjdict(G, mydict, node):
    if node in mydict:
        print "The node is already present in the adjacency matrix!", node
        pprint(mydict)
        pprint(G); nx.draw(G)
        exit()


def get_edges(curr_path):
    edge_list = []
    for i, node in enumerate(curr_path):
        if i < len(curr_path)-1:
            edge_list.append((curr_path[i], curr_path[i+1]))
    # print "Current Path: ", curr_path
    # print "Generated Edge List: ", edge_list
    return edge_list


def remove_unwanted_paths(shortest_paths, tok_id_list, tokens_dict):
    ignore_list = ["DT", "PDT", "WDT", "IN", "CC", "CD", "TO"]
    shortest_paths.pop(0, None) # Remove the path that ends in the "ROOT" token
    for tokid in tok_id_list:
        shortest_paths.pop(tokid, None)
    for target in shortest_paths.keys():
        if tokens_dict[str(target)]["pos"] in ignore_list:
            shortest_paths.pop(target, None)
    # print "Shortest Paths: "
    # pprint(shortest_paths)


def populate_mydict(my_prot_dict, features):
    for feat in features:
        if feat in my_prot_dict:
            my_prot_dict[feat] += 1
        else:
            my_prot_dict[feat] = 1


if __name__ == "__main__":
    main_body()
