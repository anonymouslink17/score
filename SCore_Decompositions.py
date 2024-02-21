import heapq
import pickle
from math import sqrt
from array import array
from collections import deque, defaultdict
from objective_function import objective_function, avg_maximin, log_objective_function

def FirmCore(multilayer_graph, nodes_iterator, layers_iterator, threshold, information, save=False):
    # degree of each node in each layer
    delta = {}

    # nodes with the same top lambda degree
    delta_map = {}

    # set of neighbors that we need to update
    neighbors = set()

    k_max = 0
    k_start = 0

    if threshold == 1:
        k_start = 1

    information.start_algorithm()

    for node in nodes_iterator:
        delta[node] = [len([neighbor for neighbor in multilayer_graph[node][layer]]) for layer in layers_iterator]
        delta_map[node] = heapq.nlargest(threshold, delta[node])[-1]

    k_max = max(list(delta_map.values()))
    # bin-sort for removing a vertex
    B = [set() for i in range(k_max + 1)]

    for node in nodes_iterator:
        B[delta_map[node]].add(node)

    print("maximum k = ", k_max)
    for k in range(k_start, k_max + 1):
        while B[k]:
            node = B[k].pop()
            delta_map[node] = k
            neighbors = set()

            for layer, layer_neighbors in enumerate(multilayer_graph[node]):
                for neighbor in layer_neighbors:
                    if delta_map[neighbor] > k:
                        delta[neighbor][layer] -= 1
                        if delta[neighbor][layer] + 1 == delta_map[neighbor]:
                            neighbors.add(neighbor)

            for neighbor in neighbors:
                B[delta_map[neighbor]].remove(neighbor)
                delta_map[neighbor] = heapq.nlargest(threshold, delta[neighbor])[-1]
                B[max(delta_map[neighbor], k)].add(neighbor)
        
    # end of the algorithm
    information.end_algorithm(max(delta_map.values()))

    if save:
        a_file = open("../output/" + information.dataset + "_" + str(threshold) + "_FirmCore_decomposition.pkl", "wb")
        pickle.dump(delta_map, a_file)
        a_file.close()

    return




def FirmDcore(multilayer_graph_in, multilayer_graph_out, nodes_iterator, layers_iterator, threshold, information, save=False, count=False):
    # in degree of each node in each layer
    delta_in = {}

    # in degree of each node in each layer
    delta_out = {}

    # nodes with the same top lambda degree
    delta_in_map = {}

    # nodes with the same top lambda degree
    delta_out_map = {}

    
    if threshold == 1:
        k_start = 1
        s_start = 1

    for node in nodes_iterator:
        delta_out[node] = [len([neighbor for neighbor in multilayer_graph_out[node][layer]]) for layer in layers_iterator]
        delta_out_map[node] = heapq.nlargest(threshold, delta_out[node])[-1]

    k_out_max = max(list(delta_out_map.values()))

    information.start_algorithm()
    
    for k in range(k_start, k_out_max):
        # calculate degrees
        for node in nodes_iterator:
            delta_in[node] = [len([neighbor for neighbor in multilayer_graph_in[node][layer]]) for layer in layers_iterator]
            delta_out[node] = [len([neighbor for neighbor in multilayer_graph_out[node][layer]]) for layer in layers_iterator]
            delta_in_map[node] = heapq.nlargest(threshold, delta_in[node])[-1]
            delta_out_map[node] = heapq.nlargest(threshold, delta_out[node])[-1]
        
        s_in_max = max(list(delta_in_map.values()))

        # build Buckets
        B_in = [set() for i in range(s_in_max + 1)]
        for node in nodes_iterator:
            B_in[delta_in_map[node]].add(node)

        
        for s in range(s_start, s_in_max):
            while B_in[s]:
                node = B_in[s].pop()
                delta_in_map[node] = s
                neighbors = set()

                for layer, layer_neighbors in enumerate(multilayer_graph_in[node]):
                    for neighbor in layer_neighbors:
                        if delta_out_map[neighbor] >= k:
                            delta_out[neighbor][layer] -= 1
                            if delta_out[neighbor][layer] + 1 == delta_out_map[neighbor]:
                                neighbors.add(neighbor)
                
                for neighbor in neighbors:
                    neighbors2 = set()
                    delta_out_map[neighbor] = heapq.nlargest(threshold, delta_out[neighbor])[-1]
                    if delta_out_map[neighbor] < k:
                        for layer, layer_neighbors in enumerate(multilayer_graph_out[neighbor]):
                            for neighbor2 in layer_neighbors:
                                if delta_in_map[neighbor2] > s:
                                    delta_in[neighbor2][layer] -= 1
                                    if delta_in[neighbor2][layer] + 1 == delta_in_map[neighbor2]:
                                        neighbors2.add(neighbor2)

                        for neighbor2 in neighbors2:
                            B_in[delta_in_map[neighbor2]].remove(neighbor2)
                            delta_in_map[neighbor2] = heapq.nlargest(threshold, delta_in[neighbor2])[-1]
                            B_in[max(delta_in_map[neighbor2], s)].add(neighbor2)

            # At the end of the for loop, here, we can save all the nodes with delta_out_map == k as S and all the nodes with delta_in_map == s as T
            if save or count:
                exist_tuples_core = []
                core = {'S':[], 'T':[]}
                for node in nodes_iterator:
                    if delta_in_map[node] >= s:
                        core['T'].append(node)
                    if delta_out_map[node] >= k:
                        core['S'].append(node)

                if save:
                    if len(core['S']) * len(core['T']) > 0:
                        a_file = open("../output/" + information.dataset + "_" + str(k) + "_" + str(s) + str(threshold) + "_FirmCore_decomposition.pkl", "wb")
                        pickle.dump(core, a_file)
                        a_file.close()
                
                if count:
                    if len(core['S']) * len(core['T']) > 0:
                        exist_tuples_core.append((k, s))
                
    # end of the algorithm
    information.end_algorithm(0)
    if count:
        with open("../output/" + information.dataset + "directed_none_empty_cores" ".txt", 'w') as f:
            for item in exist_tuples_core:
                f.write("%s " % item)
    return




############################################################
############################################################
############################################################
##################### Decomposition ########################
############################################################
############################################################
############################################################


def FirmCore_decomposition(multilayer_graph, nodes_iterator, layers_iterator, information, save=False):
    for threshold in range(1, max(layers_iterator) + 2):
        print("-------------- threshold = %d --------------"%threshold)
        # degree of each node in each layer
        delta = {}

        # nodes with the same top lambda degree
        delta_map = {}

        # set of neighbors that we need to update
        neighbors = set()

        k_max = 0
        k_start = 0

        # distinct cores
        dist_cores = 0

        if threshold == 1:
            k_start = 1

        information.start_algorithm()

        for node in nodes_iterator:
            delta[node] = [len([neighbor for neighbor in multilayer_graph[node][layer]]) for layer in layers_iterator]
            delta_map[node] = heapq.nlargest(threshold, delta[node])[-1]

        if save:
            a_file = open("../output/" + information.dataset + "_" + str(threshold) + "_FirmCore_decomposition_degree.pkl", "wb")
            pickle.dump(delta_map, a_file)
            a_file.close()

        k_max = max(list(delta_map.values()))
        # bin-sort for removing a vertex
        B = [set() for i in range(k_max + 1)]

        for node in nodes_iterator:
            B[delta_map[node]].add(node)

        print("maximum k = ", k_max)
        for k in range(k_start, k_max + 1):
            if B[k]:
                dist_cores += 1
            while B[k]:
                node = B[k].pop()
                delta_map[node] = k
                neighbors = set()

                for layer, layer_neighbors in enumerate(multilayer_graph[node]):
                    for neighbor in layer_neighbors:
                        if delta_map[neighbor] > k:
                            delta[neighbor][layer] -= 1
                            if delta[neighbor][layer] + 1 == delta_map[neighbor]:
                                neighbors.add(neighbor)

                for neighbor in neighbors:
                    B[delta_map[neighbor]].remove(neighbor)
                    delta_map[neighbor] = heapq.nlargest(threshold, delta[neighbor])[-1]
                    B[max(delta_map[neighbor], k)].add(neighbor)
            
        # end of the algorithm
        information.end_algorithm(max(delta_map.values()))
        information.print_end_algorithm()

        print("Number of Distinct cores = %s"%dist_cores)

        if save:
            a_file = open("../output/" + information.dataset + "_" + str(threshold) + "_FirmCore_decomposition.pkl", "wb")
            pickle.dump(delta_map, a_file)
            a_file.close()





def FirmDcore_decomposition(multilayer_graph_in, multilayer_graph_out, nodes_iterator, layers_iterator, information, save=False, count=False):
    for threshold in range(1, max(layers_iterator) + 2):
        print("-------------- threshold = %d --------------"%threshold)
        # in degree of each node in each layer
        delta_in = {}

        # out degree of each node in each layer
        delta_out = {}

        # nodes with the same top lambda degree
        delta_in_map = {}

        # nodes with the same top lambda degree
        delta_out_map = {}

        # non-epmty cores
        exist_tuples_core = []

        # distinct cores
        dist_cores = 0

        
        if threshold == 1:
            k_start = 1
            s_start = 1

        for node in nodes_iterator:
            delta_out[node] = [len([neighbor for neighbor in multilayer_graph_out[node][layer]]) for layer in layers_iterator]
            delta_out_map[node] = heapq.nlargest(threshold, delta_out[node])[-1]

        k_out_max = max(list(delta_out_map.values()))

        information.start_algorithm()

        
        print("maximum k = ", k_out_max)
        for k in range(k_start, k_out_max):
            flag = False
            active_S_nodes = set(nodes_iterator)
            # calculate degrees
            for node in active_S_nodes.copy():
                delta_out[node] = [len([neighbor for neighbor in multilayer_graph_out[node][layer]]) for layer in layers_iterator]
                delta_out_map[node] = heapq.nlargest(threshold, delta_out[node])[-1]
                if delta_out_map[node] < k:
                    active_S_nodes.remove(node)
               
            
            for node in nodes_iterator:
                delta_in[node] = [len([neighbor for neighbor in multilayer_graph_in[node][layer] if neighbor in active_S_nodes]) for layer in layers_iterator]
                delta_in_map[node] = heapq.nlargest(threshold, delta_in[node])[-1]
            
            s_in_max = max(list(delta_in_map.values()))

            # build Buckets
            B_in = [set() for i in range(s_in_max + 1)]
            for node in nodes_iterator:
                B_in[delta_in_map[node]].add(node)


            for s in range(s_start, s_in_max):
                if B_in[s]:
                    dist_cores += 1
                while B_in[s]:
                    node = B_in[s].pop()
                    delta_in_map[node] = s
                    neighbors = set()

                    for layer, layer_neighbors in enumerate(multilayer_graph_in[node]):
                        for neighbor in layer_neighbors:
                            if delta_out_map[neighbor] >= k:
                                delta_out[neighbor][layer] -= 1
                                if delta_out[neighbor][layer] + 1 == delta_out_map[neighbor]:
                                    neighbors.add(neighbor)
                    
                    for neighbor in neighbors:
                        neighbors2 = set()
                        delta_out_map[neighbor] = heapq.nlargest(threshold, delta_out[neighbor])[-1]
                        if delta_out_map[neighbor] < k:
                            for layer, layer_neighbors in enumerate(multilayer_graph_out[neighbor]):
                                for neighbor2 in layer_neighbors:
                                    if delta_in_map[neighbor2] > s:
                                        delta_in[neighbor2][layer] -= 1
                                        if delta_in[neighbor2][layer] + 1 == delta_in_map[neighbor2]:
                                            neighbors2.add(neighbor2)

                            for neighbor2 in neighbors2:
                                B_in[delta_in_map[neighbor2]].remove(neighbor2)
                                delta_in_map[neighbor2] = heapq.nlargest(threshold, delta_in[neighbor2])[-1]
                                B_in[max(delta_in_map[neighbor2], s)].add(neighbor2)

                # At the end of the for loop, here, we can save all the nodes with delta_out_map == k as S and all the nodes with delta_in_map == s as T
                if save or count:
                    core = {'S':[], 'T':[]}
                    for node in nodes_iterator:
                        if delta_in_map[node] >= s:
                            core['T'].append(node)
                        if delta_out_map[node] >= k:
                            core['S'].append(node)
                        
                    if save:
                        if len(core['S']) * len(core['T']) > 0:
                            flag = True
                            a_file = open("../output/" + information.dataset + "_" + str(k) + "_" + str(s) + "_" + str(threshold) + "_FirmCore_decomposition.pkl", "wb")
                            pickle.dump(core, a_file)
                            a_file.close()
        
                    if count:
                        if len(core['S']) * len(core['T']) > 0:
                            flag = True
                            exist_tuples_core.append((k, s))
            if not flag:
                break

        # end of the algorithm
        print("Number of Distinct cores = %s"%dist_cores)
        information.end_algorithm(0)
        information.print_end_algorithm()
        if count:
            with open("../output/" + information.dataset + "_" + str(threshold) + "_directed_none_empty_cores" ".txt", 'w') as f:
                for item in exist_tuples_core:
                    f.write("%s " % str(item))
















############################################################
############################################################
############################################################
######################## S-core ############################
############################################################
############################################################
############################################################

def identity(X):
    return X

def pure_s_core_decomposition(multilayer_graph, nodes_iterator, vector, layer, information, S=identity, distinct_flag=False, query_nodes=None):

    information.start_algorithm()


    # solution set
    cores = {}
    cores_order = []

    # populate current_k_core with the set of all nodes of the input multilayer_graph
    current_k_core = set(nodes_iterator)
    # instantiate current_vector_list equal to the input vector
    current_vector_list = list(vector)

    # degree of each node in the specified layer
    delta = {}
    # for each node
    for node, neighbors in enumerate(multilayer_graph.adjacency_list):
        # compute the degree in the specified layer
        delta[node] = S(len(neighbors[layer]))

    # sets of nodes divided by degree in the specified layer
    delta_sets = [set() for _ in range(max(delta.itervalues()) + 1)]
    # for each node
    for node, degree in delta.iteritems():
        # put the node in the set corresponding to its degree in the specified layer
        delta_sets[degree].add(node)

    # for each set in delta_sets
    for index, delta_set in enumerate(delta_sets):
        # while the set is not empty
        while len(delta_set) > 0:
            # remove a node from the set and from current_k_core
            node = delta_set.pop()
            current_k_core.remove(node)

            # for each neighbor in the specified layer
            for neighbor in multilayer_graph.adjacency_list[node][layer]:
                # if the neighbor is in current_k_core and its delta is more than index
                if neighbor in current_k_core and delta[neighbor] > index:
                    # update its delta_set
                    delta_sets[delta[neighbor]].remove(neighbor)
        
                    # update its delta
                    delta[neighbor] = S(len(multilayer_graph.adjacency_list[neighbor][layer]))
                    delta_sets[delta[neighbor]].add(neighbor)

        # if the core contains the query nodes or exists
        if (query_nodes is None and len(current_k_core) > 0) or (query_nodes is not None and query_nodes <= current_k_core):
            # build the core index vector of the found core
            current_vector_list[layer] = index + 1
            current_vector = tuple(current_vector_list)
            # add it to the solution set
            cores[current_vector] = array('i', current_k_core)
            cores_order.append(current_vector)
    
        elif query_nodes is not None:
            # otherwise conclude the method
            break
    
    # end of the algorithm
    information.end_algorithm()

    return cores, cores_order






def s_core(multilayer_graph, layers_iterator, vector, ancestors_intersection, information, S=identity, algorithm=None):

    information.start_algorithm()

    # solution set
    k_core = set()

    # queue of inactive nodes that have to be processed
    inactive_nodes = deque()
    # set of active nodes
    active_nodes = set(ancestors_intersection)

    # degree of each node in each layer
    delta = {}

    # for each node in ancestors_intersection
    for node in ancestors_intersection:
        # compute the degree in each layer considering only the set of active nodes
        delta[node] = [len([neighbor for neighbor in multilayer_graph.adjacency_list[node][layer] if neighbor in active_nodes]) for layer in layers_iterator]

        # check if the degree conditions over the layers are satisfied
        in_core = all(S(delta[node])[layer] >= vector[layer] for layer in layers_iterator)

        # if the node is potentially in the core
        if in_core:
            # add the node to the solution set
            k_core.add(node)
        # otherwise
        else:
            # remove the node from the set of active nodes and add it to the queue of inactive nodes
            active_nodes.remove(node)
            inactive_nodes.append(node)

            # while the queue of inactive nodes is not empty
            while len(inactive_nodes) > 0:
                # pop a node
                inactive_node = inactive_nodes.popleft()

                # for each neighbor potentially in the core
                for layer, layer_neighbors in enumerate(multilayer_graph.adjacency_list[inactive_node]):
                    for neighbor in layer_neighbors:
                        if neighbor in k_core:
                            # update its delta
                            delta[neighbor][layer] -= 1

                            # if a degree condition over one layer is no more satisfied
                            if delta[neighbor][layer] < vector[layer]:
                                # remove the node from the solution set and from the set of the active nodes, and add it to the queue of inactive nodes
                                k_core.remove(neighbor)
                                active_nodes.remove(neighbor)
                                inactive_nodes.append(neighbor)

    # if the algorithm is hybrid
    if algorithm == 'h':
        # if the core exists
        if len(k_core) > 0:
            # return the computed core and the vector of the minimum degrees
            return array('i', k_core), tuple([min([delta[node][layer] for node in k_core]) for layer in layers_iterator])

        # otherwise return the computed empty core and vector
        information.end_algorithm()
        return array('i', k_core), vector
    else:
        # otherwise return the computed core
        information.end_algorithm()
        return array('i', k_core) 
    




def s_core_single_vector(multilayer_graph, layers_iterator, vector, layer, nodes, information, query_nodes=None):
    # solution set
    cores = {}
    cores_order = []

    # populate current_k_core with the set of specified nodes
    current_k_core = set(nodes)
    # instantiate current_vector_list equal to the list of the input vector
    current_vector_list = list(vector)

    # degree of each node in each layer
    delta = {}
    # for each node
    for node in nodes:
        # compute the degree in each layer considering only the set of specified nodes
        delta[node] = [len([neighbor for neighbor in multilayer_graph.adjacency_list[node][inner_layer] if neighbor in current_k_core]) for inner_layer in layers_iterator]

    # sets of nodes divided by degree in the specified layer
    delta_sets = [set() for _ in range(max([node_degrees[layer] for node_degrees in delta.itervalues()]) + 1)]
    # for each node
    for node, node_degrees in delta.iteritems():
        # put the node in the set corresponding to its degree in the specified layer
        delta_sets[node_degrees[layer]].add(node)

    # for each set in delta_sets
    for index, delta_set in enumerate(delta_sets):
        # while the set is not empty
        while len(delta_set) > 0:
            # remove a node from the set and from current_k_core
            node = delta_set.pop()
            current_k_core.remove(node)

            # for each neighbor in current_k_core
            for inner_layer, layer_neighbors in enumerate(multilayer_graph.adjacency_list[node]):
                for neighbor in layer_neighbors:
                    if neighbor in current_k_core:
                        delta_neighbor = delta[neighbor][inner_layer]

                        # if inner_layer is the one we are performing the core decomposition on and the delta of the neighbor is more than index
                        if inner_layer == layer and delta_neighbor > index:
                            # update its delta_set
                            delta_sets[delta_neighbor].remove(neighbor)
                            delta_sets[delta_neighbor - 1].add(neighbor)

                            # update its delta
                            delta[neighbor][inner_layer] -= 1

                        # if inner_layer is another layer and the delta of the neighbor is more than the corresponding degree condition
                        elif inner_layer != layer and delta_neighbor >= vector[inner_layer]:
                            # update its delta
                            delta[neighbor][inner_layer] -= 1

                            # if the degree condition over inner_layer is no more satisfied
                            if delta[neighbor][inner_layer] < vector[inner_layer]:
                                # update its delta_set
                                delta_sets[delta[neighbor][layer]].remove(neighbor)
                                delta_sets[index].add(neighbor)

                                # update its delta
                                delta[neighbor][layer] = index

        # if the core contains the query nodes or exists
        if (query_nodes is None and len(current_k_core) > 0) or (query_nodes is not None and query_nodes <= current_k_core):
            # build the core index vector of the found core
            current_vector_list[layer] = index + 1
            current_vector = tuple(current_vector_list)
            # add it to the solution set
            cores[current_vector] = array('i', current_k_core)
            cores_order.append(current_vector)
        elif query_nodes is not None:
            # otherwise conclude the method
            break

    return cores, cores_order




def cmp(a, b):
    return (a > b) - (a < b) 




def build_ancestors_intersection(vector_ancestors, cores, descendants_count, distinct_flag, nodes_iterator=None):
    # order vector_ancestors by the length of the corresponding cores
    vector_ancestors.sort(lambda v, u: cmp(len(cores[v]), len(cores[u])))

    try:
        # initialize ancestors_intersection with the smallest core
        ancestors_intersection = set(cores[vector_ancestors[0]])
        # decrement descendants_count
        decrement_descendants_count(vector_ancestors[0], cores, descendants_count, distinct_flag)

        # for each left ancestor vector
        for ancestor in vector_ancestors[1:]:
            # intersect the corresponding core with ancestors_intersection
            ancestors_intersection &= set(cores[ancestor])
            decrement_descendants_count(ancestor, cores, descendants_count, distinct_flag)

    # if the only ancestor vector is start_vector
    except KeyError:
        ancestors_intersection = set(nodes_iterator)

    return ancestors_intersection


def decrement_descendants_count(vector, cores, descendants_count, distinct_flag):
    # decrement descendants_count
    descendants_count[vector] -= 1

    # if the vector has no more descendants in the queue
    if descendants_count[vector] == 0:
        # delete the core for the map of cores if the distinct_flag is off
        if not distinct_flag:
            del cores[vector]
        # delete vector's entry from descendants_count
        del descendants_count[vector]


def build_descendant_vector(ancestor, index):
    descendant = list(ancestor)
    descendant[index] += 1
    return tuple(descendant)


def build_ancestor_vector(descendant, index):
    ancestor = list(descendant)
    ancestor[index] -= 1
    return tuple(ancestor)


def bottom_up_visit(multilayer_graph, start_vector, end_vector, cores, distinct_flag):
    # initialize the queue of vectors
    vectors_queue = deque()
    vectors_queue.append(start_vector)
    # initialize the set of vectors
    processed_vectors = set()
    processed_vectors.add(start_vector)

    # while vectors_queue is not empty
    while len(vectors_queue) > 0:
        # remove a vector from vectors_queue (FIFO policy)
        vector = vectors_queue.popleft()

        # if the core corresponding to the vector is in not the dict of cores
        if vector not in cores:
            # add the core to the map of cores
            cores[vector] = cores[end_vector]

        # for each layer
        for index in multilayer_graph.layers_iterator:
            # if the ancestor vector is not equal to end_vector
            if vector[index] > end_vector[index] + 1:
                # compute the ancestor vector
                ancestor_vector = build_ancestor_vector(vector, index)

                # if the corresponding core has not been computed yet and it is not in vectors_queue
                if ancestor_vector not in cores and ancestor_vector not in processed_vectors:
                    # add it to the queue
                    vectors_queue.append(ancestor_vector)
                    processed_vectors.add(ancestor_vector)


def post_processing(cores, distinct_flag):
    if distinct_flag:
        # filter distinct cores
        filter_distinct_cores(cores)

    
def filter_distinct_cores(cores):
    # vectors ordered by their level
    ordered_vectors = sorted(cores.iterkeys(), key=sum)

    # for each vector in the ordered list
    for vector in ordered_vectors:
        # build the list of its descendant vectors
        descendant_vectors = [build_descendant_vector(vector, index) for index in range(len(vector))]

        # for each descendant vector
        for descendant_vector in descendant_vectors:
            # if the descendant vector core is equal to the vector core
            if descendant_vector in cores and len(cores[vector]) == len(cores[descendant_vector]) and set(cores[vector]) == set(cores[descendant_vector]):
                # delete the core and break
                del cores[vector]
                break


def filter_inner_most_cores(cores):
    # vectors ordered by their level
    ordered_vectors = sorted(cores.iterkeys(), key=sum)

    # for each vector in the ordered list
    for vector in ordered_vectors:
        # build the list of its descendant vectors
        descendant_vectors = [build_descendant_vector(vector, index) for index in range(len(vector))]

        # for each descendant vector
        for descendant_vector in descendant_vectors:
            # if the descendant vector is in the cores set
            if descendant_vector in cores:
                # delete the core and break
                del cores[vector]
                break























def hybrid(multilayer_graph, distinct_flag):
    # measures
    number_of_cores = 0
    number_of_computed_cores = 0


    # create the vector of zeros from which start the computation
    start_vector = tuple([0] * multilayer_graph.number_of_layers)

    # dict of cores
    cores = {}

    number_of_cores += 1

    # initialize the queue of vectors with the descendants of start_vector and the structure that for each vector saves its ancestor vectors
    vectors_queue = deque()
    ancestors = {}
    for index in multilayer_graph.layers_iterator:
        descendant_vector = build_descendant_vector(start_vector, index)
        vectors_queue.append(descendant_vector)
        ancestors[descendant_vector] = [start_vector]

    # initialize the dictionary that for each vector counts the number of descendants in the queue
    descendants_count = defaultdict(int)

    # compute the core decomposition layer by layer
    for layer in multilayer_graph.layers_iterator:
        cores.update(pure_s_core_decomposition(multilayer_graph, start_vector, layer, distinct_flag=distinct_flag)[0])
    number_of_computed_cores += len(cores) + multilayer_graph.number_of_layers

    # while vectors_queue is not empty
    while len(vectors_queue) > 0:
        # remove a vector from vectors_queue (FIFO policy)
        vector = vectors_queue.popleft()

        # if the number of non zero indexes of vector is equal to the number of its ancestors and is more than 1, build the intersection of its ancestor cores
        number_of_non_zero_indexes = len([index for index in vector if index > 0])
        number_of_ancestors = len(ancestors[vector])
        if number_of_non_zero_indexes == number_of_ancestors and number_of_non_zero_indexes > 1 and vector not in cores:
            ancestors_intersection = build_ancestors_intersection(ancestors[vector], cores, descendants_count, distinct_flag)

            # if the intersection of its ancestor cores is not empty
            if len(ancestors_intersection) > 0:
                # compute the core from it
                k_core, minimum_degrees_vector = s_core(multilayer_graph, multilayer_graph.layers_iterator, vector, ancestors_intersection, algorithm='h')
                number_of_computed_cores += 1

                # if the core is not empty
                if len(k_core) > 0:
                    # add the core to the dict of cores
                    cores[vector] = k_core

                    # if the vector of the minimum degrees is not equal to vector
                    if minimum_degrees_vector != vector:
                        # fill the cores that are equals
                        bottom_up_visit(multilayer_graph, minimum_degrees_vector, vector, cores, distinct_flag)
        else:
            # for each ancestor of vector
            for ancestor in ancestors[vector]:
                # decrement its number of descendants
                decrement_descendants_count(ancestor, cores, descendants_count, distinct_flag)

        # if the core corresponding to the vector is in the dict of cores
        if vector in cores:
            # increment the number of cores
            number_of_cores += 1

            # compute its descendant vectors
            for index in multilayer_graph.layers_iterator:
                descendant_vector = build_descendant_vector(vector, index)

                try:
                    # update the list of the ancestors of the descendant vector
                    ancestors[descendant_vector].append(vector)

                # if the descendant vector has not already been found
                except KeyError:
                    # add the descendant vector to the queue
                    vectors_queue.append(descendant_vector)
                    ancestors[descendant_vector] = [vector]

                # increment descendants_count
                descendants_count[vector] += 1

        # delete vector's entry from ancestors
        del ancestors[vector]

   
    # execute the post processing
    post_processing(cores, distinct_flag)


































############################################################
############################################################
############################################################
################### Densest Subgraph #######################
############################################################
############################################################
############################################################


def Densest_subgraph(multilayer_graph, nodes_iterator, layers_iterator, beta, information):
    # maximum value of objective function given threshold
    max_obj = [0, 0, 0]

    # densest subgraph
    dense_subgraph = 0

    # dense core number
    dense_core_number = (0, 0)

    information.start_algorithm()

    output = set()

    for threshold in range(1, max(layers_iterator) + 2):
        # degree of each node in each layer
        delta = {}

        # nodes with the same top lambda degree
        delta_map = {}

        # set of neighbors that we need to update
        neighbors = set()

        k_max = 0
        k_start = 0

        active_nodes = set(nodes_iterator)

        for node in nodes_iterator:
            delta[node] = [len([neighbor for neighbor in multilayer_graph[node][layer]]) for layer in layers_iterator]
            delta_map[node] = heapq.nlargest(threshold, delta[node])[-1]

        k_max = max(list(delta_map.values()))
        # bin-sort for removing a vertex
        B = [set() for i in range(k_max + 1)]

        for node in nodes_iterator:
            B[delta_map[node]].add(node)
        
        if threshold == 1:
            k_start = 1
            active_nodes.difference_update(list(B[0]))

        for k in range(k_start, k_max + 1):
            temporal_active_nodes = []
            while B[k]:
                node = B[k].pop()
                temporal_active_nodes.append(node)
                delta_map[node] = k
                neighbors = set()

                for layer, layer_neighbors in enumerate(multilayer_graph[node]):
                    for neighbor in layer_neighbors:
                        if delta_map[neighbor] > k:
                            delta[neighbor][layer] -= 1
                            if delta[neighbor][layer] + 1 == delta_map[neighbor]:
                                neighbors.add(neighbor)

                for neighbor in neighbors:
                    B[delta_map[neighbor]].remove(neighbor)
                    delta_map[neighbor] = heapq.nlargest(threshold, delta[neighbor])[-1]
                    B[max(delta_map[neighbor], k)].add(neighbor)
            
            number_nodes = len(active_nodes)
            if number_nodes > 0:
                # compute the number of edges of each layer from delta
                number_of_edges_layer_by_layer = {}
                for layer in layers_iterator:
                    number_of_edges_layer_by_layer[layer] = sum([delta[node][layer] for node in active_nodes]) / 2

                # compute core objective function
                core_objective_function = objective_function(number_nodes, number_of_edges_layer_by_layer, beta)

                if core_objective_function[0] >= max_obj[0]:
                    max_obj = core_objective_function
                    dense_subgraph = number_nodes
                    dense_core_number = (k + 1, threshold)
                    output = set(active_nodes)

            active_nodes.difference_update(temporal_active_nodes)

    # end of the algorithm
    information.end_algorithm(k_max=None)
    information.print_densest_subgraph(beta, max_obj[0], dense_subgraph, max_obj[1], dense_core_number, max_obj[2])

    a_file = open("../output/" + information.dataset + "_" + "densest_subgraph.pkl", "wb")
    pickle.dump(output, a_file)
    a_file.close()


    return 


def Directed_densest_subgraph(multilayer_graph_in, multilayer_graph_out, nodes_iterator, layers_iterator, beta, information):
    # maximum value of objective function given threshold
    max_obj = [0, 0, 0]

    # densest subgraph
    dense_subgraph = 0

    # dense core number
    dense_core_number = (0, 0, 0)

    information.start_algorithm()

    for threshold in range(1, max(layers_iterator) + 2):
        # in degree of each node in each layer
        delta_in = {}

        # in degree of each node in each layer
        delta_out = {}

        # nodes with the same top lambda degree
        delta_in_map = {}

        # nodes with the same top lambda degree
        delta_out_map = {}

        
        if threshold == 1:
            k_start = 1
            s_start = 1

        for node in nodes_iterator:
            delta_out[node] = [len([neighbor for neighbor in multilayer_graph_out[node][layer]]) for layer in layers_iterator]
            delta_out_map[node] = heapq.nlargest(threshold, delta_out[node])[-1]

        k_out_max = max(list(delta_out_map.values()))

        
        for k in range(k_start, k_out_max):
            if k > 1:
                flag = False
            else:
                flag = True
            # calculate degrees
            active_S_nodes = set(nodes_iterator)
            # calculate degrees
            for node in active_S_nodes.copy():
                delta_out[node] = [len([neighbor for neighbor in multilayer_graph_out[node][layer]]) for layer in layers_iterator]
                delta_out_map[node] = heapq.nlargest(threshold, delta_out[node])[-1]
                if delta_out_map[node] < k:
                    active_S_nodes.remove(node)
               
            
            for node in nodes_iterator:
                delta_in[node] = [len([neighbor for neighbor in multilayer_graph_in[node][layer] if neighbor in active_S_nodes]) for layer in layers_iterator]
                delta_in_map[node] = heapq.nlargest(threshold, delta_in[node])[-1]
            
            s_in_max = max(list(delta_in_map.values()))

            # build Buckets
            B_in = [set() for i in range(s_in_max + 1)]
            for node in nodes_iterator:
                B_in[delta_in_map[node]].add(node)
            
            active_nodes_S = set(nodes_iterator)
            active_nodes_T = set(nodes_iterator)

            for s in range(s_start, s_in_max):
                while B_in[s]:
                    node = B_in[s].pop()
                    delta_in_map[node] = s
                    active_nodes_T.remove(node)
                    neighbors = set()

                    for layer, layer_neighbors in enumerate(multilayer_graph_in[node]):
                        for neighbor in layer_neighbors:
                            if delta_out_map[neighbor] >= k:
                                delta_out[neighbor][layer] -= 1
                                if delta_out[neighbor][layer] + 1 == delta_out_map[neighbor]:
                                    neighbors.add(neighbor)
                    
                    for neighbor in neighbors:
                        neighbors2 = set()
                        delta_out_map[neighbor] = heapq.nlargest(threshold, delta_out[neighbor])[-1]
                        if delta_out_map[neighbor] < k:
                            active_nodes_S.remove(neighbor)
                            for layer, layer_neighbors in enumerate(multilayer_graph_out[neighbor]):
                                for neighbor2 in layer_neighbors:
                                    if delta_in_map[neighbor2] > s:
                                        delta_in[neighbor2][layer] -= 1
                                        if delta_in[neighbor2][layer] + 1 == delta_in_map[neighbor2]:
                                            neighbors2.add(neighbor2)

                            for neighbor2 in neighbors2:
                                B_in[delta_in_map[neighbor2]].remove(neighbor2)
                                delta_in_map[neighbor2] = heapq.nlargest(threshold, delta_in[neighbor2])[-1]
                                B_in[max(delta_in_map[neighbor2], s)].add(neighbor2)

                number_nodes = sqrt(len(active_nodes_S) * len(active_nodes_T))
                if number_nodes > 0:
                    flag = True
                    # compute the number of edges of each layer from delta
                    number_of_edges_layer_by_layer = {}
                    for layer in layers_iterator:
                        number_of_edges_layer_by_layer[layer] = sum([delta_out[node][layer] for node in active_nodes_S])

                    # compute core objective function
                    core_objective_function = objective_function(number_nodes, number_of_edges_layer_by_layer, beta)

                    if core_objective_function[0] > max_obj[0]:
                        max_obj = core_objective_function
                        dense_subgraph = (len(active_nodes_S), len(active_nodes_T))
                        dense_core_number = (k, s + 1, threshold)

            if not flag:
                break    

    # end of the algorithm
    information.end_algorithm(0)
    information.print_densest_subgraph(beta, max_obj[0], dense_subgraph, max_obj[1], dense_core_number, max_obj[2])

    return
