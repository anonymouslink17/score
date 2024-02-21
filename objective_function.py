from __future__ import division
import math

def objective_function(number_of_nodes, number_of_edges_layer_by_layer, beta):
    # maximum density
    maximum_density = 0.0
    # set of layers that maximizes the density
    maximum_layers = set()

    # list of the average degrees
    maximum_average_degrees = [0] * len(number_of_edges_layer_by_layer)

    # sort the layers by their number of edges
    sorted_layers = sorted(number_of_edges_layer_by_layer, key=number_of_edges_layer_by_layer.__getitem__, reverse=True)

    try:
        # for each layer
        for number_of_layers, layer in enumerate(sorted_layers):
            # compute the average degree of the selected layer
            maximum_average_degrees[layer] = number_of_edges_layer_by_layer[layer] / number_of_nodes

            # compute the density
            layer_density = maximum_average_degrees[layer] * (number_of_layers + 1) ** beta

            # if the density is greater or equal than the maximum density
            if layer_density >= maximum_density:
                # update the maximum density
                maximum_density = layer_density
                # update the set of layers that maximizes the density
                maximum_layers = set([inner_layer for inner_layer in sorted_layers[:number_of_layers + 1]])
    except ZeroDivisionError:
            pass

    return maximum_density, maximum_layers, tuple(maximum_average_degrees)










def objective_function_per_node(node_degree, beta):
    # maximum density
    maximum_density = 0.0
    # set of layers that maximizes the density
    maximum_layers = set()

    maximum_density = 0

    # sort the layers by their number of edges
    sorted_layers = sorted(node_degree, key=node_degree.__getitem__, reverse=True)

    try:
        # for each layer
        for number_of_layers, layer in enumerate(sorted_layers):

            # compute the density
            layer_density = node_degree[layer] * (number_of_layers + 1) ** beta

            # if the density is greater or equal than the maximum density
            if layer_density >= maximum_density:
                # update the maximum density
                maximum_density = layer_density
                # update the set of layers that maximizes the density
                maximum_layers = set([inner_layer for inner_layer in sorted_layers[:number_of_layers + 1]])
    except ZeroDivisionError:
            pass

    return maximum_density




def avg_maximin(number_of_nodes, nodes_degrees, beta):
    density = 0
    for node in len(nodes_degrees):
          density += objective_function_per_node(nodes_degrees[node], beta)

    return density/number_of_nodes







def log_objective_function(number_of_nodes, number_of_edges_layer_by_layer, beta):
    # maximum density
    maximum_density = 0.0
    # set of layers that maximizes the density
    maximum_layers = set()

    # list of the average degrees
    maximum_average_degrees = [0] * len(number_of_edges_layer_by_layer)

    # sort the layers by their number of edges
    sorted_layers = sorted(number_of_edges_layer_by_layer, key=number_of_edges_layer_by_layer.__getitem__, reverse=True)

    try:
        # for each layer
        for number_of_layers, layer in enumerate(sorted_layers):
            # compute the average degree of the selected layer
            maximum_average_degrees[layer] = number_of_edges_layer_by_layer[layer] / number_of_nodes

            # compute the density
            layer_density = maximum_average_degrees[layer] * math.log(number_of_layers + 1)

            # if the density is greater or equal than the maximum density
            if layer_density >= maximum_density:
                # update the maximum density
                maximum_density = layer_density
                # update the set of layers that maximizes the density
                maximum_layers = set([inner_layer for inner_layer in sorted_layers[:number_of_layers + 1]])
    except ZeroDivisionError:
            pass

    return maximum_density, maximum_layers, tuple(maximum_average_degrees)