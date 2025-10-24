# GBP-ILP
An Integer Linear Program (ILP) for the Graph Burning Problem (GBP).

The GBP-ILP is implemented in the file GBP-ILP.cpp.

To compile this file you need Gurobi installed in your system and GNU GCC. In particular, we used Gurobi 12.0.3 and GNU GCC 14.0.2.

```
sudo g++ -w -Wall GBP-ILP.cpp -o GBP-ILP -I${GUROBI_HOME}/include -L${GUROBI_HOME}/lib -lgurobi_c++ -lgurobi120
```

To run the resulting executable, you need the following input parameters: input_graph_path number_of_vertices number_of_edges upper_bound.

```
./GBP-ILP /dataset/soc-livejournal 4033137 27933062 15
```
The input graph must be in mtx format. Namely, the first line has the number of vertices, the second line has the number of edges, and the remaining lines have pairs of vertices (edges) separated by a blank space. The folder dataset contains some graphs in this format. There must be exactly one line for each edge.
