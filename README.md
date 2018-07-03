
# Info:

Efficient C code to list k-cliques such as described in:  
"Listing k-cliques in Sparse Real-World Graphs. M. Danisch, O. D. Balalau, M. Sozio. WWW2018".  
https://papers-gamma.link/paper/32

Arboricity.c is an implementation of the algorithm to list all k-cliques described in:  
"Arboricity and subgraph listing algorithms. Chiba, N., & Nishizeki, T. SIAM Journal on Computing, 14(1), 210-223, 1985".

# Arboricity.c

This program iterates over all k-cliques. This is an implementation of the 1985 paper of Chiba And Nishizeki "Arboricity and subgraph listing".

## To compile
"gcc Arboricity.c -O9 -o Arboricity".

## To execute
"./Arboricity edgelist.txt k".  
- "edgelist.txt" should contain the graph: one edge on each line separated by a space.
- Will print the number of k-cliques.

# kClist.c

This program iterates over all k-cliques. This is an improvement of the 1985 algorithm of Chiba And Nishizeki detailed in "Arboricity and subgraph listing".

## To compile
"gcc kClist.c -O9 -o kClist".

## To execute
"./kClist k edgelist.txt".  
- "edgelist.txt" should contain the graph: one edge on each line separated by a space.
- Will print the number of k-cliques.

# kClistCore.c

This program enumerates all k-cliques, computes the k-clique core decomposition and computes a k approximation of a k-cliques densest subgraph.

## To compile  
"gcc kClistCore.c -O9 -o kClistCore".

## To execute
"./kClistCore k edgelist.txt ckdeg.txt ckcore.txt ckdens.txt".  
- k of k-clique to enumerate.
- "edgelist.txt" should contain the graph: one edge on each line separated by a space.
- Will print the total number of k-cliques.
- Will print the k-clique core number of the graph
- Will print the k-clique density, the edge density and the size of the found 1/k-approximation densest subgraph
- Will write in ckdeg.txt the ID of each node followed by its k-clique degree
- Will write in ckcore.txt the ID of each node followed by its k-clique core number in a k-clique core ordering
- Will write in ckdens.txt the the k-clique density, the edge density and the size of the densest subgraph followed by the ID of each node in it

# kClistDens.c

This program iterates over all k-cliques. It find a k(1+eps) approximation of the kclique densest. It is highly scallable.

## To compile:
"gcc kClistDens.c -O9 -o kClistDens -fopenmp".

## To execute:
"./kClistDens p k eps edgelist.txt".  
- "edgelist.txt" should contain the graph: one edge on each line separated by a space.
- k is the size of the k-cliques
- p is the number of threads
- eps is the precision
- Will print the number of k-cliques and the density of the found kclique densest.


# kClistEdgeParallel.c

This program iterates over all k-cliques. This is an improvement of the 1985 algorithm of Chiba And Nishizeki detailed in "Arboricity and subgraph listing". It processes edges in parallel.

## To compile:
"gcc kClistEdgeParallel.c -O9 -o kClistEdgeParallel -fopenmp".

## To execute:
"./kClistEdgeParallel p k edgelist.txt".  
- "edgelist.txt" should contain the graph: one edge on each line separated by a space.
- k is the size of the k-cliques
- p is the number of threads
- Will print the number of k-cliques.


# kClistNodeParallel.c

This program iterates over all k-cliques. This is an improvement of the 1985 algorithm of Chiba And Nishizeki detailed in "Arboricity and subgraph listing". It processes nodes in parallel.

## To compile:
"gcc kClistNodeParallel.c -O9 -o kClistNodeParallel -fopenmp".

## To execute:
"./kClistNodeParallel p k edgelist.txt".  
- "edgelist.txt" should contain the graph: one edge on each line separated by a space.
- k is the size of the k-cliques
- p is the number of threads
- Will print the number of k-cliques.

# Initial contributors

Maximilien Danisch  
Technical consultants: Qinna Wang, Oana Balalau and Mauro Sozio  
2015 - 2018  
http://bit.ly/danisch  
maximilien.danisch@gmail.com


