// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
// Problems with connectivity: Minimum Cost Steiner Tree
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include <string>
#include <lemon/concepts/digraph.h>
#include <lemon/preflow.h>
using namespace lemon;



int main(int argc, char *argv[]) 
{
  int i,nt,n;
  Digraph g;  // graph declaration
  string graph_filename, terminals_filename;
  DNodeStringMap vname(g);  // name of graph nodes
  DNodePosMap px(g),py(g);  // xy-coodinates for each node
  DNodeColorMap vcolor(g);// color of nodes
  ArcColorMap ecolor(g); // color of edges
  ArcValueMap weight(g);   // edge weights
  srand48(1);

  if (argc!=3) {
    cout<<"Usage: "<< argv[0]<<" <steiner_filename> <number_of_random_terminals>"<<endl;
    exit(0);}
  
  graph_filename = argv[1];
  stringstream(argv[2]) >> nt; // number of terminals

  // Generate a random euclidean graph with n=500 points in the region [0,100)x[0,100)
  // GenerateTriangulatedListDigraph(g,vname,px,py,weight,500,100,100);
  // GenerateRandomEuclideanListDigraph(g,vname,px,py,weight,n,100,100);
  ReadListDigraph(graph_filename,g,vname,weight,px,py,1);
  n = countNodes(g) ;
  vector <DNode> V; // Node V[0] is the root. Nodes V[1], ... , V[nterminals-1] are the destination
  i=0;
  for (DNodeIt v(g); v!=INVALID; ++v) V.push_back(v); // put all nodes of g in V
  // Perform a permutation on the nodes. The first nt nodes are the terminals.
  for (i=0;i<n;i++) {
    DNode aux;
    int j = (int)( n* drand48() );
    aux = V[j]; 
    V[j] = V[i];
    V[i] = aux;
  }
  cout << graph_filename << " " << nt << "\n"; //    graphfilename   number_of_terminals

  for (i=0;i<nt;i++) cout << vname[V[i]] << "\n";
  return 0;
}

