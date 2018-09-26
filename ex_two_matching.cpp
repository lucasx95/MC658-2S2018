// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
// Problems with connectivity: Minimum two-matching
#include <stdio.h>
#include <string>
#include "mygraphlib.h"
#include "myutils.h"
#include <lemon/lp.h>
#include <lemon/list_graph.h>
#include <gurobi_c++.h>
using namespace lemon;
using namespace std;

int main(int argc, char *argv[]) 
{
  ListGraph g;  // graph declaration
  string graph_filename;
  NodeStringMap vname(g);  // name of graph nodes
  EdgeStringMap ename(g);  // name of graph nodes
  NodePosMap px(g),py(g);  // xy-coodinates for each node
  NodeColorMap vcolor(g);// color of nodes
  EdgeColorMap ecolor(g); // color of edges
  EdgeValueMap weight(g);   // edge weights
  srand48(1);

  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  if (argc!=2) {cout<<"Usage: "<< argv[0]<<" <graph_filename>"<< endl;exit(0);}

  graph_filename = argv[1];
  ReadListGraph(graph_filename,g,vname,weight,px,py);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a maximization problem

  /* LPI variables */
  ListGraph::EdgeMap<GRBVar> x(g); // variable for connections, 1=connected, 0=not connected
  
  GRBLinExpr expressao;
  for (EdgeIt e(g); e != INVALID; ++e) {
    // Se grafo for bipartido, basta colocar GRB_CONTINUOUS (integralidade sai por TU).
    x[e] = model.addVar(0.0, 1.0, weight[e], GRB_BINARY); 
  }
  model.update();
  
  for (NodeIt v(g); v!=INVALID; ++v) {
    GRBLinExpr expr;
    int n_edges=0;
    for (IncEdgeIt e(g,v); e != INVALID; ++e) {expr += x[e]; n_edges++;}
    model.addConstr(expr  == 2 );
    vcolor[v] = BLUE;
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  try {
    model.optimize();
    double soma=0.0;
    for (EdgeIt e(g); e!=INVALID; ++e) {
      ename[e] = DoubleToString(weight[e]);
      if (BinaryIsOne(x[e].get(GRB_DoubleAttr_X))) { soma += weight[e]; ecolor[e] = BLUE; }
      else ecolor[e] = NOCOLOR; }
    cout << "Minimum Two-Matching = " << soma << endl;

    // Esta rotina precisa do programa neato/dot do ListGraphviz 
    ViewListGraph(g,vname,ename,px,py,vcolor,ecolor,
    "Minimum weighted two-matching in graph with "+IntToString(countNodes(g))+
    	    " nodes:"+DoubleToString(soma));
  } catch(GRBException e) {
    cerr << "Nao foi possivel resolver o PLI." << endl;
    cerr << "Codigo de erro = " << e.getErrorCode() << endl;
    cerr << e.getMessage();
  }
  return 0;
}

