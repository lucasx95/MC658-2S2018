// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
// Problems with connectivity: Maximum Bipartite Matching (using digraphs as example,
// arcs have source in one part and target in the other part)
#include <stdio.h>
#include <string>
#include "mygraphlib.h"
#include "myutils.h"
#include <lemon/lp.h>
#include <lemon/list_graph.h>
#include <lemon/concepts/digraph.h>
#include <gurobi_c++.h>
using namespace lemon;
using namespace std;


int main(int argc, char *argv[]) 
{
  Digraph g;  // graph declaration
  string digraph_matching_filename;
  DNodeStringMap vname(g);  // name of graph nodes
  DNodePosMap px(g),py(g);  // xy-coodinates for each node
  DNodeColorMap vcolor(g);// color of nodes
  ArcColorMap ecolor(g); // color of edges
  ArcValueMap lpvar(g);    // used to obtain the contents of the LP variables
  ArcValueMap weight(g);   // edge weights
  srand48(1);


  // uncomment one of these lines to change default pdf reader, or insert new one
  //set_pdfreader("open");    // pdf reader for Mac OS X
  //set_pdfreader("xpdf");    // pdf reader for Linux
  //set_pdfreader("evince");  // pdf reader for Linux

      
  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  if (argc!=2) {
    cout<<endl<<"Usage: "<< argv[0]<<" <digraph_matching_filename>"<<endl<<endl;
    cout << "Example:      " << argv[0] << " digr_bipartite_100_10" << endl << endl;
    exit(0);}

  digraph_matching_filename = argv[1];
  ReadListDigraph(digraph_matching_filename,g,vname,weight,px,py,0);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE); // is a maximization problem

  /* LPI variables */
  Digraph::ArcMap<GRBVar> x(g); // variable for connections, 1=connected, 0=not connected
  
  GRBLinExpr expressao;
  for (ArcIt e(g); e != INVALID; ++e) {
    x[e] = model.addVar(0.0, 1.0, weight[e], GRB_CONTINUOUS);
    // Exercise: Using bipartite graphs, explain why we can use continuous
    // variables and still obtain integer solutions
  }
  model.update();
  
  for (DNodeIt v(g); v!=INVALID; ++v) {
    GRBLinExpr exprin, exprout;
    int n_arcs_in=0,n_arcs_out=0;
    // for each node, the number of arcs leaving is at most 1
    // remember: the graph is bipartite, with arcs going from one part to the other
    for (InArcIt e(g,v); e != INVALID; ++e) {exprin += x[e]; n_arcs_in++;}
    if (n_arcs_in > 0)  {model.addConstr(exprin  <= 1 ); vcolor[v] = BLUE;}
    
    // for each node, the number of arcs entering is at most 1 
    for (OutArcIt e(g,v); e != INVALID; ++e) {exprout += x[e]; n_arcs_out++;}
    if (n_arcs_out > 0) {model.addConstr(exprout <= 1 ); vcolor[v] = RED;}
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  try {
    model.optimize();
    double soma=0.0;
    int cor=0;
    for (ArcIt e(g); e!=INVALID; ++e) {
      lpvar[e] = x[e].get(GRB_DoubleAttr_X);
      if (BinaryIsOne(lpvar[e])) { soma += weight[e]; ecolor[e] = (cor % 8) + 2; cor++; }
      else ecolor[e] = NOCOLOR; }
    cout << "Maximum Bipartite Matching = " << soma << endl;

    // Esta rotina precisa do programa neato/dot do Graphviz 
    ViewListDigraph(g,vname,px,py,vcolor,ecolor,
    "maximum weighted matching in graph with "+IntToString(countNodes(g))+
    	    " nodes:"+DoubleToString(soma));
  } catch(GRBException e) {
    cerr << "Nao foi possivel resolver o PLI." << endl;
    cerr << "Codigo de erro = " << e.getErrorCode() << endl;
    cerr << e.getMessage();
  }
  return 0;
}

