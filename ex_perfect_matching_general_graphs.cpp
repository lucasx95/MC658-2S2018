
//----------------------------------------------------------------------
// Example of an exact program to solve the Minimum Perfect Matching in
// General Graphs, using LEMON Library and Linear Programming
// with GUROBI Solver. The formulation does not use integer variables
// and is an example of integer solution obtained by cutting planes
// obtained in polynomial time.
// The formulation has exponential number of constraints, but
// we insert only a polynomial number of them, using the fact
// that separation of a non-valid point can be made in polynomial time.
//
// Given a graph G=(V,E) with |V| even and edge weights given
// by w:E--> Reals, the formulation is: 
//
// Minimize  Sum_{e in E} w_e * x_e
//   Sum_{e in \delta(v)} x_e == 1     For each node v in V
//   Sum_{e in E[S]} x_e <= (|S|-1)/2  For each S \subseteq V with |S| odd
//   Sum_{e in \delta(S)} x_e >= 1     For each S \subseteq V with |S| odd
//   0 <= x_e <= 1    For each edge e in E,
//
// where \delta(S) is the set of edges with exactly one extremity in S
// and E[S] is the set of edges with two extremities in S
// Denote the first set of constraints by "degree constraints" and the
// second set of constraints by "odd set constraints" (or also known
// as "blossoms constraints")
// 
// Obs.: Analogously (see notes), you can replace the constraints
//   Sum_{e in E[S]} x_e <= (|S|-1)/2  For each S \subseteq V with |S| odd
// by
//   Sum_{e in \delta(S)} x_e >= 1     For each S \subseteq V with |S| odd
//
// Theorem [Edmonds'65]. The polytope given by the above formulation is integral.
//
// [Edmonds'65] J. Edmonds. Maximum matching and a polyhedron with 0,1
//              vertices. J. Res. Nat. Bur. Standards Sect. B 69, 125-130
//              
// To solve the separation problem (find a "blossom constraint" for an
// invalid point, the program uses the Gomory-Hu tree subroutine,
// available from LEMON package.
//
// Theorem. Given a point x satisfying the degree constraints, for the
//          perfect matching formulation, but not all blossoms constraints,
//          then, one of the violated blossom constraints is one of the cuts
//          presented in a Gomory-Hu cut tree.
//
// The proof of this theorem is a good exercise. Using this fact,
// the algorithm coded here is straighforward.
//
// For more details, see the slides in the link (in portuguese) for the formulation
// http://www.ic.unicamp.br/~fkm/lectures/proglin.pdf
// and my course notes (also in portuguese) used in the analysis of algorithms
// course.
//
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------
#include <gurobi_c++.h>
#include <float.h>
#include <math.h>
#include "mygraphlib.h"
#include "myutils.h"
using namespace lemon;

// there are some differences from Gurobi v. 5.5, see where the macro is used
#define GUROBI_NEWVERSION 1  

bool insert_blossom_constraints(GRBModel &model,
				ListGraph &g,
				NodeStringMap &vname,
				NodePosMap &px,  // xy-coodinates for each node
				NodePosMap &py,
				ListGraph::EdgeMap<GRBVar>& x)
{
  EdgeValueMap capacity(g);
  for (EdgeIt e(g); e!=INVALID; ++e) capacity[e] = x[e].get(GRB_DoubleAttr_X);
  GomoryHu<ListGraph, EdgeValueMap > ght(g, capacity);
  ght.run();
  
  // The Gomory-Hu tree is given as a rooted directed tree. Each node has
  // an arc that points to its father. The root node has father -1.
  // Remember that each arc in this tree represents a cut and the value of
  // the arc is the weight of the corresponding cut. So, if an arc has weight
  // less than 1 and each set obtained after removing the arc has odd number
  // of nodes, then we obtain a violated point and we insert the
  // corresponding constraint.
  //
  // To compute the number of nodes for each induced cut (when an edge (u,v) is
  // removed from T), we first make a topological sort of the given tree and
  // we compute, for each node u, the total number of nodes below a node, given
  // by CutSize(u), including the node itself. That is, when the arc u-->predNode(u)
  // is removed, CutSize(u) gives the number of nodes in the component with node u.
  // Given an arc (u-->v), if CutSize[u] is odd, the deletion of arc (u-->v) induces
  // two odd sets.

  //-----------------------------------------------------------------
  // Topological sort of the nodes in T (Gomory-Hu cut tree): arcs have the format: u--->predNode[u]
  NodeIntMap indegree(g);
  NodeIntMap CutSize(g);
  list<Node> ZeroDegreeNodes,topological_sort;
  Node v;
  // Compute the degree of each node in T
  for (NodeIt t(g); t != INVALID; ++t) { indegree[t] = 0; CutSize[t] = 1; } 
  for (NodeIt u(g); u != INVALID; ++u) {
    v = ght.predNode(u);
    if (v==INVALID) continue;// This is the root of T
    indegree[v] = indegree[v]+1; // increase the indegree of the node
  }
  // Get the nodes with indegree zero
  for (NodeIt u(g); u != INVALID; ++u) 
    if (indegree[u]==0) ZeroDegreeNodes.push_back(u);
  
  while (ZeroDegreeNodes.size()>0) {
    Node u = ZeroDegreeNodes.front();// get u, the first element of the list
    ZeroDegreeNodes.pop_front();     // remove u from the list
    topological_sort.push_back(u);   // insert u in the end of the topological sort
    v = ght.predNode(u);
    indegree[v]--;                   // decrement the indegree of v
    if (indegree[v]==0) ZeroDegreeNodes.push_back(v);
  }

  //-----------------------------------------------------------------
  // Compute CutSize[v] for each v
  while (topological_sort.size()>1) { // the last node is the root node 
    Node u = topological_sort.front();// get u, the next node in the topological order
    topological_sort.pop_front();     // remove u 
    v = ght.predNode(u);
    CutSize[v] = CutSize[v] + CutSize[u]; // update the number of nodes below v (the node v itself is already included
  }
  
  //-----------------------------------------------------------------
  // Insert a blossom cut for each violated odd set. Given arc (u,predNode(u)),
  // if CutSize[u] is odd and the corresponding cut has value < 1, we found a violated cut
  bool inserted_new_cut = false;
  for (NodeIt u(g); u != INVALID; ++u) {
    GRBLinExpr expr;
    if ((CutSize[u]%2==0)||(ght.predValue(u)>1.0-MY_EPS)) continue;  // not a violated cut

    for(GomoryHu<ListGraph,EdgeValueMap>::MinCutEdgeIt a(ght,u,ght.predNode(u));a!=INVALID;++a) expr += x[a];
    model.addConstr(expr >= 1.0 );
    inserted_new_cut = true;
  }
  return inserted_new_cut;
}



// This routine also inserts blossom constraints. It is simpler, but less efficient than the above implementation
bool insert_blossom_constraints_slow(GRBModel &model,
				     ListGraph &g,
				     NodeStringMap &vname,
				     NodePosMap &px,  // xy-coodinates for each node
				     NodePosMap &py,  // 
				     ListGraph::EdgeMap<GRBVar>& x)
{
 
  EdgeValueMap capacity(g);
  for (EdgeIt e(g); e!=INVALID; ++e) capacity[e] = x[e].get(GRB_DoubleAttr_X);
  GomoryHu<ListGraph, EdgeValueMap > ght(g, capacity);
  ght.run();
  bool inserted_new_cut = false;
  for (NodeIt u(g); u != INVALID; ++u) {
    if (ght.predNode(u)==INVALID) continue; // is root of Gomory-Hu tree
    double vcut =  ght.predValue(u);
    if (vcut >= 1.0 - MY_EPS) continue;/* is non violated cut */
    int countNodes = 0; // Count the number of nodes in the side of u
    for(GomoryHu<ListGraph,EdgeValueMap>::MinCutNodeIt a(ght,u,ght.predNode(u)); a!=INVALID; ++a) countNodes++;
    if (countNodes%2==1){ // number of nodes in the side of is odd
      GRBLinExpr expr;
      for(GomoryHu<ListGraph,EdgeValueMap>::MinCutEdgeIt e(ght,u,ght.predNode(u), true);e!=INVALID;++e) expr += x[e];
      model.addConstr(expr >= 1 );
      inserted_new_cut = true;
    }
  }
  return inserted_new_cut;
};



int main(int argc, char *argv[]) {
  ListGraph g;  // graph declaration
  string    graph_filename;
  NodeStringMap  vname(g);  // name of graph nodes
  EdgeStringMap  ename(g);  // name of graph nodes
  NodePosMap   px(g), py(g); // xy-coodinates for each node
  NodeColorMap vcolor(g);// color of nodes
  EdgeColorMap ecolor(g); // color of edges
  EdgeValueMap vx(g);    // used to obtain the contents of the LP variables
  EdgeValueMap weight(g);   // edge weights
  srand48(1);

  // uncomment one of these lines to change default pdf reader, or insert new one
  //set_pdfreader("open");    // pdf reader for Mac OS X
  //set_pdfreader("xpdf");    // pdf reader for Linux
  //set_pdfreader("evince");  // pdf reader for Linux

  if (argc!=2) {
    cout<< endl
	<< "This program computes a minimum cost perfect matching in a graph\n" 
	<< "(non-necessarily bipartite). The graph must have even number of nodes.\n\n" 
	<< "Usage: "<< argv[0]<<" <graph_filename>"<<endl << endl <<
      "Example: " << argv[0] << " gr_berlin52" << endl <<
      "         " << argv[0] << " gr_att48" << endl << endl << endl
	<< "Obs.: To generate a euclidean graph, use program generate_random_euclidean_graph.e."
	<< endl; exit(0);}
  else if (!FileExists(argv[1])) {cout<<"File "<<argv[1]<<" does not exist."<<endl; exit(0);}
  
  graph_filename = argv[1];
  ReadListGraph(graph_filename, g, vname, weight, px, py);
  if (countNodes(g) % 2) {
    cout << "\n\nError: Number of nodes in the graph "<< graph_filename
	 << " is odd.\nNumber of nodes must be even." << endl << endl;
    return 0;
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); 

#if GUROBI_NEWVERSION
  model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
#else
  model.getEnv().set(GRB_IntParam_DualReductions, 0); // Dual reductions must be disabled when using lazy constraints
#endif

  
  /* ILP variables */
  ListGraph::EdgeMap<GRBVar> x(g); // variable for connections, 1=connected, 0=not connected

  GRBLinExpr expressao;
  //Cria as variaveis e seta a funcao de otimizacao para a soma de todos os pesos das arestas escolhidas
  for (EdgeIt e(g); e != INVALID; ++e)
    x[e] = model.addVar(0.0, 1.0, weight[e], GRB_CONTINUOUS);
  model.update();

  //restricao de que a somatoria de todas as arestas ligadas a um vertice v tem que ser igual a 1
  for (NodeIt v(g); v != INVALID; ++v) {
    GRBLinExpr expr;
    for (IncEdgeIt e(g, v); e != INVALID; ++e) expr += x[e];
    model.addConstr(expr == 1.0 );
    vcolor[v] = WHITE;
  }
  try {

    model.optimize();    // Insert cutting planes until we cannot separate the node x.
    while (insert_blossom_constraints(model,g,vname,px,py,x)) {
      model.optimize();
    }
  } catch(GRBException e){cerr<<"Could not solve linear program.\n"
			      <<"Code: "<< e.getErrorCode() << " getMessage: "
			      << e.getMessage() << endl;  exit(0);}

  for (EdgeIt e(g); e!=INVALID; ++e) vx[e] = x[e].get(GRB_DoubleAttr_X);

  if (!EdgeVectorIsInteger(g,vx)) {
    cout << "Error: LP solution is not integer.\n";
    for (EdgeIt e(g); e!=INVALID; ++e) {
      ename[e] = vname[g.u(e)]+"_"+vname[g.v(e)];
      cout << "x[" << g.id(g.u(e)) << " , " << g.id(g.v(e)) << "] = " << vx[e] << endl;
    }
    cout << "\nError: LP solution is not integer.\n";
  }
  // else
  // View a graph with edge LP variables x[e] in [0,1] (may have fractional variables)
  ViewEdgeGraphLP(g,vname,px,py,
		  BLACK,   // color of nodes
		  BLUE,    // color of edges with x[e]==1
		  NOCOLOR, // color of edges with x[e]==0
		  RED,     // color of edges with 0 < x[e] < 1
		  vx,      // 0<=vx[e]<=1 for each edge e
		  "Maximum weighted matching in graph with " // message
		  + IntToString(countNodes(g))
		  + " nodes:" + DoubleToString(model.get(GRB_DoubleAttr_ObjVal)));
  
  return 0;
}

