// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
// Problems with connectivity: Minimum Cost Steiner Tree
#include <gurobi_c++.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include <string>
#include "myutils.h"
#include <lemon/concepts/digraph.h>
#include <lemon/preflow.h>
using namespace lemon;


// Para compilar:
//  g++ -lemon geompack.cpp myutils.cpp  mygraphlib.cpp generate_steiner_file.cpp -o out



int cutcount = 0;

// Steiner_Instance put all relevant information in one class.
class Steiner_Instance {
public:
  Steiner_Instance(Digraph &graph,
		   DNodeStringMap &vvname,
		   DNodePosMap &posx,
		   DNodePosMap &posy,
		   ArcValueMap &eweight,
		   int nt,
		   vector <DNode> &V); // first nt nodes are terminals
  Digraph &g;
  DNodeStringMap &vname;
  DNodePosMap &px;
  DNodePosMap &py;
  ArcValueMap &weight;
  int nt,nnodes;
  vector <DNode> &V; // Node V[0] is the root. Nodes V[1], ... , V[nt-1] are the destination
};

Steiner_Instance::Steiner_Instance(Digraph &graph,
				   DNodeStringMap &vvname,
				   DNodePosMap &posx,
				   DNodePosMap &posy,
				   ArcValueMap &eweight,
				   int nterm,
				   vector <DNode> &V):
  g(graph), vname(vvname), px(posx), py(posy), weight(eweight), V(V)
{
  nnodes = countNodes(g);
  nt = nterm;
}

// This cutting plane routine inserts finds violated cuts between the root and the other terminals.
// Any cut separating the root from the other terminals must have capacity at least 1
// This is a user cut. That is, it is called when the variables x are still fractionary

class ConnectivityCuts: public GRBCallback
{
  Steiner_Instance &T;
  ListDigraph::ArcMap<GRBVar>& x;
  double (GRBCallback::*solution_value)(GRBVar);
public:
  ConnectivityCuts(Steiner_Instance &T, Digraph::ArcMap<GRBVar>& x) : T(T),x(x)
  {    }
protected:
  void callback()
  {
    if (where==GRB_CB_MIPSOL){ solution_value = &ConnectivityCuts::getSolution;}
    else if (where==GRB_CB_MIPNODE && getIntInfo(GRB_CB_MIPNODE_STATUS)==GRB_OPTIMAL) {
      solution_value = &ConnectivityCuts::getNodeRel;
    } else return;
    try {
      Digraph &g = T.g;
      ArcValueMap capacity(g);
      DCutMap cut(g);
      double vcut;
      for (ArcIt a(g); a!=INVALID; ++a) 
	capacity[a] = (this->*solution_value)(x[a]);  // or getSolution(x[a]);
      
      for (int i=1;i< T.nt;i++) {
	GRBLinExpr expr;
	// find a mincut between root V[0] and other terminal
	vcut = DiMinCut(g,capacity, T.V[0] , T.V[i], cut);
	if (vcut >= 1.0-MY_EPS) continue;

	// found violated cut
	for (ArcIt a(g); a!=INVALID; ++a) 
	  if ((cut[g.source(a)]==cut[T.V[0]]) && (cut[g.target(a)]!=cut[T.V[0]]))
	    expr += x[a];
	addLazy( expr >= 1.0 ); // or addLazy(expr,GRB_GREATER_EQUAL,1.0);
      }
    } catch (GRBException e) {
      cout << "Error number: " << e.getErrorCode() << endl;
      cout << e.getMessage() << endl;
    } catch (...) {
      cout << "Error during callback**" << endl;
    }
  }
};


int ReadListDigraphSteiner(string steiner_filename,
			   ListDigraph &g,
			   DNodeStringMap &vname,
			   ArcValueMap    &weight,
			   DNodePosMap    &px,
			   DNodePosMap    &py,
			   const bool dupla,int &nt,
			   vector <DNode> &V)
{
  string linha,string_number_terminals,digraph_filename;
  ifstream ifSteinerfile(steiner_filename.c_str());
  getline(ifSteinerfile,linha); // 1st line has graphname and no. of terminals
  istringstream iss(linha);
  iss >> digraph_filename;
  iss >> nt;
  ReadListDigraph(digraph_filename,g,vname,weight,px,py,dupla);
  for (int i=0;i<nt;i++) {
    string nodename;
    getline(ifSteinerfile,linha);
    istringstream iss(linha);iss >> nodename;
    for (DNodeIt v(g); v != INVALID; ++v)
      if (vname[v]==nodename) {  V.push_back(v);  break;  }
  }
  ifSteinerfile.close();
  return (1);
}  


int main(int argc, char *argv[]) 
{
  int nt;
  Digraph g;  // graph declaration
  string digraph_steiner_filename;
  DNodeStringMap vname(g);  // name of graph nodes
  DNodePosMap px(g),py(g);  // xy-coodinates for each node
  DNodeColorMap vcolor(g);// color of nodes
  ArcColorMap ecolor(g); // color of edges
  ArcValueMap lpvar(g);    // used to obtain the contents of the LP variables
  ArcValueMap weight(g);   // edge weights
  Digraph::ArcMap<GRBVar> x(g); // binary variables for each arc
  vector <DNode> V;
  int seed=0;
  srand48(1);

  // uncomment one of these lines to change default pdf reader, or insert new one
  //set_pdfreader("open");    // pdf reader for Mac OS X
  //set_pdfreader("xpdf");    // pdf reader for Linux
  //set_pdfreader("evince");  // pdf reader for Linux

    
  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  if (argc!=2) {cout<< endl << "Usage: "<< argv[0]<<"  <digraph_steiner_filename>"<< endl << endl;
    cout << "Examples:      " << argv[0] << " gr_berlin52.steiner" << endl;
    cout << "               " << argv[0] << " gr_usa48.steiner" << endl << endl;
    exit(0);}

  digraph_steiner_filename = argv[1];

  //int time_limit = 3600; // solution must be obtained within time_limit seconds
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
  model.getEnv().set(GRB_IntParam_Seed, seed);
  model.set(GRB_StringAttr_ModelName, "Oriented Steiner Tree with GUROBI"); // prob. name
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a minimization problem

  ReadListDigraphSteiner(digraph_steiner_filename,g,vname,weight,px,py,1,nt,V); 
  Steiner_Instance T(g,vname,px,py,weight,nt,V);
  //for (DNodeIt v(g);v!=INVALID;++v){ if(v==T.V[0])vcolor[v]=RED; else vcolor[v]=BLUE;}
  //for (int i=1;i<T.nt;i++) vcolor[T.V[i]] = MAGENTA;
  //for (ArcIt e(g); e != INVALID; ++e) ecolor[e] = BLUE;
  //ViewListDigraph(g,vname,px,py,vcolor,ecolor,"Triangulated graph");
  
  // Generate the binary variables and the objective function
  // Add one binary variable for each edge and set its cost in the objective function
  for (ArcIt e(g); e != INVALID; ++e) {
    char name[100];
    sprintf(name,"X_%s_%s",vname[g.source(e)].c_str(),vname[g.target(e)].c_str());
    x[e] = model.addVar(0.0, 1.0, weight[e],GRB_BINARY,name); }
  model.update(); // run update to use model inserted variables
  try {
    //if (time_limit >= 0) model.getEnv().set(GRB_DoubleParam_TimeLimit,time_limit);
    //model.getEnv().set(GRB_DoubleParam_ImproveStartTime,10); //try better sol. aft. 10s
    // if (cutoff > 0) model.getEnv().set(GRB_DoubleParam_Cutoff, cutoff );

    ConnectivityCuts cb = ConnectivityCuts(T , x);
    model.setCallback(&cb);
    model.update();
    //model.write("model.lp"); system("cat model.lp");
    model.optimize();

    double soma=0.0;
    for (DNodeIt v(g);v!=INVALID;++v) vcolor[v]=GRAY; // all nodes BLUE
    for (int i=0;i<T.nt;i++) vcolor[T.V[i]]=MAGENTA; // change terminals to MAGENTA
    vcolor[T.V[0]]=RED; // change root to RED
    for (ArcIt e(g); e!=INVALID; ++e) {
      lpvar[e] = x[e].get(GRB_DoubleAttr_X);
      if (BinaryIsOne(lpvar[e])) { soma += weight[e]; ecolor[e] = RED; }
      else ecolor[e] = NOCOLOR; }
    cout << "Steiner Tree Value = " << soma << endl;
    ViewListDigraph(g,vname,px,py,vcolor,ecolor,
	"Steiner Tree cost in graph with "+IntToString(T.nnodes)+
	" nodes and "+IntToString(T.nt)+" terminals: "+DoubleToString(soma));
  } catch (...) {cout << "Error during callback..." << endl; }
  return 0;
}

