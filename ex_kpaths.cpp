// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
// Problems with connectivity: Minimum Cost k-paths (edge disjoint)
#include <gurobi_c++.h>
#include <iostream>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <queue>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include <string>
#include "myutils.h"
#include <lemon/concepts/digraph.h>
#include <lemon/preflow.h>
using namespace lemon;
using namespace std;


int cutcount = 0;

// kPaths_Instance put all relevant information in one class.
class kPaths_Instance {
public:
  kPaths_Instance(Digraph &graph,
		  DNodeStringMap &vvname,
		  DNodePosMap &posx,
		  DNodePosMap &posy,
		  ArcValueMap &eweight,
		  DNode &sourcenode,
		  DNode &targetnode,
		  int k);
  Digraph &g;
  DNodeStringMap &vname;
  DNodePosMap &px;
  DNodePosMap &py;
  ArcValueMap &weight;
  int k,nnodes;
  DNode &sourcenode;
  DNode &targetnode;
};

kPaths_Instance::kPaths_Instance(Digraph &graph,
				 DNodeStringMap &vvname,
				 DNodePosMap &posx,
				 DNodePosMap &posy,
				 ArcValueMap &eweight,
				 DNode &nsourcenode,
				 DNode &ntargetnode,
				 int nk):
  g(graph), vname(vvname), px(posx), py(posy), weight(eweight),
  sourcenode(nsourcenode), targetnode(ntargetnode)
{
  nnodes = countNodes(g);
  k = nk;
}

int main(int argc, char *argv[]) 
{
  int k,found;
  Digraph g;  // graph declaration
  string digraph_kpaths_filename, source_node_name, target_node_name;
  DNodeStringMap vname(g);  // name of graph nodes
  DNodePosMap px(g),py(g);  // xy-coodinates for each node
  DNodeColorMap vcolor(g);// color of nodes
  ArcColorMap ecolor(g); // color of edges
  ArcValueMap lpvar(g);    // used to obtain the contents of the LP variables
  ArcValueMap weight(g);   // edge weights
  Digraph::ArcMap<GRBVar> x(g); // binary variables for each arc
  vector <DNode> V;
  DNode sourcenode,targetnode;
  int seed=0;
  srand48(1);

  // uncomment one of these lines to change default pdf reader, or insert new one
  //set_pdfreader("open");    // pdf reader for Mac OS X
  //set_pdfreader("xpdf");    // pdf reader for Linux
  set_pdfreader("evince");  // pdf reader for Linux
  //set_pdfreader("open -a Skim.app");
  // double cutoff;   // used to prune non promissing branches (of the B&B tree)
  if (argc!=5) {cout<<endl<<"Usage: "<< argv[0]<<"  <digraph_kpaths_filename>  <source_node_name>  <target_node_name>  <k>"<< endl << endl;
    cout << "Example:      " << argv[0] << " digr_triang_sparse_100 12 50 5" << endl << endl;
    exit(0);}

  digraph_kpaths_filename = argv[1];
  source_node_name = argv[2];
  target_node_name = argv[3];
  k = atoi(argv[4]);

  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.getEnv().set(GRB_IntParam_Seed, seed);
  model.set(GRB_StringAttr_ModelName, "Oriented k-Paths with GUROBI"); // prob. name
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a minimization problem

  ReadListDigraph(digraph_kpaths_filename,g,vname,weight,px,py,0);
  found=0;
  for (DNodeIt v(g);v!=INVALID;++v)
    if(vname[v]==source_node_name){sourcenode=v;found=1;break;}
  if (!found) {cout<<"Could not find source node "<<source_node_name<<endl;exit(0);}
  found=0;
  for (DNodeIt v(g);v!=INVALID;++v)
    if(vname[v]==target_node_name){targetnode=v;found=1;break;}
  if (!found) {cout<<"Could not find target node "<<target_node_name<<endl;exit(0);}
    
  kPaths_Instance T(g,vname,px,py,weight,sourcenode,targetnode,k);
  
  //for (DNodeIt v(g);v!=INVALID;++v){ if(v==T.V[0])vcolor[v]=RED; else vcolor[v]=BLUE;}
  //for (int i=1;i<T.nt;i++) vcolor[T.V[i]] = MAGENTA;
  //for (ArcIt e(g); e != INVALID; ++e) ecolor[e] = BLUE;
  //ViewListDigraph(g,vname,px,py,vcolor,ecolor,"Triangulated graph");
  
  // Generate the binary variables and the objective function
  // Add one binary variable for each edge and set its cost in the objective function
  for (ArcIt e(g); e != INVALID; ++e) {
    char name[100];
    sprintf(name,"X_%s_%s",vname[g.source(e)].c_str(),vname[g.target(e)].c_str());
    x[e] = model.addVar(0.0, 1.0, weight[e],GRB_CONTINUOUS,name); }
  model.update(); // run update to use model inserted variables
  try {
    //if (time_limit >= 0) model.getEnv().set(GRB_DoubleParam_TimeLimit,time_limit);
    //model.getEnv().set(GRB_DoubleParam_ImproveStartTime,10); //try better sol. aft. 10s
    // if (cutoff > 0) model.getEnv().set(GRB_DoubleParam_Cutoff, cutoff );
    //model.write("model.lp"); system("cat model.lp");

    // Add degree constraint for each node (sum of solution edges incident to a node is 2)
    for (DNodeIt v(g); v!=INVALID; ++v) {
      GRBLinExpr exprin, exprout;
      for (InArcIt e(g,v); e != INVALID; ++e) exprin += x[e];
      for (OutArcIt e(g,v); e != INVALID; ++e) exprout += x[e];

      if (v==sourcenode)      model.addConstr(exprout - exprin == k );
      else if (v==targetnode) model.addConstr(exprin - exprout == k );
      else                model.addConstr(exprin - exprout == 0 );
    }

    model.optimize();

    double soma=0.0;
    for (DNodeIt v(g);v!=INVALID;++v) vcolor[v]=BLUE; // all nodes BLUE
    vcolor[sourcenode]=RED; // change the colors of the source node
    vcolor[targetnode]=RED; // and the target node to RED
    for (ArcIt e(g); e!=INVALID; ++e) {
      if (BinaryIsOne(x[e].get(GRB_DoubleAttr_X))) {
	soma += weight[e];
	ecolor[e] = RED;
      }else{
	ecolor[e] = GRAY;
      }
    }
    cout << "kPaths Tree Value = " << soma << endl;
    /*
    //------------------------------------------------------------------------------
    // The following code was implemented by Lucas Prado Melo: coloring paths by bfs
    bool ok = true;
    for(int i=0; i < k && ok; i++) {
      queue<DNode> q;
      DNodeArcMap pre(g, INVALID);
      q.push(sourcenode);
      while (!q.empty()) {
        DNode cur = q.front();
        q.pop();
        for(OutArcIt e(g, cur); e!=INVALID; ++e) {
          DNode nxt = g.runningNode(e);
          if ((pre[nxt]==INVALID) && (ecolor[e]==NOCOLOR) && (BinaryIsOne(lpvar[e]))){
            pre[nxt] = e;
            q.push(nxt);
          }
        }
      }
      if (pre[targetnode] == INVALID) ok = false;
      else {
        DNode x = targetnode;
        while (x != sourcenode) {
          ecolor[pre[x]] = 3+i%6; // use colors 3(RED) to 8(CYAN), see myutils.h
          x = g.oppositeNode(x, pre[x]);
        }
      }
    }
    if (!ok) {cout << "Nao eh possivel encontrar os " << k << " caminhos!" << endl;}
    //------------------------------------------------------------------------------
    */
    for(ArcIt e(g); e!=INVALID; ++e) {
      if (BinaryIsOne(lpvar[e]) && ecolor[e]==NOCOLOR ) {
        cout << "Alguma(s) aresta(s) nao pertencem a qualquer caminho!" << endl;
        break;
      }
    }
    for(ArcIt e(g); e!=INVALID; ++e)
      if (NonBinary(lpvar[e])) ecolor[e] = GRAY;

    for(ArcIt e(g); e!=INVALID; ++e) {
      if (NonBinary(lpvar[e])) {
        cout << "Alguma(s) aresta(s) possuem valor fracionario!"
	  " (marcadas com cor cinza claro)" << endl;
        break;
      }
    }
    //-----------------------------------------------------------------


    cout << "kPaths Tree Value = " << soma << endl;
    ViewListDigraph(g,vname,px,py,vcolor,ecolor,
	"minimum kPaths cost in graph with "+IntToString(T.nnodes)+
	" nodes and "+IntToString(k)+" paths: "+DoubleToString(soma));
  } catch (...) {cout << "Error during callback..." << endl; }
  return 0;
}

