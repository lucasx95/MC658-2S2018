//----------------------------------------------------------------------
// Some basic routines using non-oriented graphs
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------
#include <gurobi_c++.h>
#include <float.h>
#include <math.h>
#include <set>
#include <lemon/list_graph.h>
#include <lemon/gomory_hu.h>
#include "mygraphlib.h"
#include "myutils.h"


int main(int argc, char *argv[]) 
{
  ListGraph g;
  EdgeValueMap weight(g);
  NodeStringMap vname(g);
  EdgeStringMap ename(g);
  NodePosMap   posx(g),posy(g);
  string filename;

  // uncomment one of these lines to change default pdf reader, or insert new one
  //set_pdfreader("open");    // pdf reader for Mac OS X
  //set_pdfreader("xpdf");    // pdf reader for Linux
  set_pdfreader("evince");  // pdf reader for Linux

  if (argc!=2) {cout<< endl << "Usage: "<< argv[0]<<" <graph_filename>"<<endl << endl <<
      "Example: " << argv[0] << " gr_7" << endl <<
      "         " << argv[0] << " gr_70" << endl << endl; exit(0);}
  
  else if (!FileExists(argv[1])) {cout<<"File "<<argv[1]<<" does not exist."<<endl; exit(0);}
  filename = argv[1];
  
  // Read the graph
  if (!ReadListGraph(filename,g,vname,weight,posx,posy)) 
    {cout<<"Error reading graph file "<<argv[1]<<"."<<endl;exit(0);}

  cout << "List of Nodes\n";
  for (NodeIt v(g); v!=INVALID; ++v)
    cout << vname[v] << "  ";
  cout << "\n==============================================================\n\n";

  cout << "List of Edges\n";
  for (EdgeIt e(g); e!=INVALID; ++e) {
    ename[e] = "{"+vname[g.u(e)]+","+vname[g.v(e)]+"}_"+DoubleToString(weight[e]);
    cout << ename[e] << "  ";
  }
  cout << "\n==============================================================\n\n";


  cout << "List of Edges incident to nodes\n";
  for (NodeIt v(g); v!=INVALID; ++v) {
    cout << "Node " << vname[v] << ": ";
    for (IncEdgeIt e(g,v); e!=INVALID; ++e) cout << ename[e] << "  ";
    cout << "\n\n";
  }
  cout << "==============================================================\n\n";

  cout <<  filename << "\n";
  {
    NodeColorMap vcolor(g);
    EdgeColorMap ecolor(g);
    for (NodeIt v(g); v!=INVALID; ++v) vcolor[v] = RED;
    for (EdgeIt e(g); e!=INVALID; ++e) ecolor[e] = BLUE;
    ViewListGraph(g,vname,ename,posx,posy,vcolor,ecolor,"Grafo '"+filename+"' com "+IntToString(countNodes(g))+" vertices e "+IntToString(countEdges(g))+" arestas.");
  }

}
