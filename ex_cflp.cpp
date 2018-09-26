// Project and Analysis of Algorithms
// Flávio Keidi Miyazawa
// Problems with connectivity: Capacitated Facility Location
#include <cstdio>
#include <string>
#include <queue>
#include "mygraphlib.h"
#include "myutils.h"
#include <lemon/lp.h>
#include <lemon/list_graph.h>
#include <lemon/concepts/digraph.h>
#include <gurobi_c++.h>
using namespace lemon;
using namespace std;

#if __cplusplus >= 201103L
#include <unordered_map>
#else
#include <tr1/unordered_map>
#endif


enum node_type_t { FACILITY, CLIENT };

void PulaBrancoComentario(ifstream &ifile);


// Read an instance for the Capacitated Facility Location Problem
bool ReadCFLPInstance(string filename, ListDigraph &g,
                      DNodeStringMap & vname,
                      ArcValueMap    & edge_weight,
                      DNodeValueMap  & facility_weight,
                      DNodeIntMap    & facility_capacity,
                      DNodePosMap  & posx,
                      DNodePosMap  & posy)
{
  int i,n,m;
  double peso;
  Arc a;
  string nomeu,nomev;

#if __cplusplus >= 201103L
  std::unordered_map<string,DNode> nodemap;
#else
  std::tr1::unordered_map<string,DNode> nodemap;
#endif


  
  ifstream ifile;
  ifile.open(filename.c_str());
  if (!ifile) {cout << "File '" << filename << "' does not exist.\n"; exit(0);}
  PulaBrancoComentario(ifile);
  ifile >> n;    ifile >> m; // first line have number of nodes and number of arcs
  if (m<0||ifile.eof()){cout<<"File "<<filename<<" is not a graph given by arcs.\n";
    exit(0);}
  string STR;
  DNode u,v;
  for (i=0;i<n;i++) {
    getline(ifile,STR);
    if (ifile.eof()) {cout<<"Reached unexpected end of file "<<filename<<".\n";exit(0);}
    while (STR=="") getline(ifile,STR);
    {
      string token;
      istringstream ins; // Declare an input string stream.
      ins.str(STR);        // Specify string to read.
      int nt = 0;
      while( getline(ins, token, ' ') ) {
	// format: <node_name>  <pos_x>  <pos_y>  <facility_weight>  <facility_capacity>
	if (nt==0) {
	  auto test = nodemap.find(token);
	  if (test!=nodemap.end()){cout<<"ERROR: Repeated node: "<<nomev<<endl;exit(0);}
	  v = g.addNode();
	  nodemap[token] = v;
	  vname[v] = token;}
	else if (nt==1) { posx[v] = atof(token.c_str());}
	else if (nt==2) { posy[v] = atof(token.c_str());}
	else if (nt==3) { facility_weight[v] = atof(token.c_str());}
	else if (nt==4) { facility_capacity[v] = atoi(token.c_str());}
	nt++;
      }
    }
  }
  for (i=0;i<m;i++) {
    // format: <cliente_node>   <facility_node>   <arc_weight>
    ifile >> nomeu;  ifile >> nomev; ifile >> peso;
    if (ifile.eof()) 
      {cout << "Reached unexpected end of file " <<filename << ".\n"; exit(0);}
    auto test = nodemap.find(nomeu);
    if (test == nodemap.end()) {cout<<"ERROR: Unknown node: "<<nomeu<<endl;exit(0);}
    else u = nodemap[nomeu];
    
    test = nodemap.find(nomev);
    if (test == nodemap.end()) {cout<<"ERROR: Unknown node: "<<nomev<<endl;exit(0);}
    else v = nodemap[nomev];
    
    a = g.addArc(u,v);     edge_weight[a] = peso;
  }
  ifile.close();
  return(true);
}




int main(int argc, char *argv[])
{
  Digraph g;  // graph declaration
  string digraph_kpaths_filename;
  if ((argc!=2) && (argc!=4)) {
    cout << endl << "Sintax to read instance from a file:" << endl ;
    cout << "      " << argv[0] << " <filename>" << endl << endl;
    
    cout << "Sintax to generate random instance:" << endl;
    cout << "      " << argv[0] << " <number_of_clients>  <number_of_facilities>  <capacity>" << endl << endl;
    cout << "Examples:" << endl;
    cout << "      " << argv[0] << " 100  30  15" << endl << endl;
    cout << "      " << argv[0] << " digr_cflp_1" << endl << endl;
  exit(0);
  }
  DNodeStringMap vname(g);  // name of graph nodes
  DNodePosMap px(g),py(g);  // xy-coodinates for each node
  DNodeColorMap vcolor(g);// color of nodes
  DNodeValueMap facility_weight(g);
  DNodeIntMap   facility_capacity(g);
  ArcColorMap   ecolor(g, NOCOLOR); // color of edges
  ArcValueMap   edge_weight(g);
  Digraph::NodeMap<node_type_t> vtype(g);

  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a minimization problem

  /* LPI variables */
  Digraph::ArcMap<GRBVar> x(g); // variable for connections, 1=connected, 0=not connected
  Digraph::NodeMap<GRBVar> y(g); // variable for facilities, 1=open, 0=closed

  if (argc==2) { ReadCFLPInstance(string(argv[1]),g,vname,edge_weight,
				  facility_weight, facility_capacity,px,py);

    for (DNodeIt v(g); v!=INVALID; ++v) {
      cout << vname[v] << "  " <<  px[v] << "  " <<  py[v] << "  " <<
	facility_weight[v]  << "  " <<  facility_capacity[v] << endl;
    }
    for (ArcIt a(g); a!=INVALID; ++a) {
      cout << vname[g.source(a)] << "  " <<  vname[g.target(a)] << "  " <<  edge_weight[a] << endl;
    }
    
  }else { // argc==3
    int nC,nF,Cap;
    nC = atoi(argv[1]);
    nF = atoi(argv[2]);
    Cap = atoi(argv[3]);
    vector <DNode> Client(1000);
    vector <DNode> Facility(1000);
    for (int i=0;i<nC;i++) {
      Client[i] = g.addNode();
      px[Client[i]] = drand48()*100;
      py[Client[i]] = drand48()*100;
      vname[Client[i]] = "c"+IntToString(i+1);
    }
    for (int j=0;j<nF;j++) {
      Facility[j] = g.addNode();
      px[Facility[j]] = drand48()*100;
      py[Facility[j]] = drand48()*100;
      facility_capacity[Facility[j]] = Cap;
      facility_weight[Facility[j]] = 100;
      vname[Facility[j]] = "f"+IntToString(j+1);
    }
    for (int j=0;j<nF;j++) {
      for (int i=0;i<nC;i++) {
	Arc a;
	a = g.addArc(Client[i],Facility[j]);
	edge_weight[a] = sqrt(
	      (px[Facility[j]]-px[Client[i]])*(px[Facility[j]]-px[Client[i]])+
	      (py[Facility[j]]-py[Client[i]])*(py[Facility[j]]-py[Client[i]]));
      }
    }
  }
  int nfac=0,ncli=0;
  for(DNodeIt v(g); v != INVALID; ++v) {
    int num_out = countOutArcs(g, v);
    int num_in = countInArcs(g, v);
    if (num_out > 0 && num_in > 0 ) {
      cerr << "ERRO: Existe vertice que eh instalacao e cliente!" << endl;
      exit(1);
    } else if (num_out > 0) {
      vtype[v] = CLIENT; ncli++;
    } else {
      vtype[v] = FACILITY; nfac++;
    }
  }
  
  for(DNodeIt v(g); v != INVALID; ++v) {
    if (vtype[v] == FACILITY) {
      y[v] = model.addVar(0.0, 1.0, facility_weight[v], GRB_BINARY);
    }
  }
  for(ArcIt e(g); e != INVALID; ++e) {
    x[e] = model.addVar(0.0, 1.0, edge_weight[e], GRB_BINARY);
  }
  model.update();
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //                           Altere daqui para baixo
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // OBS.: A entrada é dada por um digrafo, onde os nós são os clientes e instalações
  //       e os arcos são as conexões.
  // * Para identificar se um nó é um cliente ou instalação, há um vetor vtype para identificá-los.
  //   Se vtype[v] é FACILITY, então o vértice é uma instalação.
  //   Se vtype[v] é CLIENT, então o vértice é um cliente.
  // 
  // * As conexões possíveis são dadas por arcos que saem dos clientes e chegam nas possíveis
  //   instalações. Dado um nó v em um digrafo do LEMON, é possível percorrer apenas os
  //   arcos que entram em v.     E também é possível percorrer apenas os que saem de v.
  //
  // * Se i eh instalacao, então sua capacidade eh de no máximo facility_capacity[i] clientes
  //   

  for(DNodeIt v(g); v != INVALID; ++v) {
    if (vtype[v] == FACILITY) {
      GRBLinExpr num_customers = 0;
      for (InArcIt e(g, v); e != INVALID; ++e) {
        num_customers += x[e];
      }
      model.addConstr(num_customers <= facility_capacity[v]);
    } else if(vtype[v] == CLIENT) {
      GRBLinExpr num_facilities = 0;
      for(OutArcIt e(g, v); e != INVALID; ++e) {
          num_facilities += x[e];
          DNode fac = g.runningNode(e);
          model.addConstr(x[e] <= y[fac]);
        }
      model.addConstr(num_facilities >= 1);
    }
  }

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //                           Altere daqui para cima
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  try {
    model.optimize();

    double total_weight = 0.0;

    /* colorindo vertices */
    for(DNodeIt v(g); v != INVALID; ++v) {
      switch(vtype[v]) {
      case FACILITY:
        if (BinaryIsOne(y[v].get(GRB_DoubleAttr_X))) {
          vcolor[v] = RED;
          total_weight += facility_weight[v];
        } else vcolor[v] = MAGENTA;
        break;
      case CLIENT:
        vcolor[v] = BLUE;
        break;
      }
    }

    /* colorindo arestas */
    for(ArcIt e(g); e != INVALID; ++e) {
      if (BinaryIsOne(x[e].get(GRB_DoubleAttr_X))) {
        total_weight += edge_weight[e];
        ecolor[e] = BLACK;
      } else ecolor[e] = NOCOLOR;
    }

    bool capacity_ok = true;
    bool clients_satisfied = true;
    bool clients_connected_to_open = true;

    for(DNodeIt v(g); v != INVALID; ++v) {
      if (vtype[v] == FACILITY) {
        int num_customers = 0;
        for (InArcIt e(g, v); e != INVALID; ++e) 
	  if (BinaryIsOne(x[e].get(GRB_DoubleAttr_X))) num_customers++;
        if (num_customers > facility_capacity[v]) {
          vcolor[v] = GREEN;
          capacity_ok = false;
        }
      } else if(vtype[v] == CLIENT) {
        int num_facilities = 0;
        for(OutArcIt e(g, v); e != INVALID; ++e) {
	  if (BinaryIsOne(x[e].get(GRB_DoubleAttr_X))) {
	      num_facilities++;
            DNode fac = g.runningNode(e);
            double facval = y[fac].get(GRB_DoubleAttr_X);
            if (BinaryIsZero(facval)) 
              clients_connected_to_open = false;
          }
        }
        if (num_facilities <= 0) {
          vcolor[v] = RED;
          clients_satisfied = false;
        }
      }
    }
    cout << "Peso da solucao encontrada: " << total_weight << endl;
    if (!capacity_ok) {
      cout << "AVISO: Existe(m) instalacao(oes) com mais clientes que o permitido (marcadas de verde)." << endl;
    }
    if (!clients_satisfied) {
      cout << "AVISO: Existe(m) cliente(s) que nao foram conectados a uma instalacao (marcados de vermelho)." << endl;
    }
    if (!clients_connected_to_open) {
      cout << "AVISO: Existe(m) cliente(s) que foram conectados a alguma instalacao nao aberta." << endl;
    }
    
    ViewListDigraph(g,vname,px,py,vcolor,ecolor,"CFLP Bipartido Euclidiano. Clientes = AZUL, Instalacoes Abertas = VERMELHO  e  Nao_Abertas = ROSA. Peso_Solucao = "+DoubleToString(total_weight)); // esta rotina gera um eps que eh transformado para pdf
  } catch(GRBException e) {
    cerr << "Nao foi possivel resolver o PLI." << endl;
    cerr << "Codigo de erro = " << e.getErrorCode() << endl;
    cerr << e.getMessage();
  }

  return 0;
}

