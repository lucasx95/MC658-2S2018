//----------------------------------------------------------------------
// Example of an exact program to solve relaxation of a packing problem
// with GUROBI Solver. See pages 105-108, 124-125 of 
// http://www.ic.unicamp.br/~fkm/lectures/proglin.pdf  (in portuguese)
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------
#include <iostream>
using namespace std;
#include <gurobi_c++.h>
int main() 
{
  int n=4;
  vector<GRBVar> x(n);
  GRBEnv env = GRBEnv();
  GRBModel model = GRBModel(env);
  cout << endl << endl <<"Relaxation of the integer linear program in page 124 of the slides in " << endl <<
    "http://www.ic.unicamp.br/~fkm/lectures/proglin.pdf" << endl << endl << endl;

  for (int i=0;i<n;i++) x[i]=model.addVar(0.0,GRB_INFINITY,1.0,GRB_CONTINUOUS,"");
  model.update(); // run update to use model inserted variables
  model.addConstr(3*x[0] + 8*x[1] +          2*x[3] >= 950);
  model.addConstr(2*x[0] +            x[2] + 3*x[3] >= 1150);
  model.addConstr(  x[0] +          3*x[2] +   x[3] >= 495);
  model.addConstr(  x[0] +          2*x[2] +   x[3] >= 450);

  model.set(GRB_StringAttr_ModelName, "Packing Example"); 
  model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); 
  model.optimize();

  for (int i=0;i<n;i++)
    cout << "x[" << i+1 << "] = "<< x[i].get(GRB_DoubleAttr_X) << endl;
}
