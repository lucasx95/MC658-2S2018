#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;
int main(int argc, char *argv[]) 
{
  int n;
  srand48(clock());
  if (argc!=2) {cout<<"Usage: "<< argv[0]<<" <number_of_nodes>"<<endl; exit(0);} 
  n = atoi(argv[1]);
  cout << n << " -1" << endl;
  for (int i=0;i<n;i++) 
    // format of each line: Node_name  x-coordinate  y-coordinate
    cout << i+1 << " " <<  ((int) (drand48()*1000)+1) << " " << ((int) (drand48()*1000)+1) << endl;
}
