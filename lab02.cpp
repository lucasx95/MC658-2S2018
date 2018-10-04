#include <gurobi_c++.h>
#include <fstream>
#include <vector>
#include <set>

using namespace std;

vector<set<int>> solveLab02(int m, int W, int n, const vector<int> &w, const vector<int> &b, const vector<int> &v);

int main(int argv, char *argc[]){

	ifstream file(argc[1], ios::out);

	int m, W, n;
	vector<int> w, b, v;
	
	file >> n;
	file >> W >> m;
	
	w = vector<int>(n);
	b = vector<int>(n);
	v = vector<int>(n);
	
	for (int i=0; i<n; i++)
		file >> w[i] >> b[i] >> v[i];
		
	// Output input data
	cout<<"\n++ ENTRADA: "<<endl;
	cout<<"Numero de propagandas: "<<n<<endl;
	cout<<"Numero de slots: "<<m<<endl;
	cout<<"Capacidade de cada slot: "<<W<<endl;
	
	for (int i=0; i<n; i++){
		cout<<"Propaganda "<<i<<":"<<endl;
		cout<<"  w = "<<w[i]<<" | b = "<<b[i]<< " | v = "<<v[i]<<endl;
	}
	
	auto sol = solveLab02(m, W, n, w, b, v);
	
	
	cout<<"\nSolução encontrada: "<<endl;
	
	if (sol.size()!=m){
		cout<<"solução com número de slots inválido!"<<endl;
		cout<<"-1"<<endl;
		exit(0);	
	}
	
	// Output da solução
	for (int s=0; s<m; s++){
		cout<<"Slot "<<s<<":"<<endl;
		for (auto i: sol[s]) cout<<i<<" ";
		cout<<endl;
	}
	
	// Checar se todos os indices estão dentro dos limites
	for (auto slot: sol)
		for (auto i: slot)
			if (i < 0 || i >= n){
				cout<<"slot contém indice de propaganda inexistente!"<<endl;
				cout<<"-1"<<endl;
				exit(0);	
			}
		
	// Checar se o número de propagandas de cada tipo está dentro do limite
	vector<int> countPropagandas(n);
	for (auto slot: sol)
		for (auto i: slot)
				countPropagandas[i]++;
				
	for (int i=0; i<n; i++) if (countPropagandas[i] > b[i]){
		cout<<"solução com número de cópias inválido de propagandas "<<i<<endl;
		cout<<"-1"<<endl;
		exit(0);	
	}
	
	// Checar se a capacidade de cada slot é respeitada
	for (auto slot: sol){
		int soma= 0;
		
		for (auto i: slot)
			soma+=w[i];
			
		if (soma >= W){
			cout<<"Slot com capacidade excedida"<<endl;
			cout<<"-1"<<endl;
			exit(0);	
		}	
	}
	
	int valorSol = 0;
	
	for (int i=0; i<n; i++)
		valorSol += countPropagandas[i]*v[i];
	
	cout<<"Valor da solução:"<<endl;
	cout<<valorSol<<endl;
	
}

// Essa função deve ser modificada para retornar uma solução ótima para o problema do Lab02.
// A solução consiste de um vetor de conjuntos.
// Cada posição do vetor é equivalente a um slot.
// Cada conjunto contém as propagandas que vão aparecer no respectivo slot.
/**
 * max Σ(i, 0 -> n)Σ(j, 0 -> m) Xij*v[i] //Maximize the whole value
 *
 * Σ(i, 0 -> n) Xij*w[i] < W, ∀(j, 0 -> m) // The weight of the adds in a slot cannot surpass W
 * Σ(j, 0 -> m) Xij < b[i], ∀(i, 0 -> n)   // The number of the adds cannot surpass b[i]
 *
 * Xij ∈ {0,1} // If a add i is in the slot j
 *
 * @param m number of slots
 * @param W max weight per slot
 * @param n number of ads
 * @param w weight of each add
 * @param b max number of each add
 * @param v value of each add
 * @return vector of slots with which add are in each
 */
vector<set<int>> solveLab02(int m, int W, int n, const vector<int> &w, const vector<int> &b, const vector<int> &v){

	vector<set<int>> sol(m);
	
	
	// Exemplo de implementação do Gurobi
	// O modelo implementado não está relacionado ao problema, é apenas um exemplo de implementação de um modelo usando o gurobi, e deve ser modificado para resolver o problema das propagandas.
	
	GRBEnv *env = new GRBEnv();
	GRBModel *model = new GRBModel(*env);
	
	GRBVar X = model->addVar(0.0, GRB_INFINITY, 25.0, GRB_CONTINUOUS);
	GRBVar Y = model->addVar(0.0, GRB_INFINITY, 12.0, GRB_CONTINUOUS);
	model->addConstr(X + Y >= 30.0);
	model->addConstr(X - Y <= 23.5);
	model->addConstr(Y - X <= 23.5);
	
	cout<<"\n+ Solving the model:\n"<<endl;
	model->optimize();
	
	cout<<"\n+ Model Solution : "<<endl;
	cout<<"X = "<<X.get(GRB_DoubleAttr_X)<<endl;
	cout<<"Y = "<<Y.get(GRB_DoubleAttr_X)<<endl;
	
	
	
	delete model;
	delete env;

	//+ Implementação de uma solução trivial para o problema:
	int slot = 0;
	int soma = 0;
	
	for (int i=0; i<n; i++) for (int b_=0; b_<b[i]; b_++) if (slot < m) if (sol[slot].count(i) == 0){
		if (soma+w[i] <= W){
			sol[slot].insert(i);
			soma+=w[i];
		}
		else{
			soma = 0;
			slot++;
		}
	}
	
	return sol;

}