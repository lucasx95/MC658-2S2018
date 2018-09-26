//----------------------------------------------------------------------
// Some basic routines using non-oriented graphs
// Send comments/corrections to Flavio K. Miyazawa.
//----------------------------------------------------------------------
#include <gurobi_c++.h>
#include <float.h>
#include <math.h>
#include <set>
#include <unordered_map>
#include <lemon/list_graph.h>
#include <lemon/gomory_hu.h>
#include "mygraphlib.h"

struct ListNode {
    Node data;
    int value = 0;
    int weight = 0;
    float valuePerWeight = 0.0;
    ListNode *next = nullptr;
    ListNode *previous = nullptr;

};

/*
 * Simple Linked List with insertOrdered option
 */
class OrderedLinkedNodeList {
private:
    ListNode *head = nullptr;
    ListNode *tail = nullptr;
    int size = 0;

public:
    OrderedLinkedNodeList() {};

    void insertOrdered(ListNode *node) {
        size++;
        // if the list is empty initialize it with both head and tail pointing to node
        if (head == nullptr) {
            head = node;
            tail = node;
            return;
        }
        ListNode *current = head;
        while (current != nullptr && current->valuePerWeight > node->valuePerWeight) {
            current = current->next;
        }

        // if its inserting in the start update head
        if (current == head) {
            node->next = head;
            head->previous = node;
            head = node;
            return;
        }

        // if its inserting in the end update tail
        if (current == nullptr) {
            tail->next = node;
            node->previous = tail;
            tail = node;
            return;
        }


        node->next = current;
        node->previous = current->previous;
        current->previous->next = node;
        current->previous = node;
    }

    // Inserts at the end. Only used when it is sure 'node' is the lowest valuable per weight
    void insert(ListNode *node) {
        size++;
        if (head == nullptr) {
            head = node;
            tail = node;
            return;
        }
        tail->next = node;
        node->previous = tail;
        tail = node;
    }

    ListNode *peak() {
        return head;
    }

    ListNode *peakBottom() {
        return tail;
    }

    ListNode *top() {
        // if the list was empty return nullptr
        if (empty()) {
            return nullptr;
        }
        size--;
        ListNode *fetch_node = head;
        // if its the last element clear pointers and returns
        if (empty()) {
            head = nullptr;
            tail = nullptr;
            return fetch_node;
        }
        head = head->next;
        head->previous = nullptr;
        fetch_node->next = nullptr;
        fetch_node->previous = nullptr;
        return fetch_node;
    }

    ListNode *bottom() {
        // if the list was empty return nullptr
        if (empty()) {
            return nullptr;
        }
        size--;
        ListNode *fetch_node = tail;
        // if its the last element clear pointers and returns
        if (empty()) {
            head = nullptr;
            tail = nullptr;
            return fetch_node;
        }
        tail = tail->previous;
        tail->next = nullptr;
        return fetch_node;
    }

    bool canInsertInSolution(ListNode *node, const ListGraph &graph, const vector<vector<int>> &edges,
                             int remaining_weight) {
        if (remaining_weight < node->weight) {
            return false;
        }
        ListNode *current = head;
        while (current != nullptr) {
            if (edges[lemon::ListGraph::id(node->data)][lemon::ListGraph::id(current->data)] == 1) {
                return false;
            }
            current = current->next;
        }
        return true;
    }

    int length() {
        return size;
    }

    bool empty() {
        return size == 0;
    }

    void clear() {
        head = tail = nullptr;
        size = 0;
    }

    // Estimate max value without looking at the edges
    int estimate(int remaining_weight) {
        int estimative = 0;
        ListNode *current = head;
        while (remaining_weight > 0 && current != nullptr) {
            // if the item does not fit completely estimate it partially
//            if (current->weight > remaining_weight) {
//                float partial_value = (float) remaining_weight / (float) current->weight;
//                estimative += (int) ceil((float) current->value * partial_value);
//                remaining_weight = 0;
//                current = current->next;
//                continue;
//            }
            estimative += current->value;
            remaining_weight -= current->weight;
            current = current->next;
        }
        return estimative;
    }

    set<Node> toSet() {
        set<Node> independentSet;
        ListNode *current = head;
        while (current != nullptr) {
            independentSet.insert(current->data);
            current = current->next;
        }
        return independentSet;
    }

    ListNode *remove(ListNode *node) {
        size--;
        if (node == head) {
            head = node->next;
        } else {
            node->previous->next = node->next;
        }
        if (node == tail) {
            tail = node->previous;
        } else {
            node->next->previous = node->previous;
        }
        node->next = nullptr;
        node->previous = nullptr;
        return node;
    }

    void copy(OrderedLinkedNodeList list) {
        head = list.head;
        tail = list.tail;
        size = list.size;
    }

    void print() {
        ListNode *current = head;
        while (current != nullptr) {
            cout << lemon::ListGraph::id(current->data) << "  ";
            current = current->next;
        }
        cout << endl;
    }

    void print(NodeStringMap &vname) {
        ListNode *current = head;
        while (current != nullptr) {
            cout << vname[current->data] << "  ";
            current = current->next;
        }
        cout << endl;
    }

};

bool ReadListGraph3(string filename,
                    ListGraph &g,
                    NodeStringMap &vname,
                    NodeIntMap &weight,
                    NodeIntMap &value,
                    NodePosMap &posx,
                    NodePosMap &posy,
                    int &C) {
    ifstream ifile;
    int n, m;
    Edge a;
    string STR;
    Node u, v;
#if __cplusplus >= 201103L
    std::unordered_map<string, Node> string2node;
#else
    std::tr1::unordered_map<string,Node> string2node;
#endif
    ifile.open(filename.c_str(), ios::in);
    if (!ifile) {
        cout << "File '" << filename << "' does not exist.\n";
        exit(0);
    }
    ifile >> n;
    ifile >> m;
    ifile >> C; // first line have number of nodes, number of arcs and capacity

    if (m < 0 || ifile.eof()) {
        cout << "File " << filename << " is not a digraph given by arcs.\n";
        exit(0);
    }

    // continue to read ifile, and insert information in g for the nex t n nodes
    for (int i = 0; i < n; i++) {

        v = g.addNode();
        string token;
        ifile >> token;
        string2node[token] = v;
        vname[v] = token;
        ifile >> weight[v];
        ifile >> value[v];
        posx[v] = 0.0;
        posy[v] = 0.0;

    }

    // continue to read ifile, and obtain edge weights
    string nomeu, nomev;
    for (int i = 0; i < m; i++) {
        // format: <node_u>   <node_v>   <edge_weight>
        ifile >> nomeu;
        ifile >> nomev;
        if (ifile.eof()) {
            cout << "Reached unexpected end of file.\n";
            exit(0);
        }
        auto test = string2node.find(nomeu);
        if (test == string2node.end()) {
            cout << "ERROR: Unknown node: " << nomeu << endl;
            exit(0);
        } else u = string2node[nomeu];

        test = string2node.find(nomev);
        if (test == string2node.end()) {
            cout << "ERROR: Unknown node: " << nomev << endl;
            exit(0);
        } else v = string2node[nomev];
        a = g.addEdge(u, v);
    }

    // if there exists some node without pre-defined position,
    // generate all node positions
    for (NodeIt v(g); v != INVALID; ++v) {
        if (posx[v] == DBL_MAX || posy[v] == DBL_MAX) {
            ListGraph t;
            EdgeValueMap temp(t);
            GenerateVertexPositions(g, temp, posx, posy);
            break;
        }
    }
    ifile.close();
    return (true);
}

set<Node>
max_ind_set(const ListGraph &g, const NodeIntMap &weight, const NodeIntMap &value, int Capacity, NodeStringMap &vname);


bool isSetIndependent(ListGraph &g, const set<Node> &indSet) {

    for (EdgeIt e(g); e != INVALID; ++e)
        if (indSet.count(g.u(e)) && indSet.count(g.v(e)))
            return false;

    return true;

}


int main(int argc, char *argv[]) {
    ListGraph g;
    NodeIntMap weight(g);
    NodeIntMap value(g);
    NodeStringMap vname(g);
    EdgeStringMap ename(g);
    NodePosMap posx(g), posy(g);
    string filename;
    int C;

    // uncomment one of these lines to change default pdf reader, or insert new one
    //set_pdfreader("open");    // pdf reader for Mac OS X
    //set_pdfreader("xpdf");    // pdf reader for Linux
    set_pdfreader("evince");  // pdf reader for Linux

    if (argc != 2) {
        cout << endl << "Usage: " << argv[0] << " <graph_filename>" << endl << endl <<
             "Example: " << argv[0] << " gr_7" << endl <<
             "         " << argv[0] << " gr_70" << endl << endl;
        exit(0);
    } else if (!FileExists(argv[1])) {
        cout << "File " << argv[1] << " does not exist." << endl;
        exit(0);
    }
    filename = argv[1];

    // Read the graph
    if (!ReadListGraph3(filename, g, vname, weight, value, posx, posy, C)) {
        cout << "Error reading graph file " << argv[1] << "." << endl;
        exit(0);
    }

    cout << "List of Nodes\n";
    for (NodeIt v(g); v != INVALID; ++v) {
        //vname[v] = vname[v] + "_" + IntToString(weight[v]) + "_" + IntToString(value[v]);
        cout << "Node " << vname[v] << " - Value " << value[v] << " - Weight " << weight[v] << endl;
    }
    cout << "\n==============================================================\n\n";

    cout << "List of Edges\n";
    for (EdgeIt e(g); e != INVALID; ++e) {
        ename[e] = "{" + vname[g.u(e)] + "," + vname[g.v(e)] + "}";
        cout << ename[e] << "  ";
    }
    cout << "\n==============================================================\n\n";


    cout << "List of Edges incident to nodes\n";
    for (NodeIt v(g); v != INVALID; ++v) {
        cout << "Node " << vname[v] << ": ";
        for (IncEdgeIt e(g, v); e != INVALID; ++e) cout << ename[e] << "  ";
        cout << "\n\n";
    }
    cout << "==============================================================\n\n";

    cout << "Graph file: " << filename << "\n\n";

    clock_t start = clock();
    auto independentSet = max_ind_set(g, weight, value, C, vname);
    cout << "time: " << (clock() - start) / (double) CLOCKS_PER_SEC << '\n';

    cout << "Independent set has vertices:" << endl;
    for (auto v: independentSet)
        cout << vname[v] << " ";

    cout << endl;

    if (!isSetIndependent(g, independentSet)) {
        cout << "Set is not independent" << endl;
        cout << -1 << endl;
        return 0;
    }

    int sum = 0;
    int weig = 0;
    for (auto v: independentSet) {
        sum += value[v];
        weig += weight[v];
    }

    if (weig > C) {
        cout << "Total weight exceeds the capacity" << endl;
        cout << -1 << endl;
        return 0;
    }
    cout << "\nSolution weight:" << endl << weig << endl;
    cout << "Solution value:" << endl << sum << endl;

    return 0;

}


// Modificar essa subrotina que retorna o conjunto independente de valor máximo
// O código a seguir é apenas um exemplo de uma solução trivial
set<Node>
max_ind_set(const ListGraph &g, const NodeIntMap &weight, const NodeIntMap &value, int Capacity, NodeStringMap &vname) {
    int max_solution = 0, remaining_weight = Capacity, current_solution = 0;
    OrderedLinkedNodeList available, solution, used;
    set<Node> independentSet;

    int numberOfNodes = 0;
    for (NodeIt v(g); v != INVALID; ++v) {
        ListNode *listNode = new ListNode;
        listNode->data = v;
        listNode->weight = weight[v];
        listNode->value = value[v];
        listNode->valuePerWeight = (float) listNode->value / (float) listNode->weight;
        available.insertOrdered(listNode);
        numberOfNodes++;
    }

    vector<vector<int>> edges(static_cast<unsigned long>(numberOfNodes), vector<int>(
            static_cast<unsigned long>(numberOfNodes)));
    for (EdgeIt e(g); e != INVALID; ++e) {
        Node v = g.v(e);
        Node u = g.u(e);
        edges[g.id(v)][g.id(u)] = 1;
        edges[g.id(u)][g.id(v)] = 1;
    }

    ListNode *clean_backtrack = nullptr;
    while (!available.empty()) {
        ListNode *candidate = available.peak();
        cout << endl << "PRE ITERECTION" << endl;
        cout << "AVAILABLE: ";
        available.print(vname);
        cout << "USED: ";
        used.print(vname);
        cout << "SOLUTION: ";
        solution.print(vname);
        cout << "CLEAN BACK: ";
        cout << (clean_backtrack == nullptr ? "-1" : vname[clean_backtrack->data])<< endl;


        while (candidate != nullptr &&
               current_solution + available.estimate(remaining_weight) > max_solution) {
            ListNode *next = candidate->next;
            if (solution.canInsertInSolution(candidate, g, edges, remaining_weight)) {
                solution.insert(available.remove(candidate));
                remaining_weight -= candidate->weight;
                current_solution += candidate->value;
            }
            candidate = next;
        }

        cout << "POST ITERECTION" << endl;
        cout << "AVAILABLE: ";
        available.print(vname);
        cout << "USED: ";
        used.print(vname);
        cout << "SOLUTION: ";
        solution.print(vname);
        cout << "CLEAN BACK: ";
        cout << (clean_backtrack == nullptr ? "-1" : vname[clean_backtrack->data]) << endl;


        if (current_solution >= max_solution) {
            max_solution = current_solution;
            independentSet = solution.toSet();
            cout << "NEW BEST SOLUTION: ";
            solution.print(vname);
            cout << "VALUE: " << current_solution << endl;
        }

        if (clean_backtrack == solution.peakBottom()) {
            if (solution.empty()) {
                return independentSet;
            }
            clean_backtrack = nullptr;
            while (!available.empty()) {
                if (used.empty()) {
                    clean_backtrack = solution.peakBottom();
                }
                used.insertOrdered(available.top());
            }
            available.copy(used);
            used.clear();
        }

        ListNode *backtracked = solution.bottom();
        current_solution -= backtracked->value;
        remaining_weight += backtracked->weight;
        if (!solution.empty()) {
            if (used.empty()) {
                clean_backtrack = solution.peakBottom();
            }
            used.insertOrdered(backtracked);
        }
    }

    return independentSet;
}
