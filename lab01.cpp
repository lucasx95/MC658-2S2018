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

/**
 * A List Node with denormalized data for easy access
 */
struct ListNode {
    Node data;
    int value = 0;
    int weight = 0;
    // Value/weight ratio
    float valuePerWeight = 0.0;
    ListNode *next = nullptr;
    ListNode *previous = nullptr;

};

/*
 * Ordered Linked List with helper methods for IndependentSet.
 */
class OrderedLinkedNodeList {
private:
    ListNode *head = nullptr;
    ListNode *tail = nullptr;
    int size = 0;

public:
    OrderedLinkedNodeList() {};

    // Insert a node ordered by valuePerWeight.
    void insertOrdered(ListNode *node) {
        size++;
        // If the list is empty initialize it with both head and tail pointing to node.
        if (head == nullptr) {
            head = node;
            tail = node;
            return;
        }
        ListNode *current = head;
        while (current != nullptr &&
               (current->valuePerWeight > node->valuePerWeight ||
                (current->valuePerWeight == node->valuePerWeight && current->weight < node->weight) ||
                (current->valuePerWeight == node->valuePerWeight && current->weight == node->weight &&
                 lemon::ListGraph::id(current->data) < lemon::ListGraph::id(node->data)))) {
            current = current->next;
        }

        // If its inserting in the start update head.
        if (current == head) {
            node->next = head;
            head->previous = node;
            head = node;
            return;
        }

        // If its inserting in the end update tail.
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

    // Inserts at the end. Only used when it is sure 'node' has the lowest valuePerWeight.
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

    // Get head without removing.
    ListNode *peak() {
        return head;
    }

    // Get tail without removing.
    ListNode *peakBottom() {
        return tail;
    }

    // Get head and remove.
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

    // Get tail and remove.
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

    // Verify if can insert the node in this Solution list. Checks if it will continue independent and with weight below capacity.
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

    // List length.
    int length() {
        return size;
    }

    // REturn true if list is empty.
    bool empty() {
        return size == 0;
    }

    // Clear list.
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
            if (current->weight > remaining_weight) {
                float partial_value = (float) remaining_weight / (float) current->weight;
                estimative += (int) ceil((float) current->value * partial_value);
                remaining_weight = 0;
                current = current->next;
                continue;
            }
            estimative += current->value;
            remaining_weight -= current->weight;
            current = current->next;
        }
        return estimative;
    }

    // Convert this list to set<Node>
    set<Node> toSet() {
        set<Node> independentSet;
        ListNode *current = head;
        while (current != nullptr) {
            independentSet.insert(current->data);
            current = current->next;
        }
        return independentSet;
    }

    // Remove node from this list.
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

    // Soft copy of this list.
    void copy(OrderedLinkedNodeList list) {
        head = list.head;
        tail = list.tail;
        size = list.size;
    }

    // Print the current nodes in this list.
    void print() {
        ListNode *current = head;
        while (current != nullptr) {
            cout << lemon::ListGraph::id(current->data) << "  ";
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
max_ind_set(const ListGraph &g, const NodeIntMap &weight, const NodeIntMap &value, int Capacity);


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

    auto independentSet = max_ind_set(g, weight, value, C);

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
max_ind_set(const ListGraph &g, const NodeIntMap &weight, const NodeIntMap &value, int Capacity) {
    // Initialize variable in empty solution state.
    int max_solution = 0, remaining_weight = Capacity, current_solution = 0;
    OrderedLinkedNodeList available, solution, used;
    set<Node> independentSet;

    int numberOfNodes = 0;
    int min_weigth = INT32_MAX;
    // Add all nodes to the availability list that is ordered by valuePerWeight.
    for (NodeIt v(g); v != INVALID; ++v) {
        ListNode *listNode = new ListNode;
        listNode->data = v;
        listNode->weight = weight[v];
        listNode->value = value[v];
        listNode->valuePerWeight = (float) listNode->value / (float) listNode->weight;
        if (listNode->weight < min_weigth) {
            min_weigth = listNode->weight;
        }
        available.insertOrdered(listNode);
        numberOfNodes++;
    }

    // Add all edges to a matrix to easy quick access.
    vector<vector<int>> edges(static_cast<unsigned long>(numberOfNodes), vector<int>(
            static_cast<unsigned long>(numberOfNodes)));
    for (EdgeIt e(g); e != INVALID; ++e) {
        Node v = g.v(e);
        Node u = g.u(e);
        edges[lemon::ListGraph::id(v)][lemon::ListGraph::id(u)] = 1;
        edges[lemon::ListGraph::id(u)][lemon::ListGraph::id(v)] = 1;
    }

    // Initialize the backtrack.

    // clean_backtrack stores the 'primary' node in the solution, the one that was the last in the solution when
    // the used list was cleared the last time, once the solution bottom is the same as clean_backtrack once again,
    // all the possibilities with this partial solution including this node has been fulfilled and the 'used' list
    // is restored as available.
    ListNode *clean_backtrack = nullptr;
    // Continue in the loop while there's available nodes to be inserted in the solution.
    while (!available.empty()) {
        // Get a candidate for inserting. And iterates to the next ones until the end or it is discarted in the
        // estimative or the remaining_weigth available is smaller the the ligther node.
        ListNode *candidate = available.peak();
        while (candidate != nullptr &&
               current_solution + available.estimate(remaining_weight) > max_solution &&
               remaining_weight >= min_weigth) {
            ListNode *next = candidate->next;
            if (solution.canInsertInSolution(candidate, g, edges, remaining_weight)) {
                solution.insert(available.remove(candidate));
                remaining_weight -= candidate->weight;
                current_solution += candidate->value;
            }
            candidate = next;
        }

        // If the solution found is better then the best known update the best known.
        if (current_solution >= max_solution) {
            max_solution = current_solution;
            independentSet = solution.toSet();
        }

        // If matches the condition described above backtracks restoring available list from used.
        bool backtrack_cleared = false;
        if (clean_backtrack == solution.peakBottom()) {
            // If at this time the solution is empty there's no other viable solution candidate so it returns.
            if (solution.empty()) {
                return independentSet;
            }
            // Restore the clean_backtrack to the empty state and merge used and available list as available.
            clean_backtrack = nullptr;
            while (!available.empty()) {
                used.insertOrdered(available.top());
            }
            available.copy(used);
            used.clear();
            backtrack_cleared = true;
        }

        ListNode *backtracked = solution.bottom();
        current_solution -= backtracked->value;
        remaining_weight += backtracked->weight;
        if (!solution.empty()) {
            // If used was empty set clean_backtrack.
            if (used.empty()) {
                clean_backtrack = solution.peakBottom();
            }
            used.insertOrdered(backtracked);
            // If it was not a complete backtrack, restore the nodes that were added after (valuePerWeight >) to the available list.
            if (!backtrack_cleared) {
                ListNode *restore = backtracked->next, *nextRestore = nullptr;
                while (restore != nullptr) {
                    nextRestore = restore->next;
                    available.insertOrdered(used.remove(restore));
                    restore = nextRestore;
                }
            }
        }
    }

    return independentSet;
}
