#include "Elm.hpp"

Elm::Elm() : ctr(3, 0.0), edges(3, -1) { }
Elm::~Elm() { }

void Elm::SetNodes(vector<unsigned long> &nodes) { swap(this->nodes, nodes); }

const vector<unsigned long>& Elm::GetNodes(void) const { return nodes; }

void Elm::SetCentroid(const double xc, const double yc, const double zc) {
    ctr[0] = xc; ctr[1] = yc; ctr[2] = zc;
}

const vector<double>& Elm::GetCentroid(void) const { return ctr;}

void Elm::SetArea(const double area) { this->area = area; }
const double Elm::Area(void) const { return area; }

void Elm::AddEdge(const unsigned long iedge, const unsigned long inode1, const unsigned long inode2) {
    //cout << "inputs: " << inode1 << "\t" << inode2 << endl;
    //cout << "elm_nodes: " << nodes[1] << "\t" << nodes[2] << "\t" << nodes[3] << endl;

    if ((nodes[0] == inode1 && nodes[1] == inode2) || (nodes[1] == inode1 && nodes[0] == inode2)) edges[0] = iedge;
    else if ((nodes[1] == inode1 && nodes[2] == inode2) || (nodes[2] == inode1 && nodes[1] == inode2)) edges[1] = iedge;
    else if ((nodes[2] == inode1 && nodes[0] == inode2) || (nodes[0] == inode1 && nodes[2] == inode2)) edges[2] = iedge;
    else {
        std::cout << "Edge not found committing self delete :(\n";
        exit(0);
    }

    //exit(0);
}
const vector<unsigned long>& Elm::Edges(void) const { return edges; }

unsigned short Elm::GetEdgeNum(const unsigned long n1, const unsigned long n2) const {
    if ((n1 == nodes[0] && n2 == nodes[1]) || (n1 == nodes[1] && n2 == nodes[0])) return 0;
    else if ((n1 == nodes[1] && n2 == nodes[2]) || (n1 == nodes[2] && n2 == nodes[1])) return 1;
    else if ((n1 == nodes[2] && n2 == nodes[0]) || (n1 == nodes[0] && n2 == nodes[2])) return 2;
    else cout << "Error in GetEdgeNum!\n"; exit(0); // panic
}
