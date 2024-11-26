#include "MshReader.hpp"

MshReader::MshReader() { }
MshReader::MshReader(string filename) { GridReader::filename = filename; }
MshReader::~MshReader() { }

void MshReader::ReadFile() {
    nDim = 2;
    ifstream file(filename);
    if (!file.is_open()) { 
        cout << "Failed to open file o7!\n";
        5 / 0;
    }

    nNode = stoi(FindTag("$Nodes", file));
    nodes.resize(nNode);
    cout << "nNode: " << nNode << endl;

    string line;
    double xyz[5];
    for (auto iLine = 0ul; iLine < nNode; iLine++) {
        getline(file, line);
        istringstream iss(line);
        double value;
        size_t counter = 0;

        while (iss >> value) xyz[counter++] = value;
        Node &node = nodes[iLine];
        node.SetXYZ(xyz[1], xyz[2]); // first value is index so skip xyz[0]
    }

    nBoundaries = 2; // Only works for airfoil. I only care about msh files for testing
    bounds.resize(nBoundaries);

    unsigned long nTotElm = stoi(FindTag("$Elements", file)); // includes boundary edges
    unsigned long nBoundEdges = 0;
    unsigned long conn[MAX_LINE_LEN];

    unsigned long nElmInt = 0;
    vector<vector<unsigned long>> elmNodes;

    vector<vector<Edge>> b_edges;
    b_edges.resize(nBoundaries);
    for (auto iLine = 0ul; iLine < nTotElm; iLine++) {
        getline(file, line);
        istringstream iss(line);
        unsigned long value;
        size_t counter = 0;

        while (iss >> value) conn[counter++] = value;
        if (conn[3] == 20001) {
            vector<unsigned long> tri_node(3, 0);
            tri_node[0] = conn[5] - 1;
            tri_node[1] = conn[6] - 1;
            tri_node[2] = conn[7] - 1;
            elmNodes.push_back(tri_node);
        }
        else {
            Edge edge;
            edge.SetNodes(conn[5] - 1, conn[6] - 1);
            if (conn[3] == 10000) b_edges[0].push_back(edge); 
            else if (conn[3] == 30000) b_edges[1].push_back(edge);
            else 5 / 0; // panic
        }

        bounds[0].SetEdges(b_edges[0]);
        bounds[1].SetEdges(b_edges[1]);
    }

    nElm = elmNodes.size();
    elms.resize(nElm);
    vtk.resize(nElm);
    for (auto iElm = 0ul; iElm < nElm; iElm++) {
        elms[iElm].SetNodes(elmNodes[iElm]);
        vtk[iElm] = 5;
    }
    cout << "nElm: " << nElm << endl;


    ofstream f("msh_check.txt");
    f << "Nodes\n";
    for (auto iNode = 0ul; iNode < nNode; iNode++) {
        Node &node = nodes[iNode];
        f << node.X() << "\t" << node.Y() << "\t" << node.Z() << endl;
    }

    f << "Elements\n";
    for (auto iElm = 0ul; iElm < nElm; iElm++) {
        Elm &elm = elms[iElm];
        const auto eN = elm.GetNodes();
        f << eN[0] << "\t" << eN[1] << "\t" << eN[2] << endl;
    }

    f << "Boundaries\n";
    for (auto iBound = 0ul; iBound < nBoundaries; iBound++) {
        f << "Boundary: " << iBound << endl;
        Bound &bound = bounds[iBound];
        const auto bedg = bound.GetEdges();
        for (auto ibedg = 0ul; ibedg < bedg.size(); ibedg++) {
            f << bedg[ibedg].Node1() << "\t" << bedg[ibedg].Node2() << endl;
        }
    }



    f.close();

}

string MshReader::FindTag(const string& tag, ifstream &file) const {
    string line;
    while (true) {
        getline(file, line);
        istringstream iss(line);
        if (line.empty()) continue;
        if (line.substr(0, tag.length()) == tag) {
            getline(file, line);
            return line;
            break;
        }
    }
}
