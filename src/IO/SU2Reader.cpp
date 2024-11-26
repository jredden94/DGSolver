#include "SU2Reader.hpp"
#include <stdexcept>
#include <string>

SU2Reader::SU2Reader() { }
SU2Reader::SU2Reader(string filename) { GridReader::filename = filename; }
SU2Reader::~SU2Reader() { }

string SU2Reader::FindTag(const string& tag, ifstream &file) const {
    string line;
    while (true) {
        getline(file, line);
        istringstream iss(line);
        if (line.empty() || line[0] == '%') continue;
        if (line.substr(0, tag.length()) == tag) {
            size_t pos = line.find('=');
            return (line.substr(pos+1));
            break;
        }
    }
}

void SU2Reader::ReadFile() {
    ifstream file(filename);
    string line;

    nDim = stoi(FindTag("NDIME", file));
    nElm = stoi(FindTag("NELEM", file));
    elms.resize(nElm);

    unsigned long conn_line[nElm][MAX_LINE_LEN];

    for (auto iLine = 0ul; iLine < nElm; iLine++) {
        getline(file, line);
        istringstream iss(line);
        unsigned long value;
        size_t counter = 0;

        while (iss >> value) conn_line[iLine][counter++] = value;
    }

    vtk.resize(nElm);
    unsigned long count = 0;
    for (auto iLine = 0; iLine < nElm; iLine++) {
        unsigned long *line = conn_line[iLine];
        vtk[iLine] = line[0];

        // triangle
        if (vtk[iLine] == 5) {
            vector<unsigned long> elm_nodes(3);
            elm_nodes[0] = line[1];
            elm_nodes[1] = line[2];
            elm_nodes[2] = line[3];
            elms[iLine].SetNodes(elm_nodes);
        }
        else { throw std::out_of_range("Only triangular elements for now.\n"); }
    }

    nNode = stoi(FindTag("NPOIN", file));
    nodes.resize(nNode);
    for (auto iNode = 0ul; iNode < nNode; iNode++) {
        getline(file, line);
        istringstream iss(line);
        size_t count = 0;
        double value;
        vector<double> xyz(3, 0.0);
        while (iss >> value) xyz[count++] = value;
        Node &n = nodes[iNode];
        n.SetXYZ(xyz[0], xyz[1]);
    }

    nBoundaries = stoi(FindTag("NMARK", file));
    bounds.resize(nBoundaries);

    for (auto iBound = 0; iBound < nBoundaries; iBound++) {
        Bound &bndry = bounds[iBound];
        string bName = FindTag("MARKER_TAG", file);
        bndry.SetName(bName);
        unsigned long nBElms = stoi(FindTag("MARKER_ELEMS", file));
        vector<Edge> bEdges(nBElms);
        for (auto iBElm = 0; iBElm < nBElms; iBElm++) {
            Edge &edge = bEdges[iBElm];
            getline(file, line);
            istringstream iss(line);
            unsigned long bLine[3];
            unsigned long value;
            size_t count = 0;
            while (iss >> value) bLine[count++] = value;
            edge.SetNodes(bLine[1], bLine[2]);
        }

        bndry.SetEdges(bEdges);
    }
}
