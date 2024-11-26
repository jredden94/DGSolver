#include "GridReader2D.hpp"

GridReader2D::GridReader2D() { }
GridReader2D::GridReader2D( string filename ) { GridReader::filename = filename; }
GridReader2D::~GridReader2D() { }

void GridReader2D::ReadFile() {
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "failed to open file :(\n";
        5 / 0; // panic;
    }

    unsigned long nTri, nQuad;
    file >> nNode >> nTri >> nQuad;

    if (nQuad != 0) {
        cout << "NO QUADS!\n";
        5 / 0; // no quads allowed
    }

    nElm = nTri + nQuad;
    elms.resize(nElm);
    vtk.resize(nElm);
    nodes.resize(nNode);

    for (auto iNode = 0ul; iNode < nNode; iNode++) {
        double x, y, z;
        file >> x >> y;
        nodes[iNode].SetXYZ(x, y, 0.0);
    }

    for (auto iTri = 0ul; iTri < nTri; iTri++) {
        unsigned long n1, n2, n3;
        file >> n1 >> n2 >> n3;
        vector<unsigned long> iNodes {n1-1, n2-1, n3-1};
        elms[iTri].SetNodes(iNodes);
    }

    for (auto iQuad = 0ul; iQuad < nQuad; iQuad--) {
        unsigned long n1, n2, n3, n4;
        file >> n1 >> n2 >> n3 >> n4;
        vector<unsigned long> iNodes {n1-1, n2-1, n3-1, n4-1};
        elms[iQuad + nTri].SetNodes(iNodes);
    }

    file >> nBoundaries;
    bounds.resize(nBoundaries);

    vector<vector<unsigned long>> bnodes;
    bnodes.resize(nBoundaries);
    for (auto iBndry = 0ul; iBndry < nBoundaries; iBndry++) {
        unsigned long nbnodes;
        file >> nbnodes;
        bnodes[iBndry].resize(nbnodes);
    }
    
    for (auto iBndry = 0ul; iBndry < nBoundaries; iBndry++) {
        const unsigned long nbnode = bnodes[iBndry].size();
        for (auto ibNode = 0ul; ibNode < nbnode; ibNode++) {
            unsigned long n1, n2;
            file >> n1;
            //file >> n1 >> n2; // bump.grid
            bnodes[iBndry][ibNode] = n1-1;
        }
    }

    for (auto iB = 0ul; iB < nBoundaries; iB++) {
        Bound &bndry = bounds[iB];
        const auto conn = bnodes[iB];
        const unsigned long nbElm = conn.size() - 1;
        unsigned long counter = 0;

        vector<Edge> bedges(nbElm);
        for (auto iElm = 0ul; iElm < nbElm; iElm++) {
            unsigned long n1, n2;
            n1 = conn[counter++];
            n2 = conn[counter];
            bedges[iElm].SetNodes(n1, n2);
        }

        bndry.SetEdges(bedges);
    }

    /*
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
    */
}
