#include "VTK.hpp"

void VTK::WriteVTK(const string &filename, const Grid &grid, const SolVector &val) {
    ofstream f(filename);

    unsigned long nNodes = grid.GetNumNodes();
    unsigned long nCells = grid.GetNumCells();
    unsigned short nDim = grid.GetNumDim();

    const vector<Elm> &elms = grid.GetElms();
    const vector<Node> &nodes = grid.GetNodes();
    const vector<unsigned short> &vtk = grid.GetVTK();

    cout << "hello\n";
    // Testing
    const vector<double> &node_values = grid.GetNodeValues();

    f << "<VTKFile type=\"UnstructuredGrid\" version=\"2.2\" byte_order=\"LittleEndian\">\n";
    f << Indent() << "<UnstructuredGrid>\n"; 
    f << Indent(2) << "<Piece NumberOfPoints=" << "\"" << nNodes << "\"" << " NumberOfCells=" 
        << "\"" << nCells << "\">\n"; 

    f << Indent(3) << "<Points>\n";
    f << Indent(4) << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (auto iNode = 0ul; iNode < nNodes; iNode++) {
        const Node &node = nodes[iNode];
        //f << Indent(5) << coord[iNode * nDim] << " " << coord[iNode * nDim + 1] << " " << 0.0 << endl;
        f << Indent(5) << node.X() << " " << node.Y() << " " << node.Z() << endl;
    }
    f << Indent(4) << "</DataArray>\n";
    f << Indent(3) << "</Points>\n";

    // Testing
    // Dummy data to check boundaries
    /*
    f << Indent(3) << "<PointData>\n";
    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"test\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
        //f << Indent(5) << sol[i] << "\n";
        f << Indent(5) << node_values[i] << endl;
    }
    f << Indent(4) << "</DataArray>\n";
    f << Indent(3) << "</PointData>\n";
    */

    /*
    f << Indent(3) << "<PointData>\n";
    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"rho\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
        f << Indent(5) << w[i][0] << endl;
    }
    f << Indent(4) << "</DataArray>\n";

    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"U\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
            for (unsigned short j = 0; j < nDim; j++) {
                f << Indent(5) << w[i][j+1] << " ";
            }
            if (nDim == 2) f << Indent(5) << 0 << endl;
            else f << endl;
    }
    f << Indent(4) << "</DataArray>\n";

    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"p\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
        f << Indent(5) << w[i][nVar-1] << endl;
    }
    f << Indent(4) << "</DataArray>\n";

    f << Indent(3) << "</PointData>\n";
    */

    f << Indent(3) <<  "<Cells>\n";
    f << Indent(4) << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (auto iCell = 0; iCell < nCells; iCell++) {
        const Elm &elm = elms[iCell];
        const auto elm_nodes = elm.GetNodes();
        f << Indent(5);
        //f << conn[iCell * 3] << " " << conn[iCell * 3 + 1] << " " << conn[iCell * 3 + 2] << endl;
        f <<  elm_nodes[0] << " " << elm_nodes[1] << " " << elm_nodes[2] << endl;
    }

    f << Indent(4) << "</DataArray>\n";
    f << Indent(4) << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    f << Indent(5);
    for (auto iCell = 0; iCell < nCells; iCell++) {
        offset += 3;
        f << offset << " ";
    }
    f << endl;

    f << Indent(4) << "</DataArray>\n";

    f << Indent(4) << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
    f << Indent(5);
    for (auto iCell = 0; iCell < nCells; iCell++) {
        f << 5 << " ";
    }
    f << endl;

    f << Indent(4) << "</DataArray>\n";
    f << Indent(3) << "</Cells>\n";
    f << Indent(3) << "<CellData>\n";
    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"wow\" NumberOfComponents=\"5\" format=\"ascii\">\n";
    for (size_t iCell = 0; iCell < nCells; iCell++) {
        const double *vel = val.GetBlock(iCell);
        f << Indent(5) << vel[0] << " " << vel[1] << " " << vel[2] << " " << 0.0 << " " << vel[3] << endl;
    }
    f << Indent(4) << "</DataArray>\n";
    f << Indent(3) << "</CellData>\n";
    f << Indent(2) << "</Piece>\n";
    f << Indent() << "</UnstructuredGrid>\n";
    f << "</VTKFile>\n";

    f.close();
}

void VTK::WriteVTK(const string &filename, const Grid &grid) {
    ofstream f(filename);

    unsigned long nNodes = grid.GetNumNodes();
    unsigned long nCells = grid.GetNumCells();
    unsigned short nDim = grid.GetNumDim();

    const vector<Elm> &elms = grid.GetElms();
    const vector<Node> &nodes = grid.GetNodes();
    const vector<unsigned short> &vtk = grid.GetVTK();

    vector<double> nval(nodes.size(), 0.0);
    const vector<Bound> &bnds = grid.GetBoundaries();

    const Bound &b1 = bnds[0];
    const auto &b1edges = b1.GetEdges();
    for (const auto b1edge : b1edges) {
        nval[b1edge.Node1()] = -100.0;
        nval[b1edge.Node2()] = -100.0;
    }

    const Bound &b2 = bnds[1];
    const auto &b2edges = b2.GetEdges();
    for (const auto b2edge : b2edges) {
        nval[b2edge.Node1()] = 100.0;
        nval[b2edge.Node2()] = 100.0;
    }



    // Testing
    f << "<VTKFile type=\"UnstructuredGrid\" version=\"2.2\" byte_order=\"LittleEndian\">\n";
    f << Indent() << "<UnstructuredGrid>\n"; 
    f << Indent(2) << "<Piece NumberOfPoints=" << "\"" << nNodes << "\"" << " NumberOfCells=" 
        << "\"" << nCells << "\">\n"; 

    f << Indent(3) << "<Points>\n";
    f << Indent(4) << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (auto iNode = 0ul; iNode < nNodes; iNode++) {
        const Node &node = nodes[iNode];
        //f << Indent(5) << coord[iNode * nDim] << " " << coord[iNode * nDim + 1] << " " << 0.0 << endl;
        f << Indent(5) << node.X() << " " << node.Y() << " " << node.Z() << endl;
    }
    f << Indent(4) << "</DataArray>\n";
    f << Indent(3) << "</Points>\n";

    // Testing
    // Dummy data to check boundaries
    f << Indent(3) << "<PointData>\n";
    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"test\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
        //f << Indent(5) << sol[i] << "\n";
        f << Indent(5) << nval[i] << endl;
    }
    f << Indent(4) << "</DataArray>\n";
    f << Indent(3) << "</PointData>\n";

    /*
    f << Indent(3) << "<PointData>\n";
    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"rho\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
        f << Indent(5) << w[i][0] << endl;
    }
    f << Indent(4) << "</DataArray>\n";

    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"U\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
            for (unsigned short j = 0; j < nDim; j++) {
                f << Indent(5) << w[i][j+1] << " ";
            }
            if (nDim == 2) f << Indent(5) << 0 << endl;
            else f << endl;
    }
    f << Indent(4) << "</DataArray>\n";

    f << Indent(4) << "<DataArray type=\"Float64\" Name=\"p\" format=\"ascii\">\n";
    for (size_t i = 0; i < nNodes; i++) {
        f << Indent(5) << w[i][nVar-1] << endl;
    }
    f << Indent(4) << "</DataArray>\n";

    f << Indent(3) << "</PointData>\n";
    */

    f << Indent(3) <<  "<Cells>\n";
    f << Indent(4) << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (auto iCell = 0; iCell < nCells; iCell++) {
        const Elm &elm = elms[iCell];
        const auto elm_nodes = elm.GetNodes();
        f << Indent(5);
        //f << conn[iCell * 3] << " " << conn[iCell * 3 + 1] << " " << conn[iCell * 3 + 2] << endl;
        f <<  elm_nodes[0] << " " << elm_nodes[1] << " " << elm_nodes[2] << endl;
    }

    f << Indent(4) << "</DataArray>\n";
    f << Indent(4) << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    int offset = 0;
    f << Indent(5);
    for (auto iCell = 0; iCell < nCells; iCell++) {
        offset += 3;
        f << offset << " ";
    }
    f << endl;

    f << Indent(4) << "</DataArray>\n";

    f << Indent(4) << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
    f << Indent(5);
    for (auto iCell = 0; iCell < nCells; iCell++) {
        f << 5 << " ";
    }
    f << endl;

    f << Indent(4) << "</DataArray>\n";
    f << Indent(3) << "</Cells>\n";
    f << Indent(2) << "</Piece>\n";
    f << Indent() << "</UnstructuredGrid>\n";
    f << "</VTKFile>\n";

    f.close();
}

void VTK::InsertIndent(ofstream &file, size_t nIndent) {
    string indent(nIndent * 2, ' ');
    file << indent;
}

std::string VTK::Indent(size_t nIndent) {
    return string(nIndent * 2, ' ');
}

int VTK::GetVTKOffset(const int vtkType) {
    switch (vtkType) {
        case (3) : return 2; break;
        case (5) : return 3; break;
        case (9) : return 4; break;
        case (10) : return 4; break;
        case (12) : return 8; break;
        default : 
            cout << "Unsupported VTK cell type " << vtkType << endl; 
            return -1; 
    }
}
