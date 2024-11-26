#include "Grid.hpp"
#include <set>

Grid::Grid() : initialized(false) { }
Grid::~Grid() { }

void Grid::Init(unsigned short nDim, 
        vector<Bound> &bounds, vector<Node> &nodes, 
        vector<Elm> &elms, vector<unsigned short> &vtk) {
    this->nDim = nDim;
    swap(this->bounds, bounds);
    swap(this->nodes, nodes);
    swap(this->elms, elms);
    swap(this->vtk, vtk);
    nElm = this->elms.size();
    nNode = this->nodes.size();
    nBound = this->bounds.size();

    cout << "nElms from grid: " << nElm << endl;
}

void Grid::ElmData(void) {
    // Centroids
    for (auto iElm = 0ul; iElm < nElm; iElm++) {
        Elm &elm = elms[iElm];
        const auto &elm_node = elm.GetNodes();
        const unsigned long &n1_ind = elm_node[0]; 
        const unsigned long &n2_ind = elm_node[1]; 
        const unsigned long &n3_ind = elm_node[2]; 

        const Node &n1 = nodes[n1_ind];
        const Node &n2 = nodes[n2_ind];
        const Node &n3 = nodes[n3_ind];
        double x1 = n1.X(); double x2 = n2.X(); double x3 = n3.X();
        double y1 = n1.Y(); double y2 = n2.Y(); double y3 = n3.Y();
        double z1 = n1.Z(); double z2 = n2.Z(); double z3 = n3.Z();

        // midpoint
        double xc = (x1 + x2 + x3) / 3.0;
        double yc = (y1 + y2 + y3) / 3.0;
        double zc = (z1 + z2 + z3) / 3.0;

        double area = 0.5 * (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));

        elm.SetCentroid(xc, yc, zc);
        elm.SetArea(area);
    }

    // Add element indices to nodes
    for (auto iElm = 0ul; iElm < nElm; iElm++) {
        Elm &elm = elms[iElm];
        const auto &elm_node = elm.GetNodes();
        const unsigned long &n1_ind = elm_node[0]; 
        const unsigned long &n2_ind = elm_node[1]; 
        const unsigned long &n3_ind = elm_node[2]; 

        Node &n1 = nodes[n1_ind];
        Node &n2 = nodes[n2_ind];
        Node &n3 = nodes[n3_ind];

        n1.AddElm(iElm);
        n2.AddElm(iElm);
        n3.AddElm(iElm);
    }
}

bool Grid::IsInternalEdge(const unsigned long node_ind1, const unsigned long node_ind2) const {
    vector<unsigned long> shared_elms = GetSharedElements(node_ind1, node_ind2);

    // Nodes will share two elements if edge is internal, one element if boundary edge
    return shared_elms.size() == 2 ? true : false;
}

vector<unsigned long> Grid::GetSharedElements(const unsigned long node_ind1, const unsigned long node_ind2) const {
    std::set<unsigned long> shared;
    const Node &n1 = nodes[node_ind1];
    const Node &n2 = nodes[node_ind2];

    const vector<unsigned long> &n1_elms = n1.GetElms();
    const vector<unsigned long> &n2_elms = n2.GetElms();

    for (const auto n1_elm : n1_elms) {
        for (const auto n2_elm : n2_elms) {
            if(n1_elm == n2_elm) shared.insert(n1_elm);
        }
    }

    vector<unsigned long> shared_elms;
    shared_elms.assign(shared.begin(), shared.end());
    return shared_elms;
}

// definitely a better way to do this
void Grid::InternalEdges(void) {
    // Count unique edges
    nEdge = 0;
    for (auto iElm = 0ul; iElm < nElm; iElm++) {
        Elm &elm = elms[iElm];
        const auto &elm_node = elm.GetNodes();
        const unsigned long &n1_ind = elm_node[0]; 
        const unsigned long &n2_ind = elm_node[1]; 
        const unsigned long &n3_ind = elm_node[2]; 

        if ((n2_ind > n1_ind) && IsInternalEdge(n1_ind, n2_ind)) nEdge++;
        if ((n3_ind > n2_ind) && IsInternalEdge(n2_ind, n3_ind)) nEdge++;
        if ((n1_ind > n3_ind) && IsInternalEdge(n3_ind, n1_ind)) nEdge++;
    }

    // Allocate edges add elm/edge indices to each other 
    edges.resize(nEdge);
    unsigned long ind_edge = 0;
    for (auto iElm = 0ul; iElm < nElm; iElm++) {
        Elm &elm = elms[iElm];
        const auto &elm_node = elm.GetNodes();
        const unsigned long &n1_ind = elm_node[0]; 
        const unsigned long &n2_ind = elm_node[1]; 
        const unsigned long &n3_ind = elm_node[2]; 

        if ((n2_ind > n1_ind) && IsInternalEdge(n1_ind, n2_ind)) {
            Edge &e = edges[ind_edge];
            e.SetNodes(n1_ind, n2_ind);
            auto shared_elms = GetSharedElements(n1_ind, n2_ind);
            //elm.AddEdge(ind_edge, n1_ind, n2_ind);
            for (auto shared_elm : shared_elms) {
                e.AddElm(shared_elm);
            }
            ind_edge++;
        };
        if ((n3_ind > n2_ind) && IsInternalEdge(n2_ind, n3_ind)) {
            Edge &e = edges[ind_edge];
            e.SetNodes(n2_ind, n3_ind);
            auto shared_elms = GetSharedElements(n2_ind, n3_ind);
            //elm.AddEdge(ind_edge, n2_ind, n3_ind);
            for (auto shared_elm : shared_elms) {
                e.AddElm(shared_elm);
            }
            ind_edge++;
        };
        if ((n1_ind > n3_ind) && IsInternalEdge(n3_ind, n1_ind)) {
            Edge &e = edges[ind_edge];
            e.SetNodes(n3_ind, n1_ind);
            auto shared_elms = GetSharedElements(n3_ind, n1_ind);
            //elm.AddEdge(ind_edge, n3_ind, n1_ind);
            for (auto shared_elm : shared_elms) {
                e.AddElm(shared_elm);
            }
            ind_edge++;
        };
    }

    /*
    for (auto iEdge = 0ul; iEdge < nEdge; iEdge++) {
        Edge &edge = edges[iEdge];
        const vector<unsigned long> &edge_elms = edge.GetElms();
        const unsigned long n1 = edge.Node1();
        const unsigned long n2 = edge.Node2();
        for (auto ielm = 0; ielm < edge_elms.size(); ielm++)  {
            elms[edge_elms[ielm]].AddEdge(iEdge, n1, n2);
        }
    }
    */

    ofstream file2("edge_check_3.txt");
    // Geometric data
    for (auto iEdge = 0ul; iEdge < nEdge; iEdge++) {
        Edge &edge = edges[iEdge];
        const unsigned long i_n1 = edge.Node1();
        const unsigned long i_n2 = edge.Node2();
        const Node &n1 = nodes[i_n1];
        const Node &n2 = nodes[i_n2];

        // compute length, norm and midpoint
        double nx = n1.Y() - n2.Y();
        double ny = n2.X() - n1.X();
        double len = sqrt(nx * nx + ny * ny);
        nx /= len;
        ny /= len;
        double mx = 0.5 * (n1.X() + n2.X());
        double my = 0.5 * (n1.Y() + n2.Y());

        edge.SetMidpoint(mx, my);
        edge.SetLength(len);
        edge.SetNorm(nx, ny);

        // Figure out left and right elements of edge
        const vector<unsigned long> &e_elms = edge.GetElms();
        //if (e_elms.size() != 2) cout << "edge has too many/few elms!\n"; exit(0); // commit self delete

        const vector<double> &e1c = elms[e_elms[0]].GetCentroid();

        double e1x = e1c[0] - mx;
        double e1y = e1c[1] - my;
        double dot = nx * e1x + ny * e1y;
        
        // norm should point from left element to right element
        // if dot is positive, normal is pointing towards the first element, make it elmR
        if (dot > 0.0) edge.SetElms(e_elms[1], e_elms[0]);
        else edge.SetElms(e_elms[0], e_elms[1]);

        // set edge numbers, indicates which side (0,1,2) of the left/right element the edge is
        // The Elm class automatically sorts edges, check iEdge against left/right element edges
        // 2  /\  1
        //   /__\  
        //    0
        Elm &elmL = elms[edge.ElmL()];
        Elm &elmR = elms[edge.ElmR()];

        const auto &ln = elmL.GetNodes();
        const auto &rn = elmR.GetNodes();

        /*
           // this data looks good
        file2 << "-----------------------------\n";
        file2 << "Edge nodes: " << edge.Node1() << "\t" << edge.Node2() << endl;
        file2 << "Left elm nodes: " << ln[0] << "\t" << ln[1] << "\t" << ln[2] << endl;
        file2 << "Right elm nodes: " << rn[0] << "\t" << rn[1] << "\t" << rn[2] << endl;
        file2 << "-----------------------------\n";
        */

        elmL.AddEdge(iEdge, i_n1, i_n2);
        elmR.AddEdge(iEdge, i_n1, i_n2);

        const auto edge_numL = elmL.GetEdgeNum(i_n1, i_n2);
        const auto edge_numR = elmR.GetEdgeNum(i_n1, i_n2);

        edge.SetEdgeNumL(edge_numL);
        edge.SetEdgeNumR(edge_numR);

        file2 << "---------------------------------------------\n";
        file2 << "edge_nodes: " << i_n1 << "\t" << i_n2 << endl;
        file2 << "edge_nums: " << edge_numL << "\t" << edge_numR << endl;
        file2 << "eL_nodes: " << ln[0] << "\t" << ln[1] << "\t" << ln[2] << endl;
        file2 << "eR_nodes: " << rn[0] << "\t" << rn[1] << "\t" << rn[2] << endl;

        /*
        if (iEdge == elmL_edges[0]) edge.SetEdgeNumL(0);
        else if (iEdge == elmL_edges[1]) edge.SetEdgeNumL(1);
        else if (iEdge == elmL_edges[2]) edge.SetEdgeNumL(2);
        //else cout << "couldn't find edge on elm!\n"; exit(0); // panic

        if (iEdge == elmR_edges[0]) edge.SetEdgeNumR(0);
        else if (iEdge == elmR_edges[1]) edge.SetEdgeNumR(1);
        else if (iEdge == elmR_edges[2]) edge.SetEdgeNumR(2);
        //else cout << "couldn't find edge on elm!\n"; exit(0);// panic
        */
    }

    file2.close();
  //  exit(0);
}

void Grid::BoundaryData(void) {

    // cyl and naca
    bounds[0].type = Bound::BCType::InviscidWall;
    bounds[1].type = Bound::BCType::Farfield;

    // bump
    /*
    bounds[0].type = Bound::BCType::InviscidWall;
    bounds[1].type = Bound::BCType::Farfield;
    bounds[2].type = Bound::BCType::Farfield;
    bounds[3].type = Bound::BCType::Farfield;
    */

    // su2 ramp
    /*
    bounds[0].type = Bound::BCType::SupersonicInlet;
    bounds[1].type = Bound::BCType::InviscidWall;
    bounds[2].type = Bound::BCType::SupersonicOutlet;
    bounds[3].type = Bound::BCType::InviscidWall;
    */

//    ofstream f("bndry_chk.txt");

    // Add elements to boundary edges
    for (auto iBound = 0ul; iBound < nBound; iBound++) {
        Bound &bound = bounds[iBound];
        vector<Edge> &bedges = bound.GetEdges();
        for (Edge &bedge : bedges) {
            const unsigned long &n1_ind = bedge.Node1();
            const unsigned long &n2_ind = bedge.Node2();
            vector<unsigned long> shared_elms = GetSharedElements(n1_ind, n2_ind);
            if (shared_elms.size() != 1) {
                cout << "bound edge has more than 1 elm!\n";
                exit(0);
            } // self delete
            bedge.AddElm(shared_elms[0]);
            bedge.SetElms(-1, shared_elms[0]); // set singled element to elmR

            const Node &n1 = nodes[n1_ind];
            const Node &n2 = nodes[n2_ind];
            double nx = n2.Y() - n1.Y();
            double ny = n1.X() - n2.X();
            double len = sqrt(nx * nx + ny * ny);
            nx /= len;
            ny /= len;
            double mx = 0.5 * (n1.X() + n2.X());
            double my = 0.5 * (n1.Y() + n2.Y());

            bedge.SetMidpoint(mx, my);
            bedge.SetLength(len);

            // Make sure normals are outward pointing
            const vector<double> &eCtr = elms[shared_elms[0]].GetCentroid();
            double ex = eCtr[0] - mx;
            double ey = eCtr[1] - my;
            double dot = nx * ex + ny * ey;

            // if dot is positive the normal is pointing inward, flip
            if (dot > 0.0) {
                nx *= -1.0;
                ny *= -1.0;
            }

            bedge.SetNorm(nx, ny);

            // set edge numbers, indicates which side (0,1,2) of the left/right element the edge is
            // The Elm class automatically sorts edges, check iEdge against left/right element edges
            // 2  /\  1
            //   /__\  
            //    0
            Elm &elmR = elms[bedge.ElmR()];
            elmR.AddEdge(-1, n1_ind, n2_ind);
            const unsigned short side = elmR.GetEdgeNum(n1_ind, n2_ind);
            bedge.SetEdgeNumR(side);

            const auto ernodes = elmR.GetNodes();


            /*
            f << "-----------------------------------\nedge_nodes: " << n1_ind << "\t" << n2_ind << "\n" << 
                "elm_nodes: " << ernodes[0] << "\t" << ernodes[1] << "\t" << ernodes[2] <<
                "\nelm: " << shared_elms[0] <<  "\tedge_num" << bedge.EdgeNumR() << endl;
                */
        }
    }

    //f.close();
}

void Grid::WriteTestData(const string filename) const {
    std::ofstream file(filename);

    for (auto iElm = 0ul; iElm < nElm; iElm++) {
        const Elm  &elm = elms[iElm];
        const auto elm_edges = elm.Edges();
        file << elm_edges[0] << "\t" << elm_edges[1] << "\t" << elm_edges[2] << endl;
    }

    /*
    for ( auto iB = 0ul; iB < nBound; iB++) {
        file << "********************************************************\n";
        file << "Boundary: " << iB << endl;;
        const vector<Edge> &bedges = bounds[iB].GetEdges();
        for (const auto &bedge : bedges) {
            file << bedge.Node1() << "\t" << bedge.Node2() << "\t";
            const Node &n1 = nodes[bedge.Node1()];
            const Node &n2 = nodes[bedge.Node2()];
            file << n1.X() << "\t" << n1.Y() << "\t" << n2.X() << "\t" << n2.Y() << "\t";
            file << bedge.Norm()[0] << "\t" << bedge.Norm()[1] <<  "\t" << bedge.Length() << endl;
        }
        file << "********************************************************\n";
    }
    */

    file.close();

    ofstream f("edge_testing2.txt");
    for (const auto edge : edges) {
        f << edge.ElmL() << "\t" << edge.ElmR() << "\t" << edge.EdgeNumL() << "\t" << edge.EdgeNumR() << endl;
    }
    f.close();
}

void Grid::BuildGridData(void) {
    cout << "Constructing Element Data...\n";
    ElmData();

    cout << "Constructing Edge Data...\n";
    InternalEdges();

    cout << "Constructing Boundary Data...\n";
    BoundaryData();
}

unsigned long Grid::GetNumNodes(void) { return nNode; }
unsigned short Grid::GetNumDim(void) { return nDim; }
unsigned long Grid::GetNumCells(void) { return nElm; }

const unsigned long Grid::GetNumNodes(void) const { return nNode; }
const unsigned short Grid::GetNumDim(void) const { return nDim; }
const unsigned long Grid::GetNumCells(void) const { return nElm; }

const vector<Node>& Grid::GetNodes(void) const { return nodes; }
const vector<Elm>& Grid::GetElms(void) const { return elms; }
const vector<Edge>& Grid::GetEdges(void) const { return edges; }
const vector<Bound>& Grid::GetBoundaries(void) const { return bounds; }
const vector<unsigned short>& Grid::GetVTK(void) const { return vtk; }

const vector<double>& Grid::GetNodeValues(void) const { return node_values; }
void Grid::FillNodeValues(void) {
    node_values.resize(nNode, 0.0);

    for (auto iBound = 0ul; iBound < nBound; iBound++) {
        Bound &b = bounds[iBound];
        const vector<Edge> &bEdge = b.GetEdges();
        for (auto iEdge = 0ul; iEdge < bEdge.size(); iEdge++) {
            auto n1 = bEdge[iEdge].Node1();
            auto n2 = bEdge[iEdge].Node2();

            node_values[n1] = 100.0;
            node_values[n2] = 100.0;
        }
        
    }
}

void Grid::ComputeJacobiansDG(void) {
    dg_jac.resize(nElm, 0.0);

    for (auto iElm = 0ul; iElm < nElm; iElm++) {
        const Elm &elm = elms[iElm];
        const vector<unsigned long> &iNodes = elm.GetNodes();
        const Node &n1 = nodes[iNodes[0]];
        const Node &n2 = nodes[iNodes[1]];
        const Node &n3 = nodes[iNodes[2]];
        double x1 = n1.X(), x2 = n2.X(), x3 = n3.X();
        double y1 = n1.Y(), y2 = n2.Y(), y3 = n3.Y();
        double z1 = n1.Z(), z2 = n2.Z(), z3 = n3.Z();

        dg_jac[iElm] = (x2-x1) * (y3-y1) - (x3-x1) * (y2-y1);
    }
}

const vector<double>& Grid::GetJacobiansDG(void) const { return dg_jac; }
