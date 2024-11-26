#pragma once

#include <memory>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "Edge.hpp"
#include "Bound.hpp"
#include "Node.hpp"
#include "Elm.hpp"

using namespace std;

class Grid {
    public:
        Grid();
        ~Grid();

        friend class GridReader;

        void Init(unsigned short nDim, vector<Bound> &bounds, 
                vector<Node> &nodes, vector<Elm> &elms, vector<unsigned short> &vtk);
        void BuildGridData(void);
        unsigned long GetNumNodes(void);
        unsigned short GetNumDim(void);
        unsigned long GetNumCells(void);
        unsigned long* GetConnectivity(void);
        unsigned short* GetVTK(void);
        double* GetCoordinates(void);
        const unsigned long GetNumNodes(void) const;
        const unsigned short GetNumDim(void) const;
        const unsigned long GetNumCells(void) const;

        const vector<Node>& GetNodes(void) const;
        const vector<Elm>& GetElms(void) const;
        const vector<Edge>& GetEdges(void) const;
        const vector<Bound>& GetBoundaries(void) const; 
        const vector<unsigned short>& GetVTK(void) const;

        void ComputeJacobiansDG(void);
        const vector<double>& GetJacobiansDG(void) const;


        // Testing
        const vector<double>& GetNodeValues(void) const;
        void FillNodeValues(void);

        void WriteTestData(const string filename) const;

    private:
        bool initialized;
        unsigned long nElm;
        unsigned long nEdge;
        unsigned short nDim;
        unsigned short nBound;
        unsigned long nNode;

        vector<Edge> edges;
        vector<Node> nodes;
        vector<Bound> bounds;
        vector<Elm> elms;
        vector<unsigned short> vtk;
        vector<double> dg_jac;

        void ElmData(void);
        void InternalEdges(void);
        void BoundaryData(void);
        bool IsInternalEdge(const unsigned long node_ind1, const unsigned long node_ind2) const;
        vector<unsigned long> GetSharedElements(const unsigned long node_ind1, const unsigned long node_ind2) const;

        // Testing
        vector<double> node_values;
};
