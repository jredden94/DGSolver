#pragma once

#include <fstream>
#include <sstream>
#include <istream>
#include <memory>
#include <exception>
#include <string>
#include <vector>
#include "../Grid/Grid.hpp"
#include "../Grid/Node.hpp"
#include "../Grid/Bound.hpp"

using namespace std;

class GridReader {
    public:
        GridReader();
        virtual ~GridReader();

        static unique_ptr<GridReader> GetReader(const string filename);

        void TransferToGrid(Grid &grid);
        virtual void ReadFile(void) = 0;

    protected:
        string filename;

        unsigned short nDim, nVar;
        unsigned long nNode, nElm, nBoundaries;

        vector<Node> nodes;
        vector<Bound> bounds;
        vector<Elm> elms;
        vector<unsigned short> vtk;

        const unsigned short MAX_LINE_LEN = 15;
        const unsigned short MAX_NODES = 8;
};
