#include "GridReader.hpp"
#include "SU2Reader.hpp"
#include "GridReader2D.hpp"
#include "MshReader.hpp"

GridReader::GridReader() { }
GridReader::~GridReader() { }

unique_ptr<GridReader> GridReader::GetReader(const string filename) { 
    size_t ext_pos = filename.find('.');
    string extension = filename.substr(ext_pos);

    if (extension == ".su2") return make_unique<SU2Reader>(filename);
    else if (extension == ".grid") return make_unique<GridReader2D>(filename);
    else if (extension == ".msh") return make_unique<MshReader>(filename);
    else return nullptr;
}

void GridReader::ReadFile(void) { }
void GridReader::TransferToGrid(Grid &grid) {
    grid.Init(nDim, bounds, nodes, elms, vtk);
}
