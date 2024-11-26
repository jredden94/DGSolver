#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include "../Grid/Grid.hpp"
#include "../LinAlg/SolVector.hpp"

namespace VTK {

    using namespace std;

    void WriteVTK(const string &filename, const Grid &grid, const SolVector &var);
    void WriteVTK(const string &filename, const Grid &grid);

    int GetVTKOffset(const int vtk);
    void InsertIndent(ofstream &file, size_t nIndents = 1);
    string Indent(size_t nIndents = 1);
}
