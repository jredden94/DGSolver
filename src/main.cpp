#include <iostream>
#include "Common/Config.hpp"
#include "Grid/Grid.hpp"
#include "IO/GridReader.hpp"
#include "Solver/Solver.hpp"
#include "IO/VTK.hpp"
#include "Basis/Quad.hpp"
#include "Basis/Basis.hpp"

int main(void) {

    Basis basis;
    Quad quad;

    unsigned short n = 3;

    Config *config = &Config::GetConfig();
    const string grid_file = config->GetGridFilename();

    unique_ptr grid_reader = GridReader::GetReader(grid_file);
    Grid grid;

    cout << "Reading file...\n";
    grid_reader->ReadFile();
    grid_reader->TransferToGrid(grid);

    cout << "Building grid...\n";
    grid.BuildGridData();
    grid.ComputeJacobiansDG();
    grid.WriteTestData("edge_testing.txt");
    //VTK::WriteVTK("out.vtu", grid);
    //return 0;

    unique_ptr<Solver> solver = Solver::CreateSolver(&grid);
    solver->Solve();

    const SolVector &u = solver->GetPrimitives();

    VTK::WriteVTK("out.vtu", grid, u);

    //grid.WriteTestData("edges.txt");
    return 0;
}
