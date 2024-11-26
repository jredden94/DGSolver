#include "Solver.hpp"
#include <cmath>
#include <fstream>


Euler::Euler() : Solver(Config::SolverType::EulerFVM) { }
Euler::~Euler() { }

void Euler::Solve() { 
    const vector<Elm> &elms = grid->GetElms();
    const vector<Edge> &edges = grid->GetEdges();
    const vector<Bound> &bounds = grid->GetBoundaries();

    cout << "nElms: " << elms.size() << endl;
    cout << "nEdges: " << edges.size() << endl;
    cout << "nBounds: " << bounds.size() << endl;

    unsigned long currentIter = 0;
    unsigned long maxIter = config->GetMaxIterations();

    const double write_increment = 0.02;
    double next_write = 0.02;
    double t = 0.0;
    const double tf = 3.0;

    while ( t < tf) {
        for (auto i = 0ul; i < nEqn; i++) waveSpeed[i] = 0.0;
        residual.SetZeroes();

        unsigned long counter = 0;

        // Internal Fluxes
        for (const Edge &edge : edges) {
            const unsigned long elmL = edge.ElmL();
            const unsigned long elmR = edge.ElmR();
            const double area = edge.Length();
            const vector<double> &norm = edge.Norm();
            const double *vL = primVar.GetBlock(elmL);
            const double *vR = primVar.GetBlock(elmR);

            /*
            cout << "-------------------------\n";
            cout << elmL << "\t" << elmR << "\t" << area << endl;
            cout << norm[0] << "\t" << norm[1] << "\t" << norm[2] << endl;
            cout << vL[0] << "\t" << vL[1] << "\t" << vL[2] << "\t" << vL[3] << endl;
            cout << vR[0] << "\t" << vR[1] << "\t" << vR[2] << "\t" << vR[3] << endl;
            if (counter++ > 8) break;
            */

            conv_flux->SetStates(vL, vR, norm.data(), area );
            conv_flux->ComputeFlux();

            /*
            const vector<double> &flux = conv_flux->Flux();
            cout << flux[0] << "\t" << flux[1] << "\t" << flux[2] << "\t" << flux[3] << endl;
            cout << "-------------------------\n";
            */

            waveSpeed[elmL] += conv_flux->MaxWaveSpeed();
            waveSpeed[elmR] += conv_flux->MaxWaveSpeed();
            residual.AddBlock(elmR, conv_flux->Flux());
            residual.SubtractBlock(elmL, conv_flux->Flux());

        }

        // Boundary Fluxes
        for (const Bound &bound : bounds) {
            switch (bound.Type()) {
                case(Bound::BCType::Farfield) : FarField(bound); break;
                case(Bound::BCType::InviscidWall) : InviscidWall(bound); break;
                case(Bound::BCType::ViscousWall) : ViscousWall(bound); break;
                case(Bound::BCType::Inlet) : Inlet(bound); break;
                case(Bound::BCType::Outlet) : Outlet(bound); break;
                case(Bound::BCType::SupersonicInlet) : SupersonicInlet(bound); break;
                case(Bound::BCType::SupersonicOutlet) : SupersonicOutlet(bound); break;
                case(Bound::BCType::Empty) : 5 / 0; break;
            }
        }

        double dt_min = MAXFLOAT;
        for ( auto iElm = 0ul; iElm < grid->GetNumCells(); iElm++) {
            double ws = waveSpeed[iElm];
            const Elm &elm = elms[iElm];
            double vol = elm.Area();
            double dt = cfl * vol / waveSpeed[iElm];
            if (dt < dt_min) dt_min = dt;
        }

        if (t + dt_min > next_write) dt_min = next_write - t;


        for ( auto iElm = 0ul; iElm < grid->GetNumCells(); iElm++) {
            double ws = waveSpeed[iElm];
            const Elm &elm = elms[iElm];
            double vol = elm.Area();
            const double *res = residual.GetBlock(iElm);
            conVar.AddBlock(iElm, res, dt_min/vol);
        }

        t += dt_min;

        UpdatePrimitiveVars();
        PrintResiduals(currentIter++, t);


        if ( abs(t - next_write) < 1e-12) {
            next_write += write_increment;
            string filename = "out_" + to_string(t) + ".vtu";
            VTK::WriteVTK(filename, *grid, primVar);
        }
    }
}


/* Note: Need to add variables to boundaries, currently always assumes freestream values at inlets, etc. */
void Euler::BoundaryFluxes(const Bound &boundary) {
    switch (boundary.Type() ) {
        case (Bound::BCType::Inlet) : Inlet(boundary); break;
        case (Bound::BCType::Outlet) : Outlet(boundary); break;
        case (Bound::BCType::InviscidWall) : InviscidWall(boundary); break;
        case (Bound::BCType::ViscousWall) : ViscousWall(boundary); break;
        case (Bound::BCType::Farfield) : FarField(boundary); break;
        case (Bound::BCType::SupersonicInlet) : SupersonicInlet(boundary); break;
        case (Bound::BCType::SupersonicOutlet) : SupersonicOutlet(boundary); break;
        case (Bound::BCType::Empty) : break;
    }
}

void Euler::InviscidWall(const Bound &wall) { 
    const vector<Edge> &edges = wall.GetEdges();
    const vector<double> &u_inf = config->GetFreestreamVelocity();
    double v_j[nVar];
    for (const Edge &edge : edges) {
        unsigned long elm = edge.GetElms()[0];

        const double *v_i = primVar.GetBlock(elm);
        const vector<double> norm = edge.Norm();
        const double area = edge.Length();


        v_j[0] = v_i[0];
        double proj_vel = 0.0;
        for (auto iDim = 0ul; iDim < nDim; iDim++) {
            proj_vel += v_i[iDim+1] * norm[iDim];
            v_j[iDim+1] = v_i[iDim+1];
        }

        for (auto iDim = 0ul; iDim < nDim; iDim++) {
            v_j[iDim+1] -= proj_vel * norm[iDim];
        }

        v_j[nVar-1] = v_i[nVar-1];

        conv_flux->SetStates(v_i, v_j, norm.data(), area);
        conv_flux->ComputeFlux();
        waveSpeed[elm] += conv_flux->MaxWaveSpeed();
        residual.SubtractBlock(elm, conv_flux->Flux());
    }
    /*
       const vector<unsigned long> &iNodes = wall.Nodes();
       const vector<vector<zdouble>> &norms = wall.NodeNorms();
       const vector<zdouble> &dual_areas = wall.DualAreas();

       for (auto i = 0ul; i < iNodes.size(); i++) {
       const unsigned long &iNode = iNodes[i];
       const vector<zdouble> &norm = norms[i];
       const zdouble &area = dual_areas[i];
       const zdouble *v_domain = primVar.GetBlock(iNode);

       v_j[0] = v_domain[0];
       zdouble proj_vel = 0.0;
       for (auto iDim = 0ul; iDim < nDim; iDim++) {
       proj_vel += v_domain[iDim+1] * norm[iDim];
       v_j[iDim+1] = v_domain[iDim+1];
       }


       for (auto iDim = 0ul; iDim < nDim; iDim++) {
       v_j[iDim+1] -= proj_vel * norm[iDim];
       }
       v_j[nVar-1] = v_domain[nVar-1];

       conv_flux->SetStates(v_domain, v_j, norm.data(), area);
       conv_flux->ComputeFlux();

       waveSpeed[iNode] += conv_flux->MaxWaveSpeed();
       residual.SubtractBlock(iNode, conv_flux->Flux());
       if (implicit) jacMat.AddBlock(iNode, iNode, conv_flux->JacI());
       }
       */
}

void Euler::ViscousWall(const Bound &wall) { 
    /*
       const vector<unsigned long> &iNodes = wall->Nodes();
       const vector<vector<zdouble>> &norms = wall->NodeNorms();
       const vector<zdouble> &area = wall->DualAreas();

    for (auto i = 0ul; i < iNodes.size(); i++) {
        const unsigned long &iNode = iNodes[i];
        zdouble *u_i = conVar.GetBlock(iNode);
        zdouble *res = residual.GetBlock(iNode);

        // Set momentums and their residual contributions to 0
        for (size_t iDim = 0; iDim < nDim; iDim++) {
            u_i[iDim+1] = 0.0;
            res[iDim+1] = 0.0;
        }

        // Set momentum contributions to jacobian (middle rows) to 0, diagonals to 1
        if (implicit) {
            zdouble *jac = jacMat.GetBlock(iNode, iNode);
            for (size_t iVar = 1; iVar < nVar-1; iVar++) {
                for (size_t jVar = 0; jVar < nVar; jVar++)
                    jac[iVar * nVar + jVar] = iVar == jVar ? 1.0 : 0.0;
            }
        }
    }
    */
}

void Euler::FarField(const Bound &far) { 
    const vector<Edge> &edges = far.GetEdges();
    const vector<double> &u_inf = config->GetFreestreamVelocity();
    double p_inf = config->GetFreestreamPressure();
    double rho_inf = config->GetFreestreamDensity();
    for (const Edge &edge : edges) {
        unsigned long elm = edge.GetElms()[0];

        const double *v_i = primVar.GetBlock(elm);
        const vector<double> norm = edge.Norm();
        const double area = edge.Length();

        double v_j[nVar];
        v_j[0] = rho_inf;
        for (size_t iDim = 0; iDim < nDim; iDim++) v_j[iDim+1] = u_inf[iDim];
        v_j[nVar-1] = p_inf;

        conv_flux->SetStates(v_i, v_j, norm.data(), area);
        conv_flux->ComputeFlux();
        waveSpeed[elm] += conv_flux->MaxWaveSpeed();
        residual.SubtractBlock(elm, conv_flux->Flux());

    /*
    const vector<unsigned long> &iNodes = far->Nodes();
    const vector<vector<zdouble>> &norms = far->NodeNorms();
    const vector<zdouble> &areas = far->DualAreas();
    const vector<zdouble> vel_inf = config->GetFreestreamVelocity();
    const zdouble rho_inf = config->GetFreestreamDensity();
    const zdouble p_inf = config->GetFreestreamPressure();

    for (auto i = 0ul; i < iNodes.size(); i++) {
        const unsigned long iNode = iNodes[i];
        const zdouble *v_domain = primVar.GetBlock(iNode);

        v_j[0] = rho_inf;
        for (auto iDim = 0ul; iDim < nDim; iDim++)
            v_j[iDim+1] = vel_inf[iDim];
        v_j[nVar-1] = p_inf;

        conv_flux->SetStates(v_domain, v_j, norms[i].data(), areas[i]);
        conv_flux->ComputeFlux();

        waveSpeed[iNode] += conv_flux->MaxWaveSpeed();
        residual.SubtractBlock(iNode, conv_flux->Flux());
        if (implicit) jacMat.AddBlock(iNode, iNode, conv_flux->JacI());
    }
    */
}
}

void Euler::Inlet(const Bound &inlet) { 
    /*
    const vector<unsigned long> &iNodes = inlet->Nodes();
    const vector<vector<zdouble>> &norms = inlet->NodeNorms();
    const vector<zdouble> &area = inlet->DualAreas();
    const vector<zdouble> vel_inf = config->GetFreestreamVelocity();
    const zdouble rho_inf = config->GetFreestreamDensity();
    const zdouble p_inf = config->GetFreestreamPressure();

    for (auto i = 0ul; i < iNodes.size(); i++) {
        const unsigned long iNode = iNodes[i];
        const zdouble *v_domain = primVar.GetBlock(iNode);

        v_j[0] = rho_inf;
        for (size_t iDim = 0; iDim < nDim; iDim++) v_j[iDim+1] = vel_inf[iDim];
        v_j[nVar-1] = v_domain[nVar-1];

        conv_flux->SetStates(v_domain, v_j, norms[i].data(), area[i]);
        conv_flux->ComputeFlux();

        waveSpeed[iNode] += conv_flux->MaxWaveSpeed();
        residual.SubtractBlock(iNode, conv_flux->Flux());
        if (implicit) jacMat.AddBlock(iNode, iNode, conv_flux->JacI());
    }
    */
}

void Euler::Outlet(const Bound &outlet) { 
    /*
    const vector<unsigned long> &iNodes = outlet->Nodes();
    const vector<vector<zdouble>> &norms = outlet->NodeNorms();
    const vector<zdouble> &area = outlet->DualAreas();
    const vector<zdouble> &vel_inf = config->GetFreestreamVelocity();
    const zdouble &rho_inf = config->GetFreestreamDensity();
    const zdouble &p_inf = config->GetFreestreamPressure();

    for (auto i = 0ul; i < iNodes.size(); i++) {
        const unsigned long iNode = iNodes[i];
        const zdouble *v_domain = primVar.GetBlock(iNode);
        zdouble v_2 = 0.0;
        for (auto iDim = 0ul; iDim < nDim; iDim++) 
            v_2 += v_domain[iDim+1] * v_domain[iDim+1];

        const zdouble soundspeed = sqrt(gamma * v_domain[nVar-1] / v_domain[0]);
        const zdouble mach = sqrt(v_2) / soundspeed;

        v_j[0] = rho_inf;
        for (size_t iDim = 0; iDim < nDim; iDim++) v_j[iDim+1] = vel_inf[iDim];

        if (mach >= 1.0) v_j[nVar-1] = v_domain[nVar-1];
        else v_j[nVar-1] = p_inf;

        conv_flux->SetStates(v_domain, v_j, norms[i].data(), area[i]);
        conv_flux->ComputeFlux();

        waveSpeed[iNode] += conv_flux->MaxWaveSpeed();
        residual.SubtractBlock(iNode, conv_flux->Flux());
        if (implicit) jacMat.AddBlock(iNode, iNode, conv_flux->JacI());
    }
    */
}

void Euler::SupersonicInlet(const Bound &inlet) { 
    const vector<Edge> &edges = inlet.GetEdges();
    const vector<double> &u_inf = config->GetFreestreamVelocity();
    double p_inf = config->GetFreestreamPressure();
    double rho_inf = config->GetFreestreamDensity();
    for (const Edge &edge : edges) {
        unsigned long elm = edge.GetElms()[0];

        const double *v_i = primVar.GetBlock(elm);
        const vector<double> norm = edge.Norm();
        const double area = edge.Length();


        double v_j[nVar];
        v_j[0] = rho_inf;
        for (auto iDim = 0; iDim < nDim; iDim++) v_j[iDim+1] = u_inf[iDim];
        v_j[nVar-1] = p_inf;

        conv_flux->SetStates(v_i, v_j, norm.data(), area);
        conv_flux->ComputeFlux();
        waveSpeed[elm] += conv_flux->MaxWaveSpeed();
        residual.SubtractBlock(elm, conv_flux->Flux());

    }
    /*
    const vector<unsigned long> &iNodes = inlet->Nodes();
    const vector<vector<zdouble>> &norms = inlet->NodeNorms();
    const vector<zdouble> &area = inlet->DualAreas();
    const vector<zdouble> vel_inf = config->GetFreestreamVelocity();
    const zdouble rho_inf = config->GetFreestreamDensity();
    const zdouble p_inf = config->GetFreestreamPressure();

    for (auto i = 0ul; i < iNodes.size(); i++) {
        const unsigned long iNode = iNodes[i];
        const zdouble *v_domain = primVar.GetBlock(iNode);

        v_j[0] = rho_inf;
        for (size_t iDim = 0; iDim < nDim; iDim++) v_j[iDim+1] = vel_inf[iDim];
        v_j[nVar-1] = p_inf;

        conv_flux->SetStates(v_domain, v_j, norms[i].data(), area[i]);
        conv_flux->ComputeFlux();

        waveSpeed[iNode] += conv_flux->MaxWaveSpeed();
        residual.SubtractBlock(iNode, conv_flux->Flux());
        if (implicit) jacMat.AddBlock(iNode, iNode, conv_flux->JacI());
    }
    */
}

void Euler::SupersonicOutlet(const Bound &outlet) { 
    const vector<Edge> &edges = outlet.GetEdges();
    const vector<double> &u_inf = config->GetFreestreamVelocity();
    double p_inf = config->GetFreestreamPressure();
    double rho_inf = config->GetFreestreamDensity();
    for (const Edge &edge : edges) {
        unsigned long elm = edge.GetElms()[0];

        const double *v_i = primVar.GetBlock(elm);
        const vector<double> norm = edge.Norm();
        const double area = edge.Length();

        double v_j[nVar];
        for (auto iVar = 0ul ;iVar < nVar; iVar++) v_j[iVar] = v_i[iVar];

        conv_flux->SetStates(v_i, v_j, norm.data(), area);
        conv_flux->ComputeFlux();
        waveSpeed[elm] += conv_flux->MaxWaveSpeed();
        residual.SubtractBlock(elm, conv_flux->Flux());

    }
    /*
    const vector<unsigned long> &iNodes = outlet->Nodes();
    const vector<vector<zdouble>> &norms = outlet->NodeNorms();
    const vector<zdouble> &area = outlet->DualAreas();

    for (auto i = 0ul; i < iNodes.size(); i++) {
        const unsigned long iNode = iNodes[i];
        const zdouble *v_domain = primVar.GetBlock(iNode);
        primVar.CopyBlock(iNode, v_j);

        conv_flux->SetStates(v_domain, v_j, norms[i].data(), area[i]);
        conv_flux->ComputeFlux();

        waveSpeed[iNode] += conv_flux->MaxWaveSpeed();
        residual.SubtractBlock(iNode, conv_flux->Flux());
        if (implicit) jacMat.AddBlock(iNode, iNode, conv_flux->JacI());
    }
    */
}
