#include "Solver.hpp"
#include <cstdlib>
#include <fstream>
#include <string>

EulerDG::EulerDG() : Solver(Config::SolverType::EulerDG) { }
EulerDG::~EulerDG() { }

void EulerDG::InitEulerDG() {

    U_coeff.Init(nEqn, n_p * nVar);
    res_dg.Init(nEqn, n_p * nVar);
    
    xr.resize(nEqn, 0.0);
    xs.resize(nEqn, 0.0);
    yr.resize(nEqn, 0.0);
    ys.resize(nEqn, 0.0);
    radii.resize(nEqn, 0.0);
    quad.GetQuadrature(r1, r2, w, s_r, w1d);
    basis.Preval(r1, r2, w, s_r, w1d);
    InitCoefficients();
    ComputeJacobians();
    InscribedCircles();
}

void EulerDG::EvalU0(double *U, int i, unsigned long iElm) {
    U[0] = 0.0;
    U[1] = 0.0;
    U[2] = 0.0;
    U[3] = 0.0;

    const double *u0 = conVar.GetBlock(iElm);

    for (auto j = 0ul; j < n_quad; j++) {
        U[0] += w[j] * u0[0] * basis.Phi(i, j);
        U[1] += w[j] * u0[1] * basis.Phi(i, j);
        U[2] += w[j] * u0[2] * basis.Phi(i, j);
        U[3] += w[j] * u0[3] * basis.Phi(i, j);
    }
}

void EulerDG::InitCoefficients() {
    double U[4];
    auto nElm = grid->GetNumCells();

    for (auto iElm = 0; iElm < nEqn; iElm++) {
        double *U_i = U_coeff.GetBlock(iElm);
        for (auto i_p = 0ul; i_p < n_p; i_p++) {
            EvalU0(U, i_p, iElm);

            for (auto n = 0; n < nVar; n++){
                U_i[n * n_p + i_p] = U[n];
            }
        }
    }
}

void EulerDG::UpdateConVars(void) {
    auto nElm = grid->GetNumCells();
    double U[nVar];
    for (auto iElm = 0ul; iElm < nElm; iElm++) {
        const double *c = U_coeff.GetBlock(iElm);

        for (auto iQuad = 0ul; iQuad < n_quad; iQuad++) {
            U[0] = 0.0;
            U[1] = 0.0;
            U[2] = 0.0;
            U[3] = 0.0;

            for (auto k = 0; k < n_p; k++) {
                U[0] += c[0*n_p + k] * basis.Phi(k, iQuad);
                U[1] += c[1*n_p + k] * basis.Phi(k, iQuad);
                U[2] += c[2*n_p + k] * basis.Phi(k, iQuad);
                U[3] += c[3*n_p + k] * basis.Phi(k, iQuad);
            }
        }

        conVar.SetBlock(iElm, U);
    }
}

void EulerDG::VolumeIntegral(const SolVector &Uc) {
    double U[nVar];
    double f_x[nVar], f_y[nVar];
    double vol_int[n_p * nVar];

    for (auto iElm = 0ul; iElm < nEqn; iElm++) {
        const double *c = Uc.GetBlock(iElm);
        for (auto i = 0; i < n_p; i++) {
            vol_int[0*n_p + i] = 0.0;
            vol_int[1*n_p + i] = 0.0;
            vol_int[2*n_p + i] = 0.0;
            vol_int[3*n_p + i] = 0.0;
        }

        const double x_r = xr[iElm];
        const double x_s = xs[iElm];
        const double y_r = yr[iElm];
        const double y_s = ys[iElm];

        for (auto iQuad = 0ul; iQuad < n_quad; iQuad++) {

            U[0] = 0.0;
            U[1] = 0.0;
            U[2] = 0.0;
            U[3] = 0.0;

            for (auto i_p = 0ul; i_p < n_p; i_p++) {
                U[0] += c[0*n_p + i_p] * basis.Phi(i_p, iQuad);
                U[1] += c[1*n_p + i_p] * basis.Phi(i_p, iQuad);
                U[2] += c[2*n_p + i_p] * basis.Phi(i_p, iQuad);
                U[3] += c[3*n_p + i_p] * basis.Phi(i_p, iQuad);
            }

            const double rho = U[0];
            const double rhou = U[1];
            const double rhov = U[2];
            const double E = U[3];
            const double u = rhou / rho;
            const double v = rhov / rho;
            const double vel_sqr = u * u + v * v;
            const double p = (E - 0.5 * rho * (u*u + v*v)) * (gamma - 1);

            // fluxes
            f_x[0] = rhou;
            f_x[1] = rhou * rhou / rho + p;
            f_x[2] = rhou * rhov / rho;
            f_x[3] = rhou * (E + p) / rho;

            f_y[0] = rhov;
            f_y[1] = rhou * rhov / rho;
            f_y[2] = rhov * rhov / rho + p;
            f_y[3] = rhov * (E + p) / rho;

            for (auto i_p = 0ul; i_p < n_p; i_p++) {
                for (auto iVar = 0ul; iVar < nVar; iVar++) {

                    vol_int[iVar*n_p + i_p] += 
                        f_x[iVar] * ( basis.PhiGradX(i_p, iQuad) * y_s
                                - basis.PhiGradY(i_p, iQuad) * y_r)
                        + f_y[iVar] * (- basis.PhiGradX(i_p, iQuad) * x_s 
                                + basis.PhiGradY(i_p, iQuad) * x_r);
                }
            }
        }

        res_dg.AddBlock(iElm, vol_int);
    }
}

void EulerDG::SurfaceIntegral(const SolVector &Uc) {

    const vector<Edge> &edges = grid->GetEdges();
    const vector<Elm> &elms = grid->GetElms();
    double uL[nVar], uR[nVar], vL[nVar], vR[nVar];
    double surf_intL[n_p * nVar];
    double surf_intR[n_p * nVar];
    for (auto iEdge = 0ul; iEdge < edges.size(); iEdge++) {
        const Edge &edge = edges[iEdge];
        const auto iElmL = edge.ElmL();
        const auto iElmR = edge.ElmR();
        const vector<double> &norm = edge.Norm();
        const double len = edge.Length();
        const auto sideL = edge.EdgeNumL();
        const auto sideR = edge.EdgeNumR();

        // Get coefficients, set surface integrals to 0 
        const double *cL = Uc.GetBlock(iElmL);
        const double *cR = Uc.GetBlock(iElmR);
        for (auto i_p = 0ul; i_p < n_p; i_p++) {
            surf_intL[0*n_p + i_p] = 0.0;
            surf_intL[1*n_p + i_p] = 0.0;
            surf_intL[2*n_p + i_p] = 0.0;
            surf_intL[3*n_p + i_p] = 0.0;
            surf_intR[0*n_p + i_p] = 0.0;
            surf_intR[1*n_p + i_p] = 0.0;
            surf_intR[2*n_p + i_p] = 0.0;
            surf_intR[3*n_p + i_p] = 0.0;
        }

        // \int_{\partial \Omega} J F psi
        for (auto iQuad = 0ul; iQuad < n_quad_1d; iQuad++) {
            EvalLeftRight(cL, cR, uL, uR, iQuad, iElmL, iElmR, sideL, sideR); 
            ConsToPrim(uR, vR);
            ConsToPrim(uL, vL);
            conv_flux->SetStates(vL, vR, norm.data(), len);
            conv_flux->ComputeFlux();
            const vector<double> &flux = conv_flux->Flux();

            // NOTE: need to multiply by line element Jacobian J = len/2 since evaluating in reference space. 
            // Flux is already scaled by len so only need to halve the terms below
            for (auto i_p = 0ul; i_p < n_p; i_p++) {
                surf_intL[0*n_p + i_p] += -0.5 * w1d[iQuad] * flux[0] * basis.PhiSide(sideL, i_p, iQuad);
                surf_intL[1*n_p + i_p] += -0.5 * w1d[iQuad] * flux[1] * basis.PhiSide(sideL, i_p, iQuad);
                surf_intL[2*n_p + i_p] += -0.5 * w1d[iQuad] * flux[2] * basis.PhiSide(sideL, i_p, iQuad);
                surf_intL[3*n_p + i_p] += -0.5 * w1d[iQuad] * flux[3] * basis.PhiSide(sideL, i_p, iQuad);
                surf_intR[0*n_p + i_p] += 0.5 * w1d[iQuad] * flux[0] * basis.PhiSide(sideR, i_p, n_quad_1d - 1 - iQuad);
                surf_intR[1*n_p + i_p] += 0.5 * w1d[iQuad] * flux[1] * basis.PhiSide(sideR, i_p, n_quad_1d - 1 - iQuad);
                surf_intR[2*n_p + i_p] += 0.5 * w1d[iQuad] * flux[2] * basis.PhiSide(sideR, i_p, n_quad_1d - 1 - iQuad);
                surf_intR[3*n_p + i_p] += 0.5 * w1d[iQuad] * flux[3] * basis.PhiSide(sideR, i_p, n_quad_1d - 1 - iQuad);
            }
        }

        // Adding to RHS. Note the negative flux from left element is accounted for
        res_dg.AddBlock(iElmR, surf_intR);
        res_dg.AddBlock(iElmL, surf_intL);
    }
}

void EulerDG::Boundaries(const SolVector &Uc) {
    const vector<Bound> &boundaries = grid->GetBoundaries();
    double u_interior[nVar], u_bndry[nVar];
    double v_interior[nVar], v_bndry[nVar];
    double surf_int[n_p*nVar];
    for (auto &bound : boundaries) {
        const Bound::BCType bType = bound.Type();
        const vector<Edge> &edges = bound.GetEdges();
        for (const auto &edge : edges) {
            const auto iElmR = edge.ElmR();
            const vector<double> &norm = edge.Norm();
            const double len = edge.Length();
            const auto sideR = edge.EdgeNumR();

            const double *c = Uc.GetBlock(iElmR);

            for (auto i_p = 0ul; i_p < n_p; i_p++) {
                surf_int[0*n_p + i_p] = 0.0;
                surf_int[1*n_p + i_p] = 0.0;
                surf_int[2*n_p + i_p] = 0.0;
                surf_int[3*n_p + i_p] = 0.0;
            }

            for (auto iQuad = 0ul; iQuad < n_quad_1d; iQuad++) {
                EvalRightU(c, iQuad, sideR, u_interior);
                ConsToPrim(u_interior, v_interior);
                EvalBoundPrim(bType, norm.data(), v_interior, v_bndry);
                conv_flux->SetStates(v_interior, v_bndry, norm.data(), len);
                conv_flux->ComputeFlux();
                const vector<double> &flux = conv_flux->Flux();

                // See note in surface integral about how jacobian is handled
                for (auto i_p = 0ul; i_p < n_p; i_p++) {
                    surf_int[0*n_p + i_p] += -0.5 * w1d[iQuad] * flux[0] * basis.PhiSide(sideR, i_p, iQuad);
                    surf_int[1*n_p + i_p] += -0.5 * w1d[iQuad] * flux[1] * basis.PhiSide(sideR, i_p, iQuad);
                    surf_int[2*n_p + i_p] += -0.5 * w1d[iQuad] * flux[2] * basis.PhiSide(sideR, i_p, iQuad);
                    surf_int[3*n_p + i_p] += -0.5 * w1d[iQuad] * flux[3] * basis.PhiSide(sideR, i_p, iQuad);
                }
            }

            res_dg.AddBlock(iElmR, surf_int);
        }
    }
}

void EulerDG::EvalLeftRight(const double *cL, const double *cR, double *uL, 
        double *uR, const unsigned long iQuad, const unsigned long elmL,
        const unsigned long elmR, const unsigned short sideL, const unsigned short sideR) {

        uL[0] = 0.0; uR[0] = 0.0;
        uL[1] = 0.0; uR[1] = 0.0;
        uL[2] = 0.0; uR[2] = 0.0;
        uL[3] = 0.0; uR[3] = 0.0;

    for (auto i = 0ul; i < n_p; i++) {
            uL[0] += cL[0*n_p + i] * basis.PhiSide(sideL, i, iQuad);
            uL[1] += cL[1*n_p + i] * basis.PhiSide(sideL, i, iQuad);
            uL[2] += cL[2*n_p + i] * basis.PhiSide(sideL, i, iQuad);
            uL[3] += cL[3*n_p + i] * basis.PhiSide(sideL, i, iQuad);

            uR[0] += cR[0*n_p + i] * basis.PhiSide(sideR, i, n_quad_1d - 1 - iQuad);
            uR[1] += cR[1*n_p + i] * basis.PhiSide(sideR, i, n_quad_1d - 1 - iQuad);
            uR[2] += cR[2*n_p + i] * basis.PhiSide(sideR, i, n_quad_1d - 1 - iQuad);
            uR[3] += cR[3*n_p + i] * basis.PhiSide(sideR, i, n_quad_1d - 1 - iQuad);
    }
}

void EulerDG::InscribedCircles() {
    const vector<Elm> elms = grid->GetElms();
    const vector<Node> nodes = grid->GetNodes();
    for (auto iElm = 0ul; iElm < nEqn; iElm++) {
        const auto &elm_nodes = elms[iElm].GetNodes();
        const Node &n1 = nodes[elm_nodes[0]];
        const Node &n2 = nodes[elm_nodes[1]];
        const Node &n3 = nodes[elm_nodes[2]];
        double a = sqrt(pow(n1.X() - n2.X(), 2.0) + pow(n1.Y() - n2.Y(), 2.0));
        double b = sqrt(pow(n2.X() - n3.X(), 2.0) + pow(n2.Y() - n3.Y(), 2.0));
        double c = sqrt(pow(n1.X() - n3.X(), 2.0) + pow(n1.Y() - n3.Y(), 2.0));
        double k = 0.5 * (a + b + c);
        radii[iElm] = sqrt(k * (k-a) * (k-b) * (k-c)) / k;
    }
}

void EulerDG::ComputeJacobians() {
    // J = [ xr, xs; yr, ys]
    const vector<Node> &nodes = grid->GetNodes();
    const vector<Elm> &elms = grid->GetElms();

    for (auto iElm = 0ul; iElm < nEqn; iElm++) {

        const Elm &elm = elms[iElm];
        const auto &elm_nodes = elm.GetNodes();
        const Node &n1 = nodes[elm_nodes[0]];
        const Node &n2 = nodes[elm_nodes[1]];
        const Node &n3 = nodes[elm_nodes[2]];

        xr[iElm] = n2.X() - n1.X();
        yr[iElm] = n2.Y() - n1.Y();
        xs[iElm] = n3.X() - n1.X();
        ys[iElm] = n3.Y() - n1.Y();
    }
}

void EulerDG::ConsToPrim(const double* u, double* v) {
            const double rho = u[0];
            const double vx = u[1] / rho;
            const double vy = u[2] / rho;
            const double E = u[3];
            const double p = (E - 0.5 * rho * (vx * vx + vy * vy)) * (gamma - 1);

            v[0] = rho; v[1] = vx; v[2] = vy; v[3] = p;
}

void EulerDG::EvalRightU(const double *cR, const unsigned long iQuad, const unsigned long sideR, double *uR) {
    uR[0] = 0.0;
    uR[1] = 0.0;
    uR[2] = 0.0;
    uR[3] = 0.0;

    for (auto i = 0ul; i < n_p; i++) {
        uR[0] += cR[0*n_p + i] * basis.PhiSide(sideR, i, iQuad);
        uR[1] += cR[1*n_p + i] * basis.PhiSide(sideR, i, iQuad);
        uR[2] += cR[2*n_p + i] * basis.PhiSide(sideR, i, iQuad);
        uR[3] += cR[3*n_p + i] * basis.PhiSide(sideR, i, iQuad);
    }
}

void EulerDG::EvalBoundPrim(const Bound::BCType bcType, const double *norm, const double *v_interior, double *v_bndry) {
    const double rho_inf = config->GetFreestreamDensity();
    const vector<double> vel_inf = config->GetFreestreamVelocity();
    const double p_inf = config->GetFreestreamPressure();
    double proj_vel = 0.0;
    switch (bcType) {
        case Bound::BCType::Inlet : break; // do these later if I feel like it
        case Bound::BCType::Outlet : break;
        case Bound::BCType::SupersonicInlet :
            v_bndry[0] = rho_inf;
            v_bndry[1] = vel_inf[0];
            v_bndry[2] = vel_inf[1];
            v_bndry[3] = p_inf;
            break;
        case Bound::BCType::SupersonicOutlet :
            v_bndry[0] = v_interior[0];
            v_bndry[1] = v_interior[1];
            v_bndry[2] = v_interior[2];
            v_bndry[3] = v_interior[3];
            break;
        case Bound::BCType::Farfield :
            v_bndry[0] = rho_inf;
            v_bndry[1] = vel_inf[0];
            v_bndry[2] = vel_inf[1];
            v_bndry[3] = p_inf;
            break;
        case Bound::BCType::InviscidWall : 
            v_bndry[0] = v_interior[0];
            v_bndry[1] = v_interior[1];
            v_bndry[2] = v_interior[2];
            proj_vel = v_interior[1] * norm[0] + v_interior[2] * norm[1];

            v_bndry[1] -= proj_vel * norm[0];
            v_bndry[2] -= proj_vel * norm[1];

            v_bndry[nVar-1] = v_interior[nVar-1];
            break;
        case Bound::BCType::ViscousWall : // do later. How to enforce hard BC with FEM/DG coefficients?
            break;
        case Bound::BCType::Empty : cout << "Empty BC!\n"; exit(1); break; // panic
    }
}

double EulerDG::EvalTimeStep(const double cfl) {
    double U[4], V[4];
    double min_rj_a = 1e10;
    for (auto iElm = 0ul; iElm < nEqn; iElm++) {
        const double *c = U_coeff.GetBlock(iElm);

        for (auto iQuad = 0ul; iQuad < n_quad; iQuad++) {
            U[0] = 0.0;
            U[1] = 0.0;
            U[2] = 0.0;
            U[3] = 0.0;

            for (auto k = 0; k < n_p; k++) {
                U[0] += c[0*n_p + k] * basis.Phi(k, iQuad);
                U[1] += c[1*n_p + k] * basis.Phi(k, iQuad);
                U[2] += c[2*n_p + k] * basis.Phi(k, iQuad);
                U[3] += c[3*n_p + k] * basis.Phi(k, iQuad);
            }
        }

        ConsToPrim(U, V);

        const double rho = V[0];
        const double u = V[1];
        const double v = V[2];
        const double p = V[3];

        const double vel_sqr = sqrt(u*u + v*v);
        const double sound_speed = sqrt(gamma * p / rho);
        const double lambda = vel_sqr + sound_speed;
        const double rj_a = radii[iElm] / lambda;
        if (rj_a < min_rj_a) min_rj_a = rj_a;
    }

    return cfl * min_rj_a * 1. / 7.;
}

void EulerDG::Update(const unsigned long max_iter) {
    double t = 0.0;
    double tf = 10.0;
    SolVector rk1, rk2, rk3, rk4, tmp;
    rk1.Init(nEqn, nVar * n_p);
    rk2.Init(nEqn, nVar * n_p);
    rk3.Init(nEqn, nVar * n_p);
    rk4.Init(nEqn, nVar * n_p);
    tmp.Init(nEqn, nVar * n_p);
    const vector<double> &jac = grid->GetJacobiansDG();

    cout << "\n--------------------------\n"
        << "Writing intial conditions.\n"
        << "--------------------------\n";
    VTK::WriteVTK("out_" + std::to_string(t) + ".vtu", *grid, primVar);
    const double write_time_inc = 0.02;
    double next_write = 0.02;
    unsigned long iter = 0;
    // on stability: https://www.sciencedirect.com/science/article/abs/pii/S0021999119308009
    // n_p = 10, rk4  =>  cfl = 0.025

    while (t < tf) {
        double dt = EvalTimeStep(0.5);
        cout << "iter: " << iter << "\ttime: " << t << "\tdt: " << dt << endl;
        if (dt > 1e-3) {
            cout << "time step big something wrong or we converged\n";
            exit(0);
        }

        if ((t+dt) > next_write) dt = next_write - t;

        /*
        res_dg.SetZeroes();
        VolumeIntegral(U_coeff);
        SurfaceIntegral(U_coeff);
        Boundaries(U_coeff);

        for (auto iElm = 0; iElm < nEqn; iElm++) {
            const double *res_i = res_dg.GetBlock(iElm);
            U_coeff.AddBlock(iElm, res_i, dt / jac[iElm]);
        }
        */

        
        // Step 1 - Evaluate at dt = 0, solve for rk1
        res_dg.SetZeroes();
        VolumeIntegral(U_coeff);
        SurfaceIntegral(U_coeff);
        Boundaries(U_coeff);

        for (auto iElm = 0; iElm < nEqn; iElm++) {
            const double *res_i = res_dg.GetBlock(iElm);
            rk1.SetBlock(iElm, res_i, dt / jac[iElm]);
        }

        // Step 2 - half time step with U_coeff + rk1, solve rk2
        U_coeff.AddVector(rk1, rk2, 0.5); // rk2 = U_coeff + 0.5 * rk1
        res_dg.SetZeroes();
        VolumeIntegral(rk2);
        SurfaceIntegral(rk2);
        Boundaries(rk2);

        for (auto iElm = 0; iElm < nEqn; iElm++) {
            const double *res_i = res_dg.GetBlock(iElm);
            rk2.SetBlock(iElm, res_i, dt / jac[iElm]);
        }

        // Step 3 - half time step with U_coeff + rk2, solve rk3
        U_coeff.AddVector(rk2, rk3, 0.5); // rk3 = U_coeff + 0.5 * rk2
        res_dg.SetZeroes();
        VolumeIntegral(rk3);
        SurfaceIntegral(rk3);
        Boundaries(rk3);

        for (auto iElm = 0; iElm < nEqn; iElm++) {
            const double *res_i = res_dg.GetBlock(iElm);
            rk3.SetBlock(iElm, res_i, dt / jac[iElm]);
        }

        // Step 4 - full time step with U_coeff + rk3, solve rk4
        U_coeff.AddVector(rk3, rk4);
        res_dg.SetZeroes();
        VolumeIntegral(rk4);
        SurfaceIntegral(rk4);
        Boundaries(rk4);

        for (auto iElm = 0; iElm < nEqn; iElm++) {
            const double *res_i = res_dg.GetBlock(iElm);
            rk4.SetBlock(iElm, res_i, dt / jac[iElm]);
        }

        // 4th order Runge-Kutta Update
        U_coeff.AddVector(rk1, 1./6.);
        U_coeff.AddVector(rk2, 1./3.);
        U_coeff.AddVector(rk3, 1./3.);
        U_coeff.AddVector(rk4, 1./6.);


        t += dt;
        iter++;

        if (abs(t - next_write) < 1e-12) {
    cout << "\n--------------------------\n"
        << "Writing solution at time: " << t << endl
        << "--------------------------\n";
            cout << "writing solution at time: " << t << endl;
            next_write += write_time_inc;
            UpdateConVars();
            UpdatePrimitiveVars();
            string filename = "out_" + std::to_string(t) + ".vtu";
            VTK::WriteVTK(filename, *grid, primVar);
        }

    }
}

void EulerDG::Solve() {
    InitEulerDG();
    conVar.SetZeroes();
    UpdateConVars();
    UpdatePrimitiveVars();

    unsigned long max_iter = 100000;
    Update(max_iter);

    UpdateConVars();
    UpdatePrimitiveVars();
}
