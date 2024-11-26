#include "Solver.hpp"
#include <memory>

Solver::Solver(Config::SolverType) { }
Solver::~Solver() { 
    delete[] v_i;
    delete[] v_j;
}

unique_ptr<Solver> Solver::CreateSolver(Grid *grid) {
    auto &config = Config::GetConfig();
    Config::SolverType solver_type = config.GetSolverType();
    unique_ptr<Solver> solver;
    switch (solver_type) {
        case Config::SolverType::EulerFVM:
            solver = make_unique<Euler>();
            solver->grid = grid;
            solver->Init();
            return solver;

        case Config::SolverType::EulerDG:
            solver = make_unique<EulerDG>();
            solver->grid = grid;
            solver->Init();
            return solver;

        default:
            solver = make_unique<Euler>();
            solver->grid = grid;
            solver->Init();
            return solver;
    }
}

void Solver::Init() {
    auto &config = Config::GetConfig();
    this->config = &config;
    config.SetNumEqn(grid->GetNumCells());
    nEqn = config.GetNumEqn();
    residual.Init(nEqn, nVar);
    conVar.Init(nEqn, nVar);
    primVar.Init(nEqn, nVar);

    v_i = new double[nVar];
    v_j = new double[nVar];
    tol = config.GetConvergenceTolerance();
    cfl = config.GetCFL();

    dt.resize(nEqn, 0.0);
    waveSpeed.resize(nEqn, 0.0);
    conv_flux = Convection::CreateConvFlux(&config);
    const vector<double> &vel_inf = config.GetFreestreamVelocity();
    const double rho_inf = config.GetFreestreamDensity();
    const double p_inf = config.GetFreestreamPressure();
    gamma = config.GetGamma();

    vector<double> u0(nVar, 0.0);
    u0[0] = rho_inf;
    double vel_sqr = 0;
    for (unsigned short i = 0; i < nDim; i++) {
        u0[i+1] = vel_inf[i] * rho_inf;
        vel_sqr += vel_inf[i] * vel_inf[i];
    }
    u0[nVar-1] = p_inf / (gamma-1) + 0.5 * rho_inf * vel_sqr;
        
    conVar.InitValues(u0);
    UpdatePrimitiveVars();

    implicit = config.IsImplicit();
    muscl = config.MUSCL();
    if (implicit) jacMat.Init(grid);
    viscous = config.IsViscous();
}

const SolVector& Solver::GetPrimitives(void) const { return primVar; }
const SolVector& Solver::GetConservatives(void) const { return conVar; }

void Solver::UpdatePrimitiveVars(void) {
    unsigned long nBlock = primVar.GetBlockCount();
    for (auto i = 0ul; i < nBlock; i++) {
        primVar[i][0] = conVar[i][0];
        double vel_sqr = 0;
        for (unsigned short j = 0; j < nDim; j++) {
            primVar[i][j+1] = conVar[i][j+1] / conVar[i][0];
            vel_sqr += primVar[i][j+1] * primVar[i][j+1];
        }
        primVar[i][nVar-1] = (conVar[i][nVar-1] - 0.5 * conVar[i][0] * vel_sqr) * (gamma - 1);

    }
}

void Solver::PrintResiduals(const unsigned long current_iter) const {
    converged = true;

    res_norm = residual.ResNorm();
    cout << current_iter << "\t";
    for (size_t i = 0; i < nVar; i++) {
        cout << res_norm[i] << "\t";
        if (res_norm[i] >= tol) converged = false;
    }
    cout << endl;
}

void Solver::PrintResiduals(const unsigned long current_iter, const double t) const {
    converged = true;

    res_norm = residual.ResNorm();
    cout << current_iter << "\t";
    for (size_t i = 0; i < nVar; i++) {
        cout << res_norm[i] << "\t";
        if (res_norm[i] >= tol) converged = false;
    }
    cout << t << endl;
}

void Solver::PrintResiduals(const unsigned long current_iter, const double cD, const double cL) const {
    converged = true;

    res_norm = residual.ResNorm();
    cout << current_iter << "\t";
    for (size_t i = 0; i < nVar; i++) {
        cout << res_norm[i] << "\t";
        if (res_norm[i] >= tol) converged = false;
    }
    cout << cL << "\t" << cD << endl;
}

void Solver::AdaptCFL() {
    if (nonlinear_res.empty()) nonlinear_res.resize(nonlinear_res_max_count);

    double new_scaled_res_sum = 0.0;
    for (unsigned short i = 0; i < nVar; i++)
        new_scaled_res_sum += log10(res_norm[i]);

    double old_scaled_res_sum = nonlinear_res_counter == 0 ? new_scaled_res_sum : nonlinear_res[nonlinear_res_counter - 1];
    nonlinear_res[nonlinear_res_counter++] = new_scaled_res_sum - old_scaled_res_sum;

    bool resetCFL, reduceCFL, increaseCFL;

    if (nonlinear_res_counter >= nonlinear_res_max_count) {
        nonlinear_res_counter = 0;
        unsigned long sign_changes = 0;
        double total_change = 0.0;
        double prev = nonlinear_res.front();
        for (const auto& val : nonlinear_res) {
            total_change += val;
            sign_changes += (prev > 0) ^ (val > 0);
            prev = val;
        }

        reduceCFL = (sign_changes > nonlinear_res_max_count / 4) || (total_change > -0.2);

        if (total_change > 1.0) {
            nonlinear_res_counter = 0;
            cfl = config->GetMinCFL();
        }
        else if (reduceCFL) cfl = max(config->GetMinCFL(), cfl * config->GetCFL_FactorDown());
        else cfl = min(cfl * config->GetCFL_FactorUp(), config->GetMaxCFL());

        if (cfl == config->GetMaxCFL()) return;
        cout << "\n--------------------------------------------\n";
        cout << "CFL updated to " << cfl << endl;
        cout << "--------------------------------------------\n";
    }
}
