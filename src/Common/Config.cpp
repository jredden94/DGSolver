#include "Config.hpp"
#include <cmath>

Config* Config::instance = nullptr;
once_flag Config::onceFlag;

Config::Config() { }
Config::~Config() { delete instance; }

Config& Config::GetConfig() { 
    call_once(onceFlag, &Config::Init);
    return *instance;
}

void Config::Init() {
    instance = new Config;
    instance->SetDefault();
}

void Config::SetDefault() {
    grid_file = "naca.grid";
    solver = Config::SolverType::EulerFVM;
    lin_solver = Config::LinearSolver::GMRES;
    convFlux = Config::ConvFlux::Roe;
    rho_inf = 1.0; //1.00042; //5.2582e-04; //1.13831;
    vel_inf.resize(3);
    aoa = 5.0;
    vel_inf[0] = 0.3 * cosf(aoa * M_PI / 180.); //0.946 * cosf(aoa * M_PI / 180.0); //170.1 * cosf(aoa * M_PI / 180.0);
    vel_inf[1] = 0.3 * sinf(aoa * M_PI / 180.); //0.946 * sinf(aoa * M_PI / 180.0); //170.1 * sinf(aoa * M_PI / 180.0);
    vel_inf[2] = 0.0;
    gamma = 1.4;
    p_inf = rho_inf / gamma; //1.0 / gamma; //100000; //1.0 / gamma; // nondimensionalization such that u_inf = M
    gamma_m1 = gamma - 1.0;

    roe_diss_coeff = 1.0;

    cfl_min = 0.9;
    cfl_start = 0.9;
    cfl_max = 10.0;
    cfl_factor_up = 1.1;
    cfl_factor_down = 0.9;
    adaptCFL = false;

    convergence_tol = 1e-15;
    max_iter = 10000;
    min_iter = 10000;
    max_lin_solver_iter = 1;
    lin_solver_tol = 1e-12;
    implicit = false;
    //limiterType = Limiter::Type::Venkatakrishnan;
    muscl = false;
    venkat_lim_coeff = 0.0001;
    isLimited = false;
    precond = Config::Preconditioner::Jacobi;

    isViscous = false;
    viscosity = 1. / 5000; // 1.78941e-05; //1.83751e-05;

    ref_area = 1.0;
}

const string& Config::GetGridFilename(void) const { return grid_file; }
const Config::ConvFlux& Config::GetConvFluxScheme(void) const { return convFlux; }
const Config::SolverType& Config::GetSolverType(void) const { return solver; }
const Config::LinearSolver& Config::GetLinearSolverType(void) const { return lin_solver; }
const vector<double>& Config::GetFreestreamVelocity(void) const { return vel_inf; }
void Config::SetFreestreamVelocity(const vector<double> &vel) { vel_inf = vel; }
const double& Config::GetFreestreamDensity(void) const { return rho_inf; }
const double& Config::GetFreestreamPressure(void) const { return p_inf; }
const double& Config::GetFreestreamMach(void) const { return M_inf; }
const double& Config::GetAngleOfAttack(void) const { return aoa; }
void Config::SetAngleOfAttack(const double aoa) { this->aoa = aoa; }
const double& Config::GetGamma(void) const { return gamma; }
const double& Config::GetGasConstant(void) const { return gasConstant; }
const double& Config::GetCFL(void) const { return cfl_start; }
const double& Config::GetMinCFL(void) const { return cfl_min; }
const double& Config::GetMaxCFL(void) const { return cfl_max; }
const double& Config::GetCFL_FactorUp(void) const { return cfl_factor_up; }
const double& Config::GetCFL_FactorDown(void) const { return cfl_factor_down; }
const bool& Config::AdaptiveCFL(void) const { return adaptCFL; }
const bool& Config::IsSteady(void) const { return steady; }
const unsigned short& Config::GetNumDims(void) const { return nDim; }
const unsigned short& Config::GetNumVars(void) const { return nVar; }
const unsigned long& Config::GetNumEqn(void) const { return nEqn; }
const double& Config::GetConvergenceTolerance(void) const { return convergence_tol; }
const unsigned long& Config::GetLinearSolverMaxIterations(void) const { return max_lin_solver_iter; }
void Config::SetLinearSolverMaxIterations(const unsigned long max_lin_iter) { this->max_lin_solver_iter = max_lin_iter; }
const unsigned long& Config::GetMaxIterations(void) const { return max_iter; }
void Config::SetMaxIterations(const unsigned long max_iter) { this->max_iter = max_iter; }
const unsigned long& Config::GetMinIterations(void) const { return min_iter; }
const double& Config::GetLinearSolverTolerance(void) const { return lin_solver_tol; };
const bool Config::IsImplicit(void) const { return implicit; }
const bool Config::MUSCL(void) const { return muscl; }
const bool Config::IsLimited(void) const { return isLimited; }
//const Limiter::Type& Config::GetLimiterType(void) const { return limiterType; }
const double& Config::GetLimiterCoefficient(void) const { return venkat_lim_coeff; }
const Config::Preconditioner Config::GetPreconditioner(void) const { return precond; }
const double& Config::GetViscosity(void) const { return viscosity; }
const bool& Config::IsViscous(void) const { return isViscous; }
const double& Config::GetRoeDissipationCoefficient(void) const { return roe_diss_coeff; }
const double& Config::GetReferenceArea(void) const { return ref_area; }

void Config::SetNumEqn(const unsigned long& nEqn) { this->nEqn = nEqn; }
