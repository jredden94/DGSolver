#pragma once

#include <vector>
#include <mutex>
#include <iostream>

using namespace std;

class Config {
    public:
        enum class ConvFlux { Roe, HLLC };
        enum class SolverType { EulerFVM, EulerDG };
        enum class OutputFile { VTK };
        enum class Preconditioner { Jacobi };
        enum class LinearSolver { BiCGSTAB, GMRES };

        static Config& GetConfig(void);

        ~Config();
        
        void SetDefault(void);

        const string& GetGridFilename(void) const;
        const ConvFlux& GetConvFluxScheme(void) const;
        const SolverType& GetSolverType(void) const;
        const LinearSolver& GetLinearSolverType(void) const;
        const vector<double>& GetFreestreamVelocity(void) const;
        void SetFreestreamVelocity(const vector<double> &vel);
        const double& GetFreestreamDensity(void) const;
        const double& GetFreestreamPressure(void) const;
        const double& GetFreestreamMach(void) const;
        const double& GetAngleOfAttack(void) const;
        void SetAngleOfAttack(const double aoa);
        const double& GetGamma(void) const;
        const double& GetGasConstant(void) const;
        const double& GetCFL(void) const;
        const double& GetMinCFL(void) const;
        const double& GetMaxCFL(void) const;
        const double& GetCFL_FactorUp(void) const;
        const double& GetCFL_FactorDown(void) const;
        const bool& AdaptiveCFL(void) const;
        const bool& IsSteady(void) const;
        const unsigned short& GetNumDims(void) const;
        const unsigned short& GetNumVars(void) const;
        const unsigned long& GetNumEqn(void) const;
        const double& GetConvergenceTolerance(void) const;
        const unsigned long& GetLinearSolverMaxIterations(void) const;
        void SetLinearSolverMaxIterations(const unsigned long max_lin_iter);
        const unsigned long& GetMaxIterations(void) const;
        void SetMaxIterations(const unsigned long max_iter);
        const unsigned long& GetMinIterations(void) const;
        const double& GetLinearSolverTolerance(void) const;
        //const Limiter::Type& GetLimiterType(void) const;
        const double& GetLimiterCoefficient(void) const;
        const bool IsImplicit(void) const;
        const bool MUSCL(void) const;
        const bool IsLimited(void) const;
        const Preconditioner GetPreconditioner(void) const;
        const double& GetViscosity(void) const;
        const bool& IsViscous(void) const;
        const double& GetRoeDissipationCoefficient(void) const;
        const double& GetReferenceArea(void) const;

        void SetNumEqn(const unsigned long&);

        Config(const Config&) = delete;
        Config& operator+(const Config&) = delete;

        // DG
        static constexpr unsigned short n_order = 3;
        static constexpr unsigned short n_polynomial = 10;
        static constexpr unsigned short n_quad = 12;
        static constexpr unsigned short n_quad_1d = 4;
        static constexpr unsigned short nDim = 2;
        static constexpr unsigned short nVar = 4;

    private:
        Config();
        static void Init();
        static Config *instance;
        static once_flag onceFlag;

        string grid_file;

        ConvFlux convFlux;
        SolverType solver;
        LinearSolver lin_solver;
        Preconditioner precond;

        //Limiter::Type limiterType;
        double venkat_lim_coeff;

        vector<double> vel_inf;
        double rho_inf, p_inf, M_inf;
        double aoa; // angle of attack
        double gamma, gamma_m1, gasConstant;
        bool steady, implicit, muscl, isLimited;
        unsigned long nEqn;

        double viscosity;
        bool isViscous;

        bool adaptCFL;
        double cfl_start, cfl_min, cfl_max, cfl_factor_up, cfl_factor_down;
        double convergence_tol, lin_solver_tol;
        unsigned long max_lin_solver_iter, max_iter, min_iter;

        double roe_diss_coeff;

        double ref_area;

};
