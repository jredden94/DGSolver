#pragma once

#include <memory>
#include "../Common/Config.hpp"
#include "../Grid/Grid.hpp"
#include "../Grid/Edge.hpp"
#include "../Grid/Node.hpp"
#include "../Numerics/Convection/Convection.hpp"
#include "../LinAlg/SolVector.hpp"
#include "../LinAlg/SolMatrix.hpp"
#include "../Numerics/Viscous/Viscous.hpp"
#include "../Basis/Basis.hpp"
#include "../Basis/Quad.hpp"
#include "../IO/VTK.hpp"

using namespace std;

class Solver {
    public:
        Solver(const Solver&) = delete;
        virtual ~Solver();
        static std::unique_ptr<Solver> CreateSolver(Grid*);
        virtual void Solve() = 0;

        Solver& operator=(const Solver&) = delete;
        void Init(void);

        const SolVector& GetPrimitives(void) const;
        const SolVector& GetConservatives(void) const;

    protected:
        Solver(Config::SolverType);
        double *v_i, *v_j;
        void UpdatePrimitiveVars(void);
        void PrintResiduals(const unsigned long current_iter) const;
        void PrintResiduals(const unsigned long current_iter, const double t) const;
        void PrintResiduals(const unsigned long current_iter, const double cD, const double cL) const;
        void AdaptCFL(void);

        Config *config;
        mutable bool converged = false;

        bool implicit, muscl;
        static constexpr unsigned short nVar = Config::nVar;;
        static constexpr unsigned short nDim = Config::nDim;
        unsigned long nEqn;

        SolVector delU, residual, conVar, primVar;
        SolMatrix jacMat;

        double cfl, gamma;
        vector<double> waveSpeed, dt;
        Grid *grid;
        unique_ptr<Convection> conv_flux;
        bool viscous;
        Viscous visc_flux;

        double tol;
        mutable vector<double> res_norm;
        vector<double> nonlinear_res;
        unsigned long nonlinear_res_counter = 0, nonlinear_res_max_count = 25;

};

class Euler : public Solver {
    public:
        Euler();
        ~Euler() override;

        void Solve() override;
        void BoundaryFluxes(const Bound &boundary);
        void InviscidWall(const Bound &wall);
        void ViscousWall(const Bound &wall);
        void FarField(const Bound &far);
        void Inlet(const Bound &inlet);
        void Outlet(const Bound &outlet);
        void SupersonicInlet(const Bound &inlet);
        void SupersonicOutlet(const Bound &outlet);
};

class EulerDG : public Solver {
    public:
        EulerDG();
        ~EulerDG() override;

        void Solve() override;

        void InitEulerDG(void);
        void InitCoefficients();

    protected:
        Basis basis;
        Quad quad;

        // Quadrature points and weights
        vector<double> r1, r2, w, s_r, w1d;

        // components of Jacobians
        vector<double> xr, yr, xs, ys;

        // diameter of inscribed circles
        vector<double> radii;

        static constexpr unsigned short n_p = Config::n_polynomial;
        static constexpr unsigned short n_quad = Config::n_quad;
        static constexpr unsigned short n_quad_1d = Config::n_quad_1d;
        SolVector res_dg, U_coeff; // residual, solution coefficients

        void UpdateConVars(void);
        void EvalU0(double *U, int i, unsigned long iBlock);

        // Evaluate Jacobian components xr, xs, yr, ys J = [ xr, xs; yr, yr]
        void ComputeJacobians(void);

        // radii of inscribed circles
        void InscribedCircles(void);

        double EvalTimeStep(const double cfl);

        // \int_{\Omega} F grad(psi)
        void VolumeIntegral(const SolVector &Uc);

        // \int_{\partial \Omega} F^* psi
        void SurfaceIntegral(const SolVector &Uc);

        // Surface integral term at boundary edges
        void Boundaries(const SolVector &Uc);

        void EvalRightU(const double *coeff, const unsigned long iQuad, const unsigned long sideR, double *u_quad);

        void EvalBoundPrim(const Bound::BCType bcType, const double *norm, const double *v_interior, double *v_bndry);

        void EvalLeftRight(const double *cL, const double *cR, double *uL, 
                double *uR, const unsigned long iQuad, const unsigned long elmL,
                const unsigned long elmR, const unsigned short sideL, const unsigned short sideR);

        void ConsToPrim(const double* u, double* v);

        void Update(const unsigned long iter);
};
