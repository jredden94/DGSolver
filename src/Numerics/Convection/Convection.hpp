#pragma once

#include <memory>
#include <vector>
#include <cmath>
#include "../../Common/Config.hpp"

using namespace std;

class Convection {
    public:
        Convection();
        Convection(Config::ConvFlux);
        virtual ~Convection();

        virtual void ComputeFlux(void) = 0;

        void SetStates(const vector<double>& v_i, const vector<double>& v_j, 
                const vector<double>& norm, const double &area_mag);

        void SetStates(const double* v_i, const double* v_j, 
                const double* areaVector, const double &area_mag);

        const vector<double>& Flux(void);
        const double* JacI(void);
        const double* JacJ(void);
        const double& MaxWaveSpeed(void);
        const double MaxLambda(void); // same as max wave speed, just not scaled by area for DG solver

        static unique_ptr<Convection> CreateConvFlux(Config *config);

        void PrintIJacobianInfo(void);
        void PrintJJacobianInfo(void);

        void ComputeJacobian(const double *vel, const double &vel_sqr, 
                const double &proj_vel, const double *norm, 
                const double &energy, double *jac, const double scale = 0.5);

    protected:

        void ComputeFlux(const double &rho, const double *vel, const double &enthalpy,
        const double *norm, const double *flux);
        void ComputePTensor(const double &rho, const double *vel,
                const double &soundSpeed, const double *norm, double *p_tensor);
        void ComputeInversePTensor(const double &rho, const double *vel, 
                const double &soundspeed, const double *norm, double *inv_p);

        Config *config;
        static constexpr unsigned short nVar=4, nDim=2;
        bool implicit;

        vector<double> norm, flux, eigen; 
        double max_lambda, maxWaveSpeed, entropy_fix_coeff; 
        double area;

        double gamma, gamma_minus_one;

        unsigned short iDim, jDim, kDim, iVar, jVar, kVar;

        // Left state values
        double *v_i, *u_i;
        vector<double> vel_i, flux_i; 
        double rho_i, p_i, enthalpy_i, energy_i;
        double vel_sqr_i, proj_vel_i, soundSpeed_i;

        // Right state values
        double *v_j, *u_j;
        vector<double> vel_j, flux_j; 
        double rho_j, p_j, enthalpy_j, energy_j;
        double vel_sqr_j, proj_vel_j, soundSpeed_j;

        // Roe Averages
        vector<double> vel_roe;
        double rhoR, rho, p, enthalpy, energy;
        double proj_vel, vel_sqr_roe, soundSpeed, soundSpeed2;

        // Jacobians
        double *jac_i, *jac_j;
};

class Roe : public Convection {
    public:
        Roe();
        ~Roe() override;

        void ComputeFlux() override;

    private:
        double *p_tensor, *p_inv, *del_u;
        double p_lam_pinv;
        double diss_coeff;
};

class HLLC : public Convection {
    public:
        HLLC();
        ~HLLC() override;

        void ComputeFlux() override;

    private:
        double sL, sR, sM;
        double rho_m, pStar, rhoSL, rhoSR;

        double *inter_state;
        double eStar, omega, omegaSM;
        double *dp_du_i, *dsM_du, *drhoStar_du, *dpStar_du, *deStar_du;
};
