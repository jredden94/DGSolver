#pragma once

#include "../../Common/Config.hpp"

class Viscous {
    public:
        Viscous(void);
        ~Viscous(void);

        void SetStates(const double *v_i, const double *v_j, const double *norm,
                const double area, const double edge_len, const double *vel_grad_i, const double *vel_grad_j);
        void ComputeResidual(void);

        const double* Flux(void) const;
        const double* JacI(void) const;
        const double* JacJ(void) const;
        const double* StressTensor(void) const;
        void ComputeStressTensor(const double *vel_grad, const double viscosity, double *tau );

    private:
        void ComputeViscousFlux(void);
        void ComputeTauJacobian(void);
        void ComputeJacobians(void);

        unsigned short nVar, nDim;

        double viscosity;
        double *tau, *visc_flux;
        double *v_i, *v_j, *v_mean, *norm, *coord_i, *coord_j;
        double area, edge_len;
        double *vel_grad_i, *vel_grad_j, *vel_grad_mean;


        bool implicit;
        double *tau_jac, *jac_i, *jac_j;
};
