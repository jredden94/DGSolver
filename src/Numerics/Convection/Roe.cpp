#include "Convection.hpp"

Roe::Roe() : Convection(Config::ConvFlux::Roe) { 
    p_tensor = new double[nVar * nVar];
    p_inv = new double[nVar * nVar];
    del_u = new double[nVar];
    entropy_fix_coeff = 0.001;
    diss_coeff = config->GetRoeDissipationCoefficient();
}
Roe::~Roe() { 
    delete[] p_tensor; 
    delete[] p_inv;
    delete[] del_u;
}

void Roe::ComputeFlux() {

    /* Velocities */
    vel_sqr_i = 0;
    proj_vel_i = 0;
    vel_sqr_j = 0;
    proj_vel_j = 0;
    vel_sqr_roe = 0;
    proj_vel = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
        vel_i[iDim] = v_i[iDim+1]; // / rho_i;
        vel_sqr_i += vel_i[iDim] * vel_i[iDim];
        proj_vel_i += vel_i[iDim] * norm[iDim];

        vel_j[iDim] = v_j[iDim+1];
        vel_sqr_j += vel_j[iDim] * vel_j[iDim];
        proj_vel_j += vel_j[iDim] * norm[iDim];

        vel_roe[iDim] = (vel_i[iDim] + rhoR * vel_j[iDim]) / (1 + rhoR);
        vel_sqr_roe += vel_roe[iDim] * vel_roe[iDim];
        proj_vel += vel_roe[iDim] * norm[iDim];
    }

    /* Left State */
    rho_i = v_i[0];
    p_i = v_i[nVar-1];
    soundSpeed_i = sqrt(gamma * p_i / rho_i);
    enthalpy_i = soundSpeed_i * soundSpeed_i 
        / gamma_minus_one + 0.5 * vel_sqr_i;
    energy_i = enthalpy_i - p_i / rho_i;

    /* Right State */
    rho_j = v_j[0];
    p_j = v_j[nVar-1];
    soundSpeed_j = sqrt(gamma * p_j / rho_j);
    enthalpy_j = soundSpeed_j * soundSpeed_j 
        / gamma_minus_one + 0.5 * vel_sqr_j;
    energy_j = enthalpy_j - p_j / rho_j;

    // Roe averages
    rhoR = sqrt(fabs(rho_j / rho_i));
    rho = rhoR * rho_i;
    enthalpy = (enthalpy_i + rhoR * enthalpy_j) / (1 + rhoR);
    soundSpeed2 = gamma_minus_one * (enthalpy - 0.5 * vel_sqr_roe);

    if (soundSpeed2 <= 0) {
        for (iVar = 0; iVar < nVar; iVar++) flux[iVar] = 0.0;

        /*
        if (implicit) {
            for (iVar = 0; iVar < nVar * nVar; iVar++) {
                jac_i[iVar] = 0.0;
                jac_j[iVar] = 0.0;
            }
        }
        */
        return;
    }

    soundSpeed = sqrt(soundSpeed2);

    /* Projected Fluxes, eigenvalues, reconstruct conservative variables */
    flux_i[0] = rho_i * proj_vel_i;
    flux_j[0] = rho_j * proj_vel_j;

    u_i[0] = rho_i;
    u_j[0] = rho_j;

    for (iDim = 0; iDim < nDim; iDim++) {
        flux_i[iDim+1] = rho_i * vel_i[iDim] * proj_vel_i + p_i * norm[iDim];
        flux_j[iDim+1] = rho_j * vel_j[iDim] * proj_vel_j + p_j * norm[iDim];

        eigen[iDim] = proj_vel;

        u_i[iDim+1] = rho_i * vel_i[iDim];
        u_j[iDim+1] = rho_j * vel_j[iDim];
    }
    flux_i[nVar-1] = rho_i * enthalpy_i * proj_vel_i;
    flux_j[nVar-1] = rho_j * enthalpy_j * proj_vel_j;

    eigen[nVar-2] = proj_vel + soundSpeed;
    eigen[nVar-1] = proj_vel - soundSpeed;

    u_i[nVar-1] = p_i / gamma_minus_one + 0.5 * rho_i * vel_sqr_i;
    u_j[nVar-1] = p_j / gamma_minus_one + 0.5 * rho_j * vel_sqr_j;

    /* Entropy Fix, delta U, solve Roe flux */
    maxWaveSpeed = fabs(proj_vel) + soundSpeed;
    for (iVar = 0; iVar < nVar; iVar++)  {
        eigen[iVar] = max(fabs(eigen[iVar]), entropy_fix_coeff * maxWaveSpeed);
        del_u[iVar] = u_j[iVar] - u_i[iVar];
        flux[iVar] = 0.5 * area * (flux_i[iVar] + flux_j[iVar]);
    }
    max_lambda = maxWaveSpeed;
    maxWaveSpeed *= area;

    ComputePTensor(rho, vel_roe.data(), soundSpeed, norm.data(), p_tensor);
    ComputeInversePTensor(rho, vel_roe.data(), soundSpeed, norm.data(), p_inv);

    /* Compute Jacobians if implicit */
    /*
    if (implicit) {
        ComputeJacobian(vel_i.data(), vel_sqr_i, proj_vel_i, norm.data(), energy_i, jac_i, 0.5 * area);
        ComputeJacobian(vel_j.data(), vel_sqr_j, proj_vel_j, norm.data(), energy_j, jac_j, 0.5 * area);
        for (iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
                p_lam_pinv = 0.0;
                for (kVar = 0; kVar < nVar; kVar++) 
                    p_lam_pinv += p_tensor[iVar * nVar + kVar] * eigen[kVar] * p_inv[kVar * nVar + jVar] * diss_coeff;

                flux[iVar] -= 0.5 * p_lam_pinv * del_u[jVar] * area;
                if (implicit) {
                    jac_i[iVar * nVar + jVar] += 0.5 * p_lam_pinv * area;
                    jac_j[iVar * nVar + jVar] -= 0.5 * p_lam_pinv * area;
                }
            }
        }
    }
    */
    //else { /* Just compute flux */
        for (iVar = 0; iVar < nVar; iVar++) {
            for (jVar = 0; jVar < nVar; jVar++) {
                p_lam_pinv = 0.0;
                for (kVar = 0; kVar < nVar; kVar++) 
                    p_lam_pinv += p_tensor[iVar * nVar + kVar] * eigen[kVar] * p_inv[kVar * nVar + jVar];

                flux[iVar] -= 0.5 * p_lam_pinv * del_u[jVar] * area;
     //       }
        }
    }
}
