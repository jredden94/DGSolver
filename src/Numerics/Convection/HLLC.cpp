#include "Convection.hpp"

HLLC::HLLC() : Convection(Config::ConvFlux::HLLC) { 
    inter_state = new double[nVar];
    if (implicit) {
        dp_du_i = new double[nVar];
        dsM_du = new double[nVar];
        drhoStar_du = new double[nVar];
        dpStar_du = new double[nVar];
        deStar_du = new double[nVar];
    }
}
HLLC::~HLLC() { 
    delete[] inter_state;
    if (implicit) {
        delete[] dp_du_i;
        delete[] dsM_du;
        delete[] drhoStar_du;
        delete[] dpStar_du;
        delete[] deStar_du;
    }
}

void HLLC::ComputeFlux() { 
    rho_i = v_i[0];
    rho_j = v_j[0];
    rhoR = sqrt(rho_j / rho_i);

    vel_sqr_i = proj_vel_i = 0; 
    vel_sqr_j = proj_vel_j = 0;
    vel_sqr_roe = proj_vel = 0;
    for (iDim = 0; iDim < nDim; iDim++) {
        vel_i[iDim] = v_i[iDim+1];
        vel_sqr_i += vel_i[iDim] * vel_i[iDim];
        proj_vel_i += vel_i[iDim] * norm[iDim];
        vel_j[iDim] = v_j[iDim+1];
        vel_sqr_j += vel_j[iDim] * vel_j[iDim];
        proj_vel_j += vel_j[iDim] * norm[iDim];
        vel_roe[iDim] = (vel_i[iDim] + rhoR * vel_j[iDim]) / (1 + rhoR);
        vel_sqr_roe += vel_roe[iDim] * vel_roe[iDim];
        proj_vel += vel_roe[iDim] * norm[iDim];
    }

    p_i = v_i[nVar-1];
    p_j = v_j[nVar-1];
    soundSpeed_i = sqrt(gamma * p_i / rho_i);
    soundSpeed_j = sqrt(gamma * p_j / rho_j);
    enthalpy_i = soundSpeed_i * soundSpeed_i / gamma_minus_one + 0.5 * vel_sqr_i;
    enthalpy_j = soundSpeed_j * soundSpeed_j / gamma_minus_one + 0.5 * vel_sqr_j;
    energy_i = enthalpy_i - p_i / rho_i;
    energy_j = enthalpy_j - p_j / rho_j;

    enthalpy = (enthalpy_i + rhoR * enthalpy_j) / (1 + rhoR);
    soundSpeed = sqrt(gamma_minus_one * (enthalpy - 0.5 * vel_sqr_roe));

    maxWaveSpeed = abs(proj_vel) + soundSpeed;
    max_lambda = maxWaveSpeed;
    maxWaveSpeed *= area;

    sL = min(proj_vel - soundSpeed, proj_vel_i - soundSpeed_i);
    sR = max(proj_vel + soundSpeed, proj_vel_j + soundSpeed_j);
    rho_m = rho_j * (sR - proj_vel_j) - rho_i * (sL - proj_vel_i);
    sM = ( p_i - p_j - rho_i * proj_vel_i * (sL - proj_vel_i) + rho_j * proj_vel_j * (sR - proj_vel_j)) / rho_m;
    pStar = rho_j * (proj_vel_j - sR) * (proj_vel_j - sM) + p_j;

    /* Flux */
    if (sM > 0.0) {
        if (sL > 0.0) {
            flux[0] = rho_i * proj_vel_i;
            for (iDim = 0; iDim < nDim; iDim++)
                flux[iDim+1] = rho_i * vel_i[iDim] * proj_vel_i + p_i * norm[iDim];
            flux[nVar-1] = enthalpy_i * rho_i * proj_vel_i;
        }
        else {
            rhoSL = (sL - proj_vel_i) / (sL - sM);

            inter_state[0] = rhoSL * rho_i;
            flux[0] = sM * inter_state[0];
            for (iDim = 0; iDim < nDim; iDim++) {
                inter_state[iDim+1] = rhoSL * (rho_i * vel_i[iDim] + (pStar - p_i) / (sL - proj_vel_i) * norm[iDim]);
                flux[iDim+1] = sM * inter_state[iDim+1] + pStar * norm[iDim];
            }
            inter_state[nVar-1] = rhoSL * (rho_i * energy_i - (p_i * proj_vel_i - pStar * sM) / (sL - proj_vel_i));
            flux[nVar-1] = sM * (inter_state[nVar-1] + pStar);
        }
    }
    else {
        if (sR < 0.0) {
            flux[0] = rho_j * proj_vel_j;
            for (iDim = 0; iDim < nDim; iDim++)
                flux[iDim+1] = rho_j * vel_j[iDim] * proj_vel_j + p_j * norm[iDim];
            flux[nVar-1] = enthalpy_j * rho_j * proj_vel_j;
        }
        else {
            rhoSR = (sR - proj_vel_j) / (sR - sM);

            inter_state[0] = rhoSR * rho_j;
            flux[0] = sM * inter_state[0];
            for (iDim = 0; iDim < nDim; iDim++) {
                inter_state[iDim+1] = rhoSR * (rho_j * vel_j[iDim] + (pStar - p_j) / (sR - proj_vel_j) * norm[iDim]);
                flux[iDim+1] = sM * inter_state[iDim+1] + pStar * norm[iDim];
            }
            inter_state[nVar-1] = rhoSR * (rho_j * energy_j - (p_j * proj_vel_j - pStar * sM) / (sR - proj_vel_j));
            flux[nVar-1] = sM * (inter_state[nVar-1] + pStar);
        }
    }

    for (iVar = 0; iVar < nVar; iVar++) flux[iVar] *= area;

    /* Jacobians - Less readable, tighter version of SU2's implementation */
#define JAC_I(I, J) jac_i[I * nVar + J]
#define JAC_J(I, J) jac_j[I * nVar + J]
    if (implicit) {
        if (sM > 0.0) {
            if (sL > 0.0) {
                for (iVar = 0; iVar < nVar * nVar; iVar++) jac_j[iVar] = 0;
                ComputeJacobian(vel_i.data(), vel_sqr_i, proj_vel_i, norm.data(), energy_i, jac_i);
            }
            else {
                /*--- Compute Jacobian based on Left Star State ---*/

                eStar = inter_state[nVar-1];
                omega = 1/(sL-sM);
                omegaSM = omega * sM;


                /*--------- Left Jacobian ---------*/

                /*--- Computing derivatives ---*/
                dsM_du[0] = ( - proj_vel_i * proj_vel_i + sM * sL + dp_du_i[0] ) / rho_m;
                dp_du_i[0] = 0.5 * gamma_minus_one * vel_sqr_i;
                drhoStar_du[0] = omega * ( sL + inter_state[0] * dsM_du[0] );
                for (iDim = 0; iDim < nDim; iDim++) {
                    dp_du_i[iDim+1] = - gamma_minus_one * vel_i[iDim];
                    dsM_du[iDim+1] = ( norm[iDim] * ( 2 * proj_vel_i - sL - sM ) + dp_du_i[iDim+1] ) / rho_m;
                    drhoStar_du[iDim+1] = omega * ( - norm[iDim] + inter_state[0] * dsM_du[iDim+1] );
                }
                dp_du_i[nVar-1] = gamma_minus_one;
                dsM_du[nVar-1] = dp_du_i[nVar-1] / rho_m;
                drhoStar_du[nVar-1] = omega * inter_state[0] * dsM_du[nVar-1];

                for (iVar = 0; iVar < nVar; iVar++) {
                    dpStar_du[iVar] = rho_i * (sR - proj_vel_j) * dsM_du[iVar];
                    deStar_du[iVar] = omega * ( sM * dpStar_du[iVar] + ( eStar + pStar ) * dsM_du[iVar] );
                    jac_i[iVar] = sM * drhoStar_du[iVar] + inter_state[0] * dsM_du[iVar];
                }

                deStar_du[0] += omega * proj_vel_i * ( enthalpy_i - dp_du_i[0] );
                for (iDim = 0; iDim < nDim; iDim++)
                    deStar_du[iDim+1] += omega * ( - norm[iDim] * enthalpy_i - proj_vel_i * dp_du_i[iDim+1] );
                deStar_du[nVar-1] += omega * ( sL - proj_vel_i - proj_vel_i * dp_du_i[nVar-1] );

                /*--- jac_i ---*/
                for (iDim = 0; iDim < nDim; iDim++) {
                    JAC_I(iDim+1, 0) += omegaSM * vel_i[iDim] * proj_vel_i;
                    for (iVar = 0; iVar < nVar; iVar++) {
                        JAC_I(iDim+1, iVar) = ( omegaSM + 1 ) * ( norm[iDim] * dpStar_du[iVar] + inter_state[iDim+1] * dsM_du[iVar] );
                        JAC_I(iDim+1, iVar) -= omegaSM * dp_du_i[iVar] * norm[iDim];
                    }
                    JAC_I(iDim+1, iDim+1) += omegaSM * (sL - proj_vel_i);

                    for (jDim = 0; jDim < nDim; jDim++)
                        JAC_I(iDim+1, jDim+1) -= omegaSM * vel_i[iDim] * norm[jDim];
                }

                for (iVar = 0; iVar < nVar; iVar++)
                    JAC_I(nVar-1, iVar) = sM * ( deStar_du[iVar] + dpStar_du[iVar] ) + ( eStar + pStar ) * dsM_du[iVar];


                /*--------- jac_j ---------*/

                dsM_du[0] = ( proj_vel_j * proj_vel_j - sM * sR - 0.5 * gamma_minus_one * vel_sqr_j ) / rho_m;
                for (iDim = 0; iDim < nDim; iDim++)
                    dsM_du[iDim+1] = - ( norm[iDim] * ( 2 * proj_vel_j - sR - sM) - gamma_minus_one * vel_j[iDim] ) / rho_m;
                dsM_du[nVar-1]  = - gamma_minus_one / rho_m;

                for (iVar = 0; iVar < nVar; iVar++) {
                    dpStar_du[iVar] = rho_j * (sL - proj_vel_i) * dsM_du[iVar];
                    deStar_du[iVar] = omega * ( sM * dpStar_du[iVar] + ( eStar + pStar ) * dsM_du[iVar] );
                    JAC_J(0, iVar) = inter_state[0] * ( omegaSM + 1.0 ) * dsM_du[iVar];
                    JAC_J(nVar-1, iVar) = sM * (deStar_du[iVar] + dpStar_du[iVar]) + (eStar + pStar) * dsM_du[iVar];
                }

                for (iDim = 0; iDim < nDim; iDim++) {
                    for (jVar = 0; jVar < nVar; jVar++)
                        JAC_J(iDim+1, jVar) = ( omegaSM + 1.0 ) * ( inter_state[iDim+1] * dsM_du[jVar] + norm[iDim] * dpStar_du[jVar] );
                }
            }
        }
        else {
            if (sR < 0.0) {
                for (iVar = 0; iVar < nVar * nVar; iVar++) jac_i[iVar] = 0;
                ComputeJacobian(vel_j.data(), vel_sqr_j, proj_vel_j, norm.data(), energy_j, jac_j);
            }
            else {
                /*--- Compute Jacobian based on Right Star State ---*/

                eStar = inter_state[nVar-1];
                omega = 1/(sR-sM);
                omegaSM = omega * sM;

                /*--------- jac_i ---------*/

                dsM_du[0] = ( - proj_vel_i * proj_vel_i + sM * sL + 0.5 * gamma_minus_one * vel_sqr_i ) / rho_m;
                for (iDim = 0; iDim < nDim; iDim++)
                    dsM_du[iDim+1] = ( norm[iDim] * ( 2 * proj_vel_i - sL - sM ) - gamma_minus_one * vel_i[iDim] ) / rho_m;
                dsM_du[nVar-1] = gamma_minus_one / rho_m;

                for (iVar = 0; iVar < nVar; iVar++) {
                    dpStar_du[iVar] = rho_i * (sR - proj_vel_j) * dsM_du[iVar];
                    deStar_du[iVar] = omega * ( sM * dpStar_du[iVar] + ( eStar + pStar ) * dsM_du[iVar] );
                    JAC_I(0, iVar) = inter_state[0] * ( omegaSM + 1 ) * dsM_du[iVar];
                }

                for (iDim = 0; iDim < nDim; iDim++) {
                    for (jVar = 0; jVar < nVar; jVar++)
                        JAC_I(iDim+1, jVar) = (omegaSM + 1) * ( inter_state[iDim+1] * dsM_du[jVar] + norm[iDim] * dpStar_du[jVar] );
                }

                for (iVar = 0; iVar < nVar; iVar++)
                    JAC_I(nVar-1, iVar) = sM * (deStar_du[iVar] + dpStar_du[iVar]) + (eStar + pStar) * dsM_du[iVar];


                /*--------- Right Jacobian ---------*/

                /*--- Computing derivatives ---*/
                dp_du_i[0] = 0.5 * gamma_minus_one * vel_sqr_j;
                dsM_du[0] = - ( - proj_vel_j * proj_vel_j + sM * sR + dp_du_i[0] ) / rho_m;
                drhoStar_du[0] = omega * ( sR + inter_state[0] * dsM_du[0] );
                deStar_du[0] += omega * proj_vel_j * ( enthalpy_j - dp_du_i[0] );
                for (iDim = 0; iDim < nDim; iDim++) {
                    dp_du_i[iDim+1] = - gamma_minus_one * vel_j[iDim];
                    dsM_du[iDim+1] = - ( norm[iDim] * ( 2 * proj_vel_j - sR - sM) + dp_du_i[iDim+1] ) / rho_m;
                    drhoStar_du[iDim+1] = omega * ( - norm[iDim] + inter_state[0] * dsM_du[iDim+1] );
                    deStar_du[iDim+1] += omega * ( - norm[iDim] * enthalpy_j - proj_vel_j * dp_du_i[iDim+1] );
                }
                dp_du_i[nVar-1] = gamma_minus_one;
                dsM_du[nVar-1]  = - dp_du_i[nVar-1] / rho_m;
                drhoStar_du[nVar-1] = omega * inter_state[0] * dsM_du[nVar-1];
                deStar_du[nVar-1] += omega * ( sR - proj_vel_j - proj_vel_j * dp_du_i[nVar-1] );

                for (iVar = 0; iVar < nVar; iVar++) {
                    dpStar_du[iVar] = rho_j * (sL - proj_vel_i) * dsM_du[iVar];
                    deStar_du[iVar] = omega * ( sM * dpStar_du[iVar] + ( eStar + pStar ) * dsM_du[iVar] );
                }


                /* jac_j */
                for (iVar = 0; iVar < nVar; iVar++) {
                    JAC_J(0, iVar) = sM * drhoStar_du[iVar] + inter_state[0] * dsM_du[iVar];
                    JAC_J(nVar-1, iVar) = sM * ( deStar_du[iVar] + dpStar_du[iVar] ) + ( eStar + pStar ) * dsM_du[iVar];
                }

                for (jDim = 0; jDim < nDim; jDim++) {
                    for (iVar = 0; iVar < nVar; iVar++)
                        JAC_J(jDim+1, iVar) = ( omegaSM + 1 ) * ( norm[jDim] * dpStar_du[iVar] + inter_state[jDim+1] * dsM_du[iVar] );

                    JAC_J(jDim+1, 0) += omegaSM * vel_j[jDim] * proj_vel_j;
                    JAC_J(jDim+1, jDim+1) += omegaSM * (sR - proj_vel_j);

                    for (iDim = 0; iDim < nDim; iDim++)
                        JAC_J(jDim+1, iDim+1) -= omegaSM * vel_j[jDim] * norm[iDim];

                    for (iVar = 0; iVar < nVar; iVar++)
                        JAC_J(jDim+1, iVar) -= omegaSM * dp_du_i[iVar] * norm[jDim];
                }

            }
        }

        /* Scale Jacobians by half the area */
        for (iVar = 0; iVar < nVar * nVar; iVar++) {
            jac_i[iVar] *= 0.5 * area;
            jac_j[iVar] *= 0.5 * area;
        }
    }

#undef JAC_I
#undef JAC_J
}
