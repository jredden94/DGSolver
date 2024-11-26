#include "Convection.hpp"
#include <stdexcept>

Convection::Convection() { }
Convection::~Convection() { 
    delete[] u_i;
    delete[] u_j;
    delete[] v_i;
    delete[] v_j;
    if (implicit) {
        delete[] jac_i;
        delete[] jac_j;
    }
}
Convection::Convection(Config::ConvFlux) { 
    config = &Config::GetConfig();
    gamma = config->GetGamma();
    gamma_minus_one = this->gamma - 1;
    //nVar = config->GetNumVars();
    //nDim = config->GetNumDims();
    implicit = config->IsImplicit();
    u_i = new double[nVar];
    u_j = new double[nVar];
    v_i = new double[nVar];
    v_j = new double[nVar];
    eigen.resize(nVar);
    flux_i.resize(nVar);
    flux_j.resize(nVar);
    flux.resize(nVar);
    vel_i.resize(nDim);
    vel_j.resize(nDim);
    vel_roe.resize(nDim);
    norm.resize(3);

    if (implicit) {
        jac_i = new double[nVar * nVar];
        jac_j = new double[nVar * nVar];
    }
}

void Convection::SetStates(const vector<double> &U_i, 
        const vector<double> &U_j, 
        const vector<double> &areaVector, const double &area_mag) {
    for (unsigned short i = 0; i < nVar; i++) {
        u_i[i] = U_i[i];
        u_j[i] = U_j[i];
    }
    norm = areaVector;
    area = area_mag;
}

void Convection::SetStates(const double* v_i, const double* v_j, 
        const double* areaVector, const double &area_mag) {

    for (unsigned short i = 0; i < nVar; i++) {
        this->v_i[i] = v_i[i];
        this->v_j[i] = v_j[i];
    }

    norm[0] = areaVector[0];
    norm[1] = areaVector[1];
    norm[2] = areaVector[2];
    area = area_mag;
}

const vector<double>& Convection::Flux() { return flux; }
const double* Convection::JacI() { return jac_i; }
const double* Convection::JacJ() { return jac_j; }
const double& Convection::MaxWaveSpeed() { return maxWaveSpeed; }
const double Convection::MaxLambda() { return max_lambda; }

unique_ptr<Convection> Convection::CreateConvFlux(Config *config) {
    const Config::ConvFlux &fluxScheme = config->GetConvFluxScheme();
    switch (fluxScheme) {
        case (Config::ConvFlux::Roe) : return make_unique<Roe>(); break;
        case (Config::ConvFlux::HLLC) : return make_unique<HLLC>(); break;
        default : throw std::length_error("Unrecognized flux scheme");
    }
}

void Convection::ComputeFlux(const double &rho, const double *vel, const double &enthalpy,
        const double *norm, const double *flux) {
}

void Convection::ComputeJacobian(const double *vel, const double &vel_sqr, 
        const double &proj_vel, const double *norm, const double &energy, double *jac, const double scale) {

    double phi = 0.5 * gamma_minus_one * vel_sqr;
    double a1 = gamma * energy - phi;

    if (nDim == 3) {
        jac[0] = 0.0;
        jac[1] = scale * norm[0];
        jac[2] = scale * norm[1];
        jac[3] = scale * norm[2];
        jac[4] = 0.0;

        jac[5] = scale * (norm[0] * phi - vel[0] * proj_vel);
        jac[6] = scale * (vel[0] * norm[1] - gamma_minus_one * vel[1] * norm[0] + proj_vel); // + scale * proj_vel;
        jac[7] = scale * (vel[0] * norm[2] - gamma_minus_one * vel[2] * norm[0]);
        jac[8] = scale * (vel[0] * norm[3] - gamma_minus_one * vel[3] * norm[0]);
        jac[9] = scale * gamma_minus_one * norm[0];

        jac[10] = scale * (norm[1] * phi - vel[1] * proj_vel);
        jac[11] = scale * (vel[1] * norm[1] - gamma_minus_one * vel[1] * norm[1]);
        jac[12] = scale * (vel[1] * norm[2] - gamma_minus_one * vel[2] * norm[1] + proj_vel); // + scale * proj_vel;
        jac[13] = scale * (vel[1] * norm[3] - gamma_minus_one * vel[3] * norm[1]);
        jac[14] = scale * gamma_minus_one * norm[1];

        jac[15] = scale * (norm[2] * phi - vel[2] * proj_vel);
        jac[16] = scale * (vel[2] * norm[1] - gamma_minus_one * vel[1] * norm[2]);
        jac[17] = scale * (vel[2] * norm[2] - gamma_minus_one * vel[2] * norm[2]);
        jac[18] = scale * (vel[2] * norm[3] - gamma_minus_one * vel[3] * norm[2] + proj_vel); // + scale * proj_vel;
        jac[19] = scale * gamma_minus_one * norm[2];

        jac[20] = scale * proj_vel * (phi - a1);
        jac[21] = scale * (norm[0] * a1 - gamma_minus_one * vel[0] * proj_vel);
        jac[22] = scale * (norm[1] * a1 - gamma_minus_one * vel[1] * proj_vel);
        jac[23] = scale * (norm[2] * a1 - gamma_minus_one * vel[2] * proj_vel);
        jac[24] = scale * gamma * proj_vel;
    }
    else {
        jac[0] = 0.0;
        jac[1] = scale * norm[0];
        jac[2] = scale * norm[1];
        jac[3] = 0.0;

        jac[4] = scale * (norm[0] * phi - vel[0] * proj_vel);
        jac[5] = scale * (vel[0] * norm[1] - gamma_minus_one * vel[1] * norm[0] + proj_vel); // + scale * proj_vel;
        jac[6] = scale * (vel[0] * norm[2] - gamma_minus_one * vel[2] * norm[0]);
        jac[7] = scale * gamma_minus_one * norm[0];

        jac[8] = scale * (norm[1] * phi - vel[1] * proj_vel);
        jac[9] = scale * (vel[1] * norm[1] - gamma_minus_one * vel[1] * norm[1]);
        jac[10] = scale * (vel[1] * norm[2] - gamma_minus_one * vel[2] * norm[1] + proj_vel); // + scale * proj_vel;
        jac[11] = scale * gamma_minus_one * norm[1];

        jac[12] = scale * proj_vel * (phi - a1);
        jac[13] = scale * (norm[0] * a1 - gamma_minus_one * vel[0] * proj_vel);
        jac[14] = scale * (norm[1] * a1 - gamma_minus_one * vel[1] * proj_vel);
        jac[15] = scale * gamma * proj_vel;
    }
}

void Convection::ComputePTensor(const double &rho, const double *vel,
        const double &soundSpeed, const double *norm, double *p_tensor) {
    double sqvel, rhooc, rhoxc;
    rhooc = rho / soundSpeed;
    rhoxc = rho * soundSpeed;

    /*
    if (nDim == 3) {
        sqvel = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];

        p_tensor[0] = norm[0];
        p_tensor[1] = norm[1];
        p_tensor[2] = norm[2];
        p_tensor[3] = 0.5 * rhooc;
        p_tensor[4] = 0.5 * rhooc;

        p_tensor[5] = vel[0] * norm[0];
        p_tensor[6] = vel[0] * norm[1] - rho * norm[2];
        p_tensor[7] = vel[0] * norm[2] + rho * norm[1];
        p_tensor[8] = 0.5 * (vel[0] * rhooc + rho * norm[0]);
        p_tensor[9] = 0.5 * (vel[0] * rhooc - rho * norm[0]);

        p_tensor[10] = vel[1] * norm[0] + rho * norm[2];
        p_tensor[11] = vel[1] * norm[1];
        p_tensor[12] = vel[1] * norm[2] - rho * norm[0];
        p_tensor[13] = 0.5 * (vel[1] * rhooc + rho * norm[1]);
        p_tensor[14] = 0.5 * (vel[1] * rhooc - rho * norm[1]);

        p_tensor[15] = vel[2] * norm[0] - rho * norm[1];
        p_tensor[16] = vel[2] * norm[1] + rho * norm[0];
        p_tensor[17] = vel[2] * norm[2];
        p_tensor[18] = 0.5 * (vel[2] * rhooc + rho * norm[2]);
        p_tensor[19] = 0.5 * (vel[2] * rhooc - rho * norm[2]);

        p_tensor[20] = 0.5 * sqvel * norm[0] + rho * vel[1] * norm[2] - rho * vel[2] * norm[1];
        p_tensor[21] = 0.5 * sqvel * norm[1] - rho * vel[0] * norm[2] + rho * vel[2] * norm[0];
        p_tensor[22] = 0.5 * sqvel * norm[2] + rho * vel[0] * norm[1] - rho * vel[1] * norm[0];
        p_tensor[23] = 0.5 * (0.5 * sqvel * rhooc + rho * (vel[0] * norm[0] + vel[1] * norm[1] + vel[2] * norm[2])
                + rhoxc / gamma_minus_one);
        p_tensor[24] = 0.5 * (0.5 * sqvel * rhooc - rho * (vel[0] * norm[0] + vel[1] * norm[1] + vel[2] * norm[2])
                + rhoxc / gamma_minus_one);
    }
    else {
    */
        sqvel = vel[0] * vel[0] + vel[1] * vel[1];

        p_tensor[0] = 1.0;
        p_tensor[1] = 0.0;
        p_tensor[2] = 0.5 * rhooc;
        p_tensor[3] = 0.5 * rhooc;

        p_tensor[4] = vel[0];
        p_tensor[5] = rho * norm[1];
        p_tensor[6] = 0.5 * (vel[0] * rhooc + norm[0] * rho);
        p_tensor[7] = 0.5 * (vel[0] * rhooc - norm[0] * rho);

        p_tensor[8] = vel[1];
        p_tensor[9] = -rho * norm[0];
        p_tensor[10] = 0.5 * (vel[1] * rhooc + norm[1] * rho);
        p_tensor[11] = 0.5 * (vel[1] * rhooc - norm[1] * rho);

        p_tensor[12] = 0.5 * sqvel;
        p_tensor[13] = rho * vel[0] * norm[1] - rho * vel[1] * norm[0];
        p_tensor[14] = 0.5 * (0.5 * sqvel * rhooc + rho * vel[0] * norm[0] + rho * vel[1] * norm[1] + rhoxc / gamma_minus_one);
        p_tensor[15] = 0.5 * (0.5 * sqvel * rhooc - rho * vel[0] * norm[0] - rho * vel[1] * norm[1] + rhoxc / gamma_minus_one);
    //}
}

void Convection::ComputeInversePTensor(const double &rho, const double *vel, 
        const double &soundspeed, const double *norm, double *inv_p) {
    double rhoxc, c2, gm1, k0orho, k1orho, gm1_o_c2, gm1_o_rhoxc, sqvel;
    rhoxc = rho * soundspeed;
    c2 = soundspeed * soundspeed;
    gm1 = gamma_minus_one;
    k0orho = norm[0] / rho;
    k1orho = norm[1] / rho;
    gm1_o_c2 = gm1 / c2;
    gm1_o_rhoxc = gm1 / rhoxc;

    /*
    if (nDim == 3) {
        sqvel = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];

        inv_p[0] = norm[0] - norm[2] * vel[1] / rho + norm[1] * vel[2] / rho - norm[0] * 0.5 * gm1 * sqvel / c2;
        inv_p[1] = norm[0] * gm1 * vel[0] / c2;
        inv_p[2] = norm[2] / rho + norm[0] * gm1 * vel[1] / c2;
        inv_p[3] = -norm[1] / rho + norm[0] * gm1 * vel[2] / c2;
        inv_p[4] = -norm[0] * gm1 / c2;

        inv_p[5] = norm[1] + norm[2] * vel[0] / rho - norm[0] * vel[2] / rho - norm[1] * 0.5 * gm1 * sqvel / c2;
        inv_p[6] = -norm[2] / rho + norm[1] * gm1 * vel[0] / c2;
        inv_p[7] = norm[1] * gm1 * vel[1] / c2;
        inv_p[8] = norm[0] / rho + norm[1] * gm1 * vel[2] / c2;
        inv_p[9] = -norm[1] * gm1 / c2;

        inv_p[10] = norm[2] - norm[1] * vel[0] / rho + norm[0] * vel[1] / rho - norm[2] * 0.5 * gm1 * sqvel / c2;
        inv_p[11] = norm[1] / rho + norm[2] * gm1 * vel[0] / c2;
        inv_p[12] = -norm[0] / rho + norm[2] * gm1 * vel[1] / c2;
        inv_p[13] = norm[2] * gm1 * vel[2] / c2;
        inv_p[14] = -norm[2] * gm1 / c2;

        inv_p[15] = -(norm[0] * vel[0] + norm[1] * vel[1] + norm[2] * vel[2]) / rho + 0.5 * gm1 * sqvel / rhoxc;
        inv_p[16] = norm[0] / rho - gm1 * vel[0] / rhoxc;
        inv_p[17] = norm[1] / rho - gm1 * vel[1] / rhoxc;
        inv_p[18] = norm[2] / rho - gm1 * vel[2] / rhoxc;
        inv_p[19] = gm1 / rhoxc;

        inv_p[20] = (norm[0] * vel[0] + norm[1] * vel[1] + norm[2] * vel[2]) / rho + 0.5 * gm1 * sqvel / rhoxc;
        inv_p[21] = -norm[0] / rho - gm1 * vel[0] / rhoxc;
        inv_p[22] = -norm[1] / rho - gm1 * vel[1] / rhoxc;
        inv_p[23] = -norm[2] / rho - gm1 * vel[2] / rhoxc;
        inv_p[24] = gm1 / rhoxc;
    }
    else {
    */
        sqvel = vel[0] * vel[0] + vel[1] * vel[1];

        inv_p[0] = 1.0 - 0.5 * gm1_o_c2 * sqvel;
        inv_p[1] = gm1_o_c2 * vel[0];
        inv_p[2] = gm1_o_c2 * vel[1];
        inv_p[3] = -gm1_o_c2;

        inv_p[4] = -k1orho * vel[0] + k0orho * vel[1];
        inv_p[5] = k1orho;
        inv_p[6] = -k0orho;
        inv_p[7] = 0.0;

        inv_p[8] = -k0orho * vel[0] - k1orho * vel[1] + 0.5 * gm1_o_rhoxc * sqvel;
        inv_p[9] = k0orho - gm1_o_rhoxc * vel[0];
        inv_p[10] = k1orho - gm1_o_rhoxc * vel[1];
        inv_p[11] = gm1_o_rhoxc;

        inv_p[12] = k0orho * vel[0] + k1orho * vel[1] + 0.5 * gm1_o_rhoxc * sqvel;
        inv_p[13] = -k0orho - gm1_o_rhoxc * vel[0];
        inv_p[14] = -k1orho - gm1_o_rhoxc * vel[1];
        inv_p[15] = gm1_o_rhoxc;
    //}
}
