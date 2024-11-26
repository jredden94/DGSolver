#pragma once
#include <cmath>
#include <vector>
#include <fstream>

using namespace std;

class Basis {
    public:
        Basis();
        ~Basis();

        void Preval(const vector<double> &r1, const vector<double> &r2, 
                const vector<double> &w, const vector<double> &r1d, const vector<double> &w1d);


        inline double Phi(const unsigned short i_p, const unsigned short iQuad){
            return basis[n_quad * i_p + iQuad];
        }

        inline double PhiGradX(const unsigned short i_p, const unsigned short iQuad) {
            return basis_x[n_quad * i_p + iQuad];
        }

        inline double PhiGradY(const unsigned short i_p, const unsigned short iQuad) {
            return basis_y[n_quad * i_p + iQuad];
        }

        inline double PhiSide(const unsigned short side, const unsigned short i_p, const unsigned short iQuad) {
            return basis_side[side * n_quad_1d * n_p + i_p * n_quad_1d + iQuad];
        }

    private:
        unsigned short n = 3;
        unsigned short n_p = (n+1) * (n+2) / 2;
        unsigned short n_quad = 12;
        unsigned short n_quad_1d = 4;

        double *basis, *basis_x, *basis_y, *basis_side;

        double EvalPhi(double r, double s, int n) const;
        double Phi_x(double r, double s, int n) const;
        double Phi_y(double r, double s, int n) const;
};
