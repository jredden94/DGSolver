#include "Quad.hpp"

Quad::Quad() { }
Quad::~Quad() { }

// 12 points
double quad_2d_degree6[] = {0.873821971016996,0.063089014491502,0.050844906370207,
    0.063089014491502,0.873821971016996,0.050844906370207,
    0.063089014491502,0.063089014491502,0.050844906370207,
    0.501426509658179,0.249286745170910,0.116786275726379,
    0.249286745170910,0.501426509658179,0.116786275726379,
    0.249286745170910,0.249286745170910,0.116786275726379,
    0.636502499121399,0.310352451033784,0.082851075618374,
    0.310352451033784,0.636502499121399,0.082851075618374,
    0.636502499121399,0.053145049844816,0.082851075618374,
    0.310352451033784,0.053145049844816,0.082851075618374,
    0.053145049844816,0.310352451033785,0.082851075618374,
    0.053145049844816,0.636502499121399,0.082851075618374};

double quad_1d_degree4[] = {-8.611363115940526e-01,3.478548451374537e-01, 
    -3.399810435848563e-01, 6.521451548625464e-01,  
    3.399810435848563e-01, 6.521451548625464e-01, 
    8.611363115940526e-01, 3.478548451374537e-01};

void Quad::GetQuadrature(vector<double> &r1, 
        vector<double> &r2, vector<double> &w, 
        vector<double> &r1d, vector<double> &w1d) {
    r1.resize(12, 0.0);
    r2.resize(12, 0.0);
    w.resize(12, 0.0);
    r1d.resize(4, 0.0);
    w1d.resize(4, 0.0);

    for (auto i = 0ul; i < 12; i++) {
        r1[i] = quad_2d_degree6[3 * i];
        r2[i] = quad_2d_degree6[3 * i + 1];
        w[i] = quad_2d_degree6[3 * i + 2] / 2.0;
    }

    for (auto i = 0ul; i < 4; i++) {
        r1d[i] = quad_1d_degree4[2 * i];
        w1d[i] = quad_1d_degree4[2 * i + 1];
    }
}
