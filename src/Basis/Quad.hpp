#pragma once

#include <vector>

using namespace std;

class Quad {
    public:
        Quad();
        ~Quad();

        void GetQuadrature(vector<double> &r1, 
                vector<double> &r2, vector<double> &w, 
                vector<double> &r1d, vector<double> &w1d);

    private:
};
