#pragma once

#include <vector>

using namespace std;

class Node {
    public:
        Node();
        ~Node();

        double X() const;
        double Y() const;
        double Z() const;
        void SetXYZ(const double x, const double y, const double z = 0.0);
        void SetXYZ(vector<double> &xyz);
        void AddElm(const unsigned long elm_ind);
        const vector<unsigned long>& GetElms(void) const;

    private:
        double x, y, z;
        vector<unsigned long> elms;
};
