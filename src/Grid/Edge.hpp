#pragma once

#include <vector>

using namespace std;

class Edge {
    public:
        Edge();
        ~Edge();

        void SetNodes(const unsigned long n1, const unsigned long n2);
        void AddElm(const unsigned long ielm);
        unsigned long Node1(void) const;
        unsigned long Node2(void) const;
        unsigned long ElmL(void) const;
        unsigned long ElmR(void) const;
        const vector<double> Norm(void) const;
        const vector<double> Midpoint(void) const;
        double Length(void) const;
        const vector<unsigned long>& GetElms(void) const;
        void SetMidpoint(const double mx, const double my);
        void SetLength(const double len);
        void SetNorm(const double nx, const double ny);
        void SetElms(const unsigned long elmL, const unsigned long elmR);
        unsigned long EdgeNumL(void) const;
        unsigned long EdgeNumR(void) const;
        void SetEdgeNumL(unsigned long edge_num);
        void SetEdgeNumR(unsigned long edge_num);

    private:
        unsigned long n1, n2, elmL, elmR;
        unsigned short edge_numL, edge_numR; // 0, 1, 2 edge number for respective element
        vector<unsigned long> elms;

        // Unit Normal from elmL to elmR
        double nx, ny, mx, my;
        double len;
};
