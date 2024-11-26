#pragma once

#include <vector>
#include <iostream>

using namespace std;

class Elm {
    public:
        Elm();
        ~Elm();

        void SetNodes(vector<unsigned long> &nodes);
        void SetCentroid(const double xc, const double yc, const double zc);
        const vector<unsigned long>& GetNodes(void) const;
        const vector<double>& GetCentroid(void) const;
        void AddEdge(const unsigned long iedge, const unsigned long inode1, const unsigned long inode2);
        void SetArea(const double area);
        const double Area(void) const;
        const vector<unsigned long>& Edges(void) const;

        unsigned short GetEdgeNum(const unsigned long n1, const unsigned long n2) const;

    private:
        vector<unsigned long> nodes, edges;
        vector<double> ctr;
        double area;
};
