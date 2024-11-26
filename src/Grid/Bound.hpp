#pragma once

#include <string>
#include <vector>
#include "Edge.hpp"

using namespace std;

class Bound {
    public:
        friend class Grid;

        enum class BCType { Inlet, Outlet, InviscidWall, ViscousWall, SupersonicInlet, SupersonicOutlet, Farfield, Empty };

        Bound();
        ~Bound();

        void SetName(const string name);
        void SetEdges(vector<Edge> edges);

        string Name(void) const;
        BCType Type(void) const;
        const vector<Edge>& GetEdges(void) const;
        vector<Edge>& GetEdges(void);

    private:
        string name;
        BCType type;
        vector<Edge> edges;
};
