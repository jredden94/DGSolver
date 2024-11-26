#include "Node.hpp"

Node::Node() : x(0.0), y(0.0), z(0.0) { }
Node::~Node() { }

double Node::X() const { return x; }
double Node::Y() const { return y; }
double Node::Z() const { return z; }

void Node::SetXYZ(const double x, const double y, const double z) {
    this->x = x; this->y = y; this->z = z;
}

void Node::SetXYZ(vector<double> &xyz) {
    this->x = xyz[0]; this->y = xyz[1]; this->z = xyz[2];
}

void Node::AddElm(const unsigned long elm_ind) { elms.push_back(elm_ind); }

const vector<unsigned long>& Node::GetElms(void) const { return elms; }
