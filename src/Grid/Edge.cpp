#include "Edge.hpp"

Edge::Edge() { }
Edge::~Edge() { }

void Edge::SetNodes(const unsigned long n1, const unsigned long n2) {
    this->n1 = n1; this->n2 = n2;
}

void Edge::AddElm(const unsigned long ielm) { elms.push_back(ielm); }
const vector<unsigned long>& Edge::GetElms(void) const { return elms; }

unsigned long Edge::Node1(void) const { return n1; }
unsigned long Edge::Node2(void) const { return n2; }
unsigned long Edge::ElmL(void) const { return elmL; }
unsigned long Edge::ElmR(void) const { return elmR; }

void Edge::SetMidpoint(const double mx, const double my) { this->mx = mx; this->my = my; }
void Edge::SetLength(const double len) { this->len = len; }
void Edge::SetNorm(const double nx, const double ny) { this->nx = nx; this->ny = ny; }
void Edge::SetElms(const unsigned long elmL, const unsigned long elmR) {
    this->elmL = elmL;
    this->elmR = elmR;
}

double Edge::Length(void) const { return len; }
const vector<double> Edge::Norm(void) const { return vector<double> {nx, ny, 0.0}; }
const vector<double> Edge::Midpoint(void) const { return vector<double> {mx, my, 0.0}; }
unsigned long Edge::EdgeNumL(void) const { return edge_numL; }
unsigned long Edge::EdgeNumR(void) const { return edge_numR; }
void Edge::SetEdgeNumL(unsigned long edge_num) { edge_numL = edge_num; }
void Edge::SetEdgeNumR(unsigned long edge_num) { edge_numR = edge_num; }
