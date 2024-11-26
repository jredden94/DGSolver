#include "Bound.hpp"

Bound::Bound() { }
Bound::~Bound() { }

void Bound::SetName(const string name) { this->name = name; }
void Bound::SetEdges(vector<Edge> edges) { swap(this->edges, edges); }

string Bound::Name(void) const { return name; }
Bound::BCType Bound::Type(void) const { return type; }
const vector<Edge>& Bound::GetEdges(void) const { return edges; }
vector<Edge>& Bound::GetEdges(void) { return edges; }
