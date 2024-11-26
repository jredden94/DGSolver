#include "SolVector.hpp"
#include <fstream>

SolVector::SolVector() : initialized(true) { }
SolVector::SolVector(const SolVector& vec) : initialized(true) { Init(vec.nBlock, vec.nVar); }
SolVector::~SolVector() { 
    if (initialized) {
        delete[] val;
        val = nullptr;
        initialized = false;
    }
}

void SolVector::operator=(const SolVector &vec) { 
    for (auto i = 0ul; i < len; i++) val[i] = vec.val[i];
}

double* SolVector::operator[](const unsigned long &iBlock) { return &val[iBlock * nVar]; }
const double* const SolVector::operator[](const unsigned long &iBlock) const { return &val[iBlock * nVar]; }
void SolVector::operator+=(const SolVector &vec) { AddVector(vec); }
void SolVector::operator-=(const SolVector &vec) { SubtractVector(vec); }
void SolVector::operator*=(const double scalar) {
    for (auto i = 0ul; i < len; i++) val[i] *= scalar;
}
void SolVector::operator/=(const double scalar) {
    double inv_scalar = 1.0 / scalar;
    for (auto i = 0ul; i < len; i++) val[i] *= inv_scalar;
}
void SolVector::PlusEqualScalarMult(const SolVector &vec, const double scalar) {
    for (auto i = 0ul; i < len; i++) val[i] += scalar * vec.val[i];
}
void SolVector::MinusEqualScalarMult(const SolVector &vec, const double scalar) {
    for (auto i = 0ul; i < len; i++) val[i] -= scalar * vec.val[i];
}

void SolVector::Init(const unsigned long nBlock, const unsigned short nVar) { 
    if (initialized) delete[] val;
    this->nVar = nVar;
    this->nBlock = nBlock;
    len = nVar * nBlock;
    val = new double[len]; 
    initialized = true;
}

void SolVector::AddBlock(const unsigned long& iNode, const vector<double>& block, const double& alpha) {
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] += block[i] * alpha;
}

void SolVector::AddBlock(const unsigned long& iNode, const double* block, const double& alpha) {
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] += block[i] * alpha;
}

void SolVector::SetBlock(const unsigned long& iNode, const vector<double>& block, const double& alpha) {
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] = block[i] * alpha;
}
void SolVector::SetBlock(const unsigned long& iNode, const double* block, const double& alpha) {
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] = block[i] * alpha;
}

const double * const SolVector::GetBlock(const unsigned long& iNode) const { return &val[iNode * nVar]; }
double* SolVector::GetBlock(const unsigned long& iNode) { return &val[iNode * nVar]; }
double* SolVector::GetVal() { return val; }
unsigned long SolVector::GetLength() { return len; }
void SolVector::CopyBlock(const unsigned long& iBlock, double *ret_block) { 
    const double *block = GetBlock(iBlock);
    for (unsigned short i = 0; i < nVar; i++) ret_block[i] = block[i];
}

unsigned long SolVector::GetBlockCount(void) const { return nBlock; }

void SolVector::SubtractBlock(const unsigned long& iNode, const vector<double>& block, const double& alpha) { 
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] -= block[i] * alpha;
}

void SolVector::SubtractBlock(const unsigned long& iNode, const double *block, const double& alpha) { 
    for (unsigned short i = 0; i < nVar; i++) val[iNode * nVar + i] -= block[i] * alpha;
}

void SolVector::SetZeroes(void) { 
    for (auto i = 0ul; i < len; i++) val[i] = 0;
}

double SolVector::DotProd(const SolVector& vec) const { 
    double dot = 0;
    for (auto i = 0ul; i < len; i++) dot += val[i] * vec.val[i];
    return dot;
}

void SolVector::InitValues(const vector<double>& x0) {
    cout << "nVar: " << nVar << "\tinit_size: " << x0.size() << endl;
    assert(nVar == x0.size());
    for (auto i = 0ul; i < nBlock; i++) SetBlock(i, x0);
}

vector<double> SolVector::ResNorm(void) const {
    vector<double> res(nVar, 0);
    for (auto i = 0ul; i < nBlock; i++) {
        for (unsigned short j = 0; j < nVar; j++) {
            res[j] += pow(val[i * nVar + j], 2);
        }
    }

    for (unsigned short i = 0; i < nVar; i++) {
        res[i] = res[i] /  nBlock;
    }

    return res;
}
void SolVector::AddVector(const SolVector &vec, SolVector &result, const double scale) const {
    for (auto i = 0ul; i < len; i++) result.val[i] = val[i] + scale * vec.val[i];
}

void SolVector::AddVector(const SolVector &vec, const double scale) const {
    for (auto i = 0ul; i < len; i++) val[i] += scale * vec.val[i];
}

void SolVector::SubtractVector(const SolVector &vec, SolVector &result, const double scale) const {
    const double *vec_val = vec.val;
    for (auto i = 0ul; i < len; i++) result.val[i] = val[i] - scale * vec_val[i];
}

void SolVector::SubtractVector(const SolVector &vec, const double scale) {
    for (auto i = 0ul; i < len; i++) val[i] -= scale * vec.val[i];
}

double SolVector::ComputeMagnitude() const {
    double mag = 0;
    for (auto i = 0ul; i < len; i++) mag += val[i] * val[i];
    return sqrt(mag);
}

double SolVector::SquaredNorm(void) const { return DotProd(*this); }

double SolVector::Norm(void) const { return sqrt(SquaredNorm() ); }

void SolVector::WriteValues(string filename) const {
    ofstream f(filename);
    for (auto i = 0ul; i < len; i++) f << val[i] << "\n";
    f.close();
}
