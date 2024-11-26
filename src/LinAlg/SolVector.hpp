#pragma once

#include <cassert>
#include <cmath>
#include "../Common/Config.hpp"

class SolVector {
    public:
        SolVector();
        SolVector(const SolVector&);
        ~SolVector();

        void operator=(const SolVector&);
        double* operator[](const unsigned long&);
        const double* const operator[](const unsigned long&) const;
        void operator+=(const SolVector &vec);
        void operator-=(const SolVector &vec);
        void operator*=(const double scalar);
        void operator/=(const double scalar);
        void PlusEqualScalarMult(const SolVector &vec, const double scalar);
        void MinusEqualScalarMult(const SolVector &vec, const double scalar);
        void swap(SolVector &vec);

        void Init(const unsigned long nBlock, const unsigned short nVar);
        void SetBlock(const unsigned long&, const vector<double>&, const double& = 1);
        void SetBlock(const unsigned long&, const double*, const double& = 1);
        void AddBlock(const unsigned long&, const vector<double>&, const double& = 1);
        void AddBlock(const unsigned long&, const double*, const double& = 1);
        void SubtractBlock(const unsigned long&, const vector<double>&, const double& = 1);
        void SubtractBlock(const unsigned long&, const double*, const double& = 1);
        const double* const GetBlock(const unsigned long&) const;
        double* GetBlock(const unsigned long&);
        double* GetVal(void);
        unsigned long GetLength(void);
        void CopyBlock(const unsigned long& iBlock, double *block);
        unsigned long GetBlockCount(void) const;
        void SetZeroes(void);
        double DotProd(const SolVector&) const;
        void InitValues(const vector<double>& x0);
        vector<double> ResNorm(void) const;
        void AddVector(const SolVector &vec, SolVector &result, const double scale = 1) const;
        void AddVector(const SolVector &vec, const double scale = 1) const;
        void SubtractVector(const SolVector &vec, SolVector &result, const double scale = 1) const;
        void SubtractVector(const SolVector &vec, const double scale = 1);
        double ComputeMagnitude(void) const;
        double SquaredNorm(void) const;
        double Norm(void) const;

        void WriteValues(string filename) const;

    private:
        double *val = nullptr;
        unsigned short nVar;
        unsigned long nBlock;
        unsigned long len;
        bool initialized;
};
