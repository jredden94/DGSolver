#pragma once

#include "GridReader.hpp"

class SU2Reader : public GridReader {
    public:
        SU2Reader();
        SU2Reader(string filename);
        ~SU2Reader() override;

        void ReadFile() override;

    private:
        string FindTag(const string& tag, ifstream &file) const;
};
