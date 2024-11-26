#pragma once

#include "GridReader.hpp"

class MshReader : public GridReader {
    public:
        MshReader();
        MshReader(string filename);
        ~MshReader() override;

        void ReadFile() override;

    private:
        string FindTag(const string& tag, ifstream &file) const;
};
