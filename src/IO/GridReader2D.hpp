#pragma once

#include "GridReader.hpp"

// Reads .grid files
class GridReader2D : public GridReader {
    public:
        GridReader2D();
        GridReader2D(string filename);
        ~GridReader2D() override;

        void ReadFile() override;

    private:
};
