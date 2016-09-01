#ifndef EFANNA_MATRIX_H
#define EFANNA_MATRIX_H

#include <cstddef>
#include <algorithm>
#include <stdexcept>
#include <vector>
#include <iostream>
#include "distance.hpp"
namespace efanna {

template <typename T>
class Matrix {
public:
    /**
     * Create a matrix using number of rows and cols with rawdata.
     * rawdata should be compact without any alignments.
     * Notice that efnn::Matrix itself just keep a reference of data and
     * does not make a copy of data, so keep it available.
     */
    Matrix(const size_t rows, const size_t cols, const void* data):
        rows_(rows), cols_(cols) {
        size_t sse_cols = (cols + 3)/4*4;//re align to sse format
        for (size_t i = 0; i < rows; i++) {
            row_pointers_.push_back(reinterpret_cast<const T*>(data) + (sse_cols * i));
        }
    }

    size_t get_cols() const {
        return cols_;
    }

    size_t get_rows() const {
        return rows_;
    }

    const T* get_row(const size_t index) const {
        if (index >= rows_) {
            throw std::runtime_error("index out of range");
        }
        return row_pointers_[index];
    }

    // Debug usage only
    std::vector<std::pair<T, size_t> > brute_force_search(size_t idx, size_t k, Distance<T>* distance) const {
        printf("idx: %lu\n", idx);
        std::vector<std::pair<T, size_t> > result;
        for (size_t i = 0; i < rows_; i++) {
            result.push_back(std::make_pair(
                    distance->compare(get_row(i), get_row(idx), cols_),
                    i));
        }
        std::partial_sort(result.begin(), result.begin() + k, result.end());
        result.resize(k);
        return result;
    }
private:
    size_t rows_, cols_;
    std::vector<const T*> row_pointers_;
};

}
#endif
