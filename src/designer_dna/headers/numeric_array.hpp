/*
 * BSD 3-Clause License
 *
 * Copyright (c) 2025, Spill-Tea
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this
 *    list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef NUMERIC_ARRAY_HPP
#define NUMERIC_ARRAY_HPP

#include <cstddef>
#include <stdexcept>
#include <algorithm>

template <typename T>
class NumericArray {
public:
    NumericArray(std::size_t size);
    ~NumericArray();

    T& operator[](std::size_t index);
    T get(std::size_t index);
    void set(std::size_t index, T& value);
    void fill(T& value);
    NumericArray<T> slice(std::size_t start, std::size_t end);
    std::size_t length() const;

private:
    T* data;
    std::size_t size_;
};

template <typename T>
NumericArray<T>::NumericArray(std::size_t size)
    : size_(size) {
    if (size_ == 0) throw std::invalid_argument("Size must be > 0");
    data = new T[size_];
}

template <typename T>
NumericArray<T>::~NumericArray() {
    delete[] data;
}

template <typename T>
T& NumericArray<T>::operator[](std::size_t index) {
    if (index >= size_) throw std::out_of_range("Index out of range");
    return data[index];
}

template <typename T>
T NumericArray<T>::get(std::size_t index) {
    if (index >= size_) throw std::out_of_range("Index out of range");
    return data[index];
}

template <typename T>
void NumericArray<T>::set(std::size_t index, T& value) {
    if (index >= size_) throw std::out_of_range("Index out of range");
    data[index] = value;
}

template <typename T>
void NumericArray<T>::fill(T& value) {
    for (std::size_t i = 0; i < size_; ++i)
        data[i] = value;
}

template <typename T>
NumericArray<T> NumericArray<T>::slice(std::size_t start, std::size_t end) {
    if (start > end || end > size_)
        throw std::out_of_range("Invalid slice range");
    NumericArray<T> result(end - start);
    for (std::size_t i = start; i < end; ++i)
        result.set(i - start, data[i]);
    return result;
}

template <typename T>
std::size_t NumericArray<T>::length() const {
    return size_;
}

#endif // NUMERIC_ARRAY_HPP
