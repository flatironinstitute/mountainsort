/*
 * Copyright 2016-2017 Flatiron Institute, Simons Foundation
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#ifndef MLVECTOR_H
#define MLVECTOR_H

#include <vector>
#include <QString>

template <typename T>
class MLVector : public std::vector<T> {
public:
    using std::vector<T>::vector;
    MLVector& operator<<(const T& val)
    {
        this->emplace_back(val);
        return *this;
    }
    qintptr count() const { return (qintptr)std::vector<T>::size(); } // breaks for large values
    T operator[](qintptr i) const { return std::vector<T>::operator[]((std::size_t)i); }
    T& operator[](qintptr i) { return std::vector<T>::operator[]((std::size_t)i); }
    bool isEmpty() const { return std::vector<T>::size() == 0; }
    void append(const MLVector<T>& other) { this->insert(this->end(), other.begin(), other.end()); }
    /// Witold, when you have a chance, I need the proper way to set constBegin=begin and constEnd=end for consistency with QVector. The following did not compile:
    //std::vector<T>::const_iterator constBegin() const {return begin();}
    //std::vector<T>::const_iterator constEnd() const {return end();}
    T value(qintptr i) const
    {
        if (i < 0)
            return 0;
        if (i >= this->count())
            return 0;
        return (*this)[i];
    }
};

#endif // MLVECTOR_H
