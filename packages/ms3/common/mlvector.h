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
