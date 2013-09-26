#ifndef _RANGE_H
#define _RANGE_H

#include <iostream>
#include <iterator>

//for using python style range() constructs.  From
//https://github.com/mfontanini/Programs-Scripts/blob/master/c%2B%2B_stuff/range.cpp

template<class T>
class range_iterator : public std::iterator<std::input_iterator_tag, T>{
public:
    range_iterator(const T &item) : item(item) {}
    
    // Dereference, returns the current item.
    const T &operator*() {
        return item;
    }
    
    // Prefix.
    range_iterator<T> &operator++() {
        ++item;
        return *this;
    }
    
    // Postfix.
    range_iterator<T> &operator++(int) {
        range_iterator<T> range_copy(*this);
        ++item;
        return range_copy;
    }
    
    // Compare internal item
    bool operator==(const range_iterator<T> &rhs) {
        return item == rhs.item;
    }
    
    // Same as above
    bool operator!=(const range_iterator<T> &rhs) {
        return !(*this == rhs);
    }
private:
    T item;
};

template<class T>
class range_wrapper {
public:
    range_wrapper(const T &r_start, const T &r_end) 
    : r_start(r_start), r_end(r_end) {}

    range_iterator<T> begin() {
        return {r_start};
    }
    
    range_iterator<T> end() {
        return {r_end};
    }
private:
    T r_start, r_end;
};

// Returns a range_wrapper<T> containing the range [start, end)
template<class T>
range_wrapper<T> range(const T &start, const T &end) {
    return {start, end};
}

// Returns a range_wrapper<T> containing the range [T(), end)
template<class T>
range_wrapper<T> range(const T &end) {
    return {T(), end};
}

#endif
