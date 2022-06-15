#include "iterators.h"
#include <iostream>
int main() {
auto Ag = _AtomGroupIterator(10);

auto Ar = _ArrayIterator(10);

int64_t a[10] = {0};

for( int i=0; i <10; i++){
    a[i] = i;
    std::cout << a[i] << " \n";
}

Ag.copy_ix(a);

std::cout << " ix \n";

for( int i=0; i < 10; i++){
    std::cout << Ag.ix[i] << " \n";
}

}