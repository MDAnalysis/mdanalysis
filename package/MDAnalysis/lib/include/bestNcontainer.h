#ifndef DISTOPIA_BESTNCONTAINER_H
#define DISTOPIA_BESTNCONTAINER_H

#include <map>
#include <utility>  // pair
#include <float.h>  // FLT_MAX

/*
 * Stores the N nearest *j* to *i*
 */
template<typename T, typename J>
struct BestNContainer {
    int n_results;
    std::map<std::pair<J, T>, T> results;
    std::map<T, J> existing;

    BestNContainer() {};
    BestNContainer(int n) : n_results(n) {};

    int size() const {
      return results.size();
    }

    J current_max () const {
      typename std::map<std::pair<J, T>, T>::const_iterator it;

      if (size() < n_results) {
        return FLT_MAX;
      }

      it = results.end();
      it--;

      return it->first.first;
    }

    // report a hit between *i* and *j* at distance *d*
    // returns if this hit was added to the best
    bool add_hit(J d, T i, T j) {
      // if this hit is worse than our current max, ignore it
      if (d > current_max())
        return false;

      // check if we've previously seen a hit involving *j*
      typename std::map<T, J>::iterator existing_it = existing.find(j);
      if (existing_it != existing.end()) {
        // if this hit is worse than the existing hit for *j*, ignore
        if (d > existing_it->second)
          return false;
        // otherwise
        results.erase({existing_it->second, j});
        existing.erase(existing_it);
      }

      results.insert({{d, j}, i});
      existing.insert({j, d});

      if (results.size() > n_results) {
        // too many results, chuck the largest
        typename std::map<std::pair<J, T>, T>::const_iterator result_it;
        result_it = results.end();
        result_it--;
        existing.erase(result_it->first.second);
        results.erase(result_it);
      }

      return true;
    }
};

#endif //DISTOPIA_BESTNCONTAINER_H
