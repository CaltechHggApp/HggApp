#ifndef _comparator_hh_included_
#define _comparator_hh_included_

// Some comparators to be used as predicates for std::sort.
  template<typename T>
  bool PtComparator(const T& v1, const T& v2) {return v1.Pt() < v2.Pt();}; //< Returns v1.Pt() < v2.Pt().
  template<typename T>
  bool EtComparator(const T& v1, const T& v2) {return v1.Et() < v2.Et();}; //< Return v1.Et() < v2.Et().
  template<typename T>
  bool MassComparator(const T& v1, const T& v2) {return v1.M() < v2.M();}; //< Return v1.M() < v2.M().
  template<typename T>
  bool Mod2Comparator(const T& v1, const T& v2) {return v1.Mod2() < v2.Mod2();}; //< Return v1.Mod2() < v2.Mod2().
#endif
