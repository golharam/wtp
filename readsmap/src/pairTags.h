#ifndef __pairTags_h__
#define __pairTags_h__ 1

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "zutil.h"
#include <math.h>
#include <zlib.h>
#include <stdarg.h>

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using std::string;
using std::vector;
using std::ofstream;
using std::ostream;
using std::ifstream;

using std::cout;
using std::cerr;

#include "zcompress.h"

//If the pairing is out of range, use this value.
#define PAIR_DIST_NA -1999999999

class PairingDistanceHistogram {
 public:
  PairingDistanceHistogram();
  PairingDistanceHistogram(int hist_min, int hist_max);
  void init(int hist_min, int hist_max);

  int addToHistogram(const vector<int> &distances);
  //distances is the list of pairing distances that come from a single beadid
  void bin_histogram(const size_t bin_size);

  void printHistogram( ostream &out ) const;
  static void printSetHistograms( const vector<PairingDistanceHistogram*> &histSet, ostream &out, const vector<int> &sumLocations );

 private:
  vector<int> clone_sizes;
  vector<int> frequencies;

  int min_pair_size;
  //This bin holds # of occurances where pair is on different chomosomes.
  size_t nan_dist_bin;
};

class PairingDistanceStat {
 public:
  PairingDistanceStat();
  PairingDistanceStat(int hist_min,int hist_max) {init(hist_min,hist_max);};
  void init(int hist_min, int hist_max);
  void set_bin_size( size_t bin_size_in) { bin_size = bin_size_in; };
  ~PairingDistanceStat();

  void add_R3_distances(const vector<int> &R3_distances);
  void add_F3_distances(const vector<int> &F3_distances);
  void add_nonRescue_distances(const vector<int> &nr_distances);
  void add_3types_distances(const vector<int> &r3_distances,
			    const vector<int> &f3_distances,
			    const vector<int> &nr_distances);
  void bin_histograms();
  void bin_histograms(const size_t bin_size);

  void addHeaderArgLine(int argc, char **argv);
  void printHistograms( ostream &out ) const;
  void setReqUniqPairing(bool value) { reqUniqPairing=value; };

 private:
  vector<PairingDistanceHistogram*> histograms;
  string headerArgLine;
  int hist_min;
  int hist_max;
  size_t bin_size;
  bool reqUniqPairing;
};

#endif
