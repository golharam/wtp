#include "pairTags.h"


PairingDistanceHistogram::PairingDistanceHistogram(int hist_min, int hist_max) {
  init(hist_min,hist_max);
}

void PairingDistanceHistogram::init(int hist_min, int hist_max) {
  const int startFreq = 0;
  int numSlots = hist_max - hist_min + 1;
  frequencies.reserve(numSlots);
  clone_sizes.reserve(numSlots);


  //Now populate clone_sizes and dist_frequencies vectors
  int pairDist;
  for(pairDist=hist_min; pairDist<hist_max+1; pairDist++) {
    clone_sizes.push_back(pairDist);
    frequencies.push_back(startFreq);
  }
  min_pair_size = hist_min;

  nan_dist_bin = clone_sizes.size();
  clone_sizes.push_back(PAIR_DIST_NA );
  frequencies.push_back(startFreq);
}

int PairingDistanceHistogram::addToHistogram(const vector<int> &distances) {

  //if(distances.size()==0) return -1; //Nothing to add.

  int dist_sum = sum(distances,PAIR_DIST_NA); //Any value that is PAIR_DIST_NA won't be summed
  size_t num_out_of_range = count(distances.begin(),distances.end(),PAIR_DIST_NA);

  if(num_out_of_range==distances.size()) {
    frequencies[nan_dist_bin]++;
    return PAIR_DIST_NA; //all distances are na (different chroms)
  }

  //Even if some distances are out of range, will ignore those and go with
  //other distances for the histogram.
  int pairDist_toAdd;

  if ( (dist_sum%distances.size()) == 0) {
    pairDist_toAdd = dist_sum / (distances.size() - num_out_of_range);
  } else {
    double average = dist_sum / (distances.size() - num_out_of_range);
    pairDist_toAdd = static_cast<int>(floor (average) );
  }
  int bin_num = (pairDist_toAdd - min_pair_size);
  if( bin_num<0 || (unsigned int)bin_num>nan_dist_bin)
    bin_num = nan_dist_bin;

  /*
  cout << "[(" << pairDist_toAdd << " - "
       << min_pair_size << ") = " << pairDist_toAdd - min_pair_size << "("
       << bin_num << "), " << clone_sizes[bin_num] << "]";
  cerr << " Adding 1 to freq[(" << pairDist_toAdd << " - "
       << min_pair_size << ") = " << pairDist_toAdd - min_pair_size << "("
       << bin_num << "), " << clone_sizes[bin_num] << "]\n";
  */

  frequencies[bin_num]++;
  return pairDist_toAdd;
}

void PairingDistanceHistogram::bin_histogram( const size_t bin_size ) {

  vector<int> new_dists;
  vector<int> new_freq;
  vector<int>::const_iterator distItr = clone_sizes.begin();
  vector<int>::const_iterator freqItr = frequencies.begin();

  size_t cur_bin = 1;
  int bin_min_dist_size;
  int freq_sum = 0;
  for(;freqItr!= frequencies.end();) {

    if((*distItr)==PAIR_DIST_NA) {
      new_dists.push_back(bin_min_dist_size);
      new_freq.push_back(freq_sum);
      new_dists.push_back(PAIR_DIST_NA);
      new_freq.push_back(*freqItr);
      if (distItr+1!=clone_sizes.end()) {
	cerr << "ERROR in bin_histogram.\n";
	exit(1);
      }
      break;
    }
    if(cur_bin==1)
      bin_min_dist_size = (*distItr);
    freq_sum += (*freqItr);

    if(cur_bin == bin_size) {
      new_dists.push_back(bin_min_dist_size);
      new_freq.push_back(freq_sum);
      cur_bin = 0;
      freq_sum = 0;
    }
    cur_bin++;

    distItr++;
    freqItr++;
  }
  clone_sizes = new_dists;
  frequencies = new_freq;
}

void PairingDistanceHistogram::printHistogram( ostream &out ) const {

  vector<int>::const_iterator distItr = clone_sizes.begin();
  vector<int>::const_iterator freqItr = frequencies.begin();

  for(;freqItr!= frequencies.end();) {
    out << (*distItr) << "\t" << (*freqItr) << "\n";
    distItr++;
    freqItr++;
  }
}

void PairingDistanceHistogram::printSetHistograms( const vector<PairingDistanceHistogram*> &histSet, ostream &out, const vector<int> &sumLocations ) {
  vector<PairingDistanceHistogram*>::const_iterator histItr = histSet.begin();

  int min_pair_size = (*histItr)->min_pair_size;
  size_t num_bins= (*histItr)->clone_sizes.size();

  //make sure min_pair_sizes and number of bins is the same for all elements
  for(;histItr!=histSet.end();histItr++) {
    if (min_pair_size != (*histItr)->min_pair_size ) {
      cerr << "ERROR: In PairingDistanceHistogram class, min_pair_size are not the same.\n";
      exit(1);
    } else if (num_bins != (*histItr)->clone_sizes.size()) {
      cerr << "ERROR: In PairingDistanceHistogram class, num_bins are not the same.\n";
      exit(1);
    }
  }

  size_t i;

  size_t totValColumns = histSet.size() + sumLocations.size() + 1;
  vector<int> totals(totValColumns,0);

  for(i=0;i<num_bins;i++) {
    histItr=histSet.begin();
    if((*histItr)->clone_sizes[i]==PAIR_DIST_NA) {
      out << "N/A";
    } else {
      out << (*histItr)->clone_sizes[i];
    }
    out << "\t";
    int setNum = 1;
    int setSum = 0;
    int colNum = 0;
    for(;histItr!=histSet.end();histItr++) {
      int value = (*histItr)->frequencies[i];
      setSum += value;
      out << value;
      totals[colNum++] += value;
      if(histItr<histSet.end()-1) {
        out << "\t";
        if(find(sumLocations.begin(),sumLocations.end(),setNum++)!=sumLocations.end()) {
          out << setSum << "\t";
          totals[colNum++] += setSum;
          setSum = 0;
        }
      }
    }

    out << "\t" << setSum;
    totals[colNum++] += setSum;
    out << "\n";
  }
  out << "Totals:\t";
  printVector(totals,"\t",out);
  out << "\n";
}

PairingDistanceStat::PairingDistanceStat() {
  bin_size = 0;
  reqUniqPairing=false;
}

PairingDistanceStat::~PairingDistanceStat() {
  vector<PairingDistanceHistogram*>::iterator histItr = histograms.begin();
  for(;histItr<histograms.end();histItr++) {
    delete (*histItr);
  }
}

void PairingDistanceStat::init(int hist_min, int hist_max) {
  int numSets = 4;
  //The 4 sets are
  //R3, F3, Non-res, and Comb for all tags considered

  int setNum;
  for(setNum=0;setNum<numSets;setNum++) {
    PairingDistanceHistogram *pdh = new PairingDistanceHistogram(hist_min, hist_max);
    histograms.push_back(pdh);
  }
}

void PairingDistanceStat::addHeaderArgLine(int argc, char **argv) {
  headerArgLine = "# ";
  int i;
  for (i = 0; i <argc; i++) {
    headerArgLine += argv[i];
    if(i<argc-1) headerArgLine += " ";
  }
  headerArgLine += "\n";
}

void PairingDistanceStat::printHistograms( ostream &out) const {
  vector<int> sumLoc;  //Empty list means will just total all columns

  //Header lines
  out << headerArgLine;
  out << "#pair_dist\t" << "R3_rescue_dist\t" << "F3_rescue_dist\t"
      << "non_rescue_dist\t" << "Combination\t" << "Total\n";

  PairingDistanceHistogram::printSetHistograms(histograms,out,sumLoc);
}

void PairingDistanceStat::add_R3_distances(const vector<int> &r3_distances) {
  histograms[0]->addToHistogram(r3_distances);
}

void PairingDistanceStat::add_F3_distances(const vector<int> &f3_distances) {
  histograms[1]->addToHistogram(f3_distances);
}

void PairingDistanceStat::add_nonRescue_distances(const vector<int> &nr_distances) {
  histograms[2]->addToHistogram(nr_distances);
}

void PairingDistanceStat::add_3types_distances(const vector<int> &r3_distances,
					       const vector<int> &f3_distances,
					       const vector<int> &nr_distances) {
  size_t r3s = r3_distances.size();
  size_t f3s = f3_distances.size();
  size_t nrs = nr_distances.size();  //not rescued

  if(reqUniqPairing) {
    if(r3s==1 && f3s==0 && nrs==0) {
      histograms[0]->addToHistogram(r3_distances);
    } else if(r3s==0 && f3s==1 && nrs==0) {
      histograms[1]->addToHistogram(f3_distances);
    } else if(r3s==0 && f3s==0 && nrs==1) {
      histograms[2]->addToHistogram(nr_distances);
    } else if(r3s==0 && f3s==0 && nrs==0) {
      //do nothing
    } else {
      //Combinations by definition don't have unique pairings
    }
  } else {  //This was the only option originally (CL v4.0r2.0)
            //Most non-uniqs here will be in the comb bucket,
            //but if they are all one type, they'll occur in other buckets as well.
    if(r3s>0 && f3s==0 && nrs==0) {
      histograms[0]->addToHistogram(r3_distances);
    } else if(r3s==0 && f3s>0 && nrs==0) {
      histograms[1]->addToHistogram(f3_distances);
    } else if(r3s==0 && f3s==0 && nrs>0) {
      histograms[2]->addToHistogram(nr_distances);
    } else if(r3s==0 && f3s==0 && nrs==0) {
      //do nothing
    } else {
      vector<int> comb;
      comb.insert(comb.end(),r3_distances.begin(),r3_distances.end());
      comb.insert(comb.end(),f3_distances.begin(),f3_distances.end());
      comb.insert(comb.end(),nr_distances.begin(),nr_distances.end());
      histograms[3]->addToHistogram(comb);
    }
  }
}

void PairingDistanceStat::bin_histograms() {

  const size_t estNumBins = 500;
  if(bin_size==0) {
    //automatically determins bin_size if not defined
    unsigned int pairingRange = abs(hist_max - hist_min);
    bin_size = pairingRange / estNumBins;
    bin_size++;
  }
  bin_histograms(bin_size);
}

void PairingDistanceStat::bin_histograms(const size_t bin_size_in) {
  vector<PairingDistanceHistogram*>::iterator histItr = histograms.begin();
  for(;histItr<histograms.end();histItr++) {
    (*histItr)->bin_histogram(bin_size_in);
  }
}

int main_old_zzzz() {
  //test function
  //PairingDistanceHistogram ph(0,120),ph2(0,120);
  PairingDistanceStat ph(0,121);
  vector<int> d1,d2;
  d1.push_back(101);d1.push_back(102);d1.push_back(103);
  d1.push_back(104);d1.push_back(105);

  d2.push_back(1); d2.push_back(3); d2.push_back(9);

  ph.add_R3_distances(d1);
  ph.add_R3_distances(d2);
  //ph.printHistogram(cout);

  vector<int> d3,d4,d5,d6,d7;
  d3.push_back(0);d3.push_back(3); d3.push_back(9);
  d4.push_back(2);d4.push_back(5);
  d5.push_back(119);d5.push_back(PAIR_DIST_NA);
  d6.push_back(PAIR_DIST_NA);
  d7.push_back(120);

  ph.add_F3_distances(d3);
  ph.add_F3_distances(d4);
  ph.add_nonRescue_distances(d5);
  ph.add_nonRescue_distances(d6);
  ph.add_nonRescue_distances(d7);

  //vector<PairingDistanceHistogram*> histVect;
  //histVect.push_back(&ph);
  //histVect.push_back(&ph2);

  //vector<int> sumLoc;
  //sumLoc.push_back(1);
  //PairingDistanceHistogram::printSetHistograms(histVect,cout,sumLoc);

  //ph.bin_histograms(5);
  ph.printHistograms(cout);

  return 0;
}
