/* 
 * File:   SearchResult.h
 * Author: Eric
 *
 * Created on December 1, 2014, 2:56 PM
 */

#ifndef SEARCHRESULT_H
#define	SEARCHRESULT_H
#include<vector>

#include "Nucleotide.h"
using namespace std;
class SearchResult {
public:
    SearchResult();
    vector<int> loci;
    vector<Nucleotide_t> motif;
    double score;
private:

};

#endif	/* SEARCHRESULT_H */

