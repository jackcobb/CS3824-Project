/* 
 * File:   RandomEnumerator.h
 * Author: Eric
 *
 * Created on October 29, 2014, 2:59 PM
 */

#ifndef RANDOMENUMERATOR_H
#define	RANDOMENUMERATOR_H
#include <vector>
class RandomEnumerator {
public:
    RandomEnumerator();
    RandomEnumerator(std::vector<double> & sizes);
    RandomEnumerator(const RandomEnumerator& orig);
    bool InitializeDistribution(std::vector<double> & sizes);
    int NumberDistrubtionBuckets();
    std::vector<double> & GetDistrubtionBuckets();
    int EnumerateRandomVar();
    ~RandomEnumerator();
private:
    bool ValidateSizes(std::vector<double> & sizes);
    std::vector<double> buckets;
};
double UniformOneZero();

#endif	/* RANDOMENUMERATOR_H */

