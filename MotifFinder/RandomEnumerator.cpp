/* 
 * File:   RandomEnumerator.cpp
 * Author: Eric
 * 
 * Created on October 29, 2014, 2:59 PM
 */
#include <iterator>
#include <stdlib.h>
#include <time.h>
#include "RandomEnumerator.h"

RandomEnumerator::RandomEnumerator() {
    buckets = std::vector<double>(1,0);
}

RandomEnumerator::RandomEnumerator(const RandomEnumerator& orig) {
    
    std::copy(orig.buckets.begin(), orig.buckets.end(), buckets.begin());    
}

RandomEnumerator::~RandomEnumerator() {
}

bool RandomEnumerator::InitializeDistribution(std::vector<double> & input)
{
    if(!ValidateSizes(input))
        return false;
    buckets.resize(input.size(), 0);
    double currentEndpoint = 0;
    for(int i = 0; i < input.size(); i++)
    {
        currentEndpoint += input[i];
        buckets[i] = currentEndpoint;
    }
    return true;
}

bool RandomEnumerator::ValidateSizes(std::vector<double>& input)
{
    double sum = 0;
    for(int i = 0; i < input.size(); i++)
    {
        sum += input[i];
    }
    return sum == 1;
}

std::vector<double> & RandomEnumerator::GetDistrubtionBuckets()
{
    return buckets;
}

int RandomEnumerator::EnumerateRandomVar()
{
    double toEnum = UniformOneZero();
    for(int i = 0; i < buckets.size(); i++)
        if(toEnum < buckets[i])
            return i;
    return -1;
}

double UniformOneZero()
{
    int i = rand();
    double random = (i / (double)(RAND_MAX));
    return random;
}