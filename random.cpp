#include "random.h"



RandomNumbers::RandomNumbers(unsigned long int s)
:seed(s)
{
	if (s==0)
	{
		std::random_device rd;
		seed = rd();
	}
	rng = std::mt19937(seed);
}

// méthodes de distributions aléatoires classsiques

void RandomNumbers::uniform_double(std::vector<double>& vec, double lower, double upper)
{  
	std::uniform_real_distribution<> d_uni(lower,upper); 
	for (size_t n=0; n<vec.size(); n++) 
	{
		vec[n] = d_uni(rng);
	}
}
	
double RandomNumbers::uniform_double(double lower, double upper)
{
	std::uniform_real_distribution<> d_uni(lower, upper); 
	return d_uni(rng);
}
		
		
		
		
		
		
void RandomNumbers::normal(std::vector<double>& vec, double mean, double sd)
{
	std::normal_distribution<double> distribution(mean,sd);
	for (size_t n=0; n<vec.size(); n++) 
	{
		vec[n] = distribution(rng);
	}
}

double RandomNumbers::normal(double mean, double sd)
{
	std::normal_distribution<double> distribution(mean,sd);
	return distribution(rng);
}

void RandomNumbers::poisson(std::vector<int>& vec, double mean)
{
	std::poisson_distribution<int> distribution(mean);
	for (size_t n=0; n<vec.size(); n++) 
	{
		vec[n] = distribution(rng);
	}
}
	

int RandomNumbers::poisson(double mean)
{
	std::poisson_distribution<int> distribution(mean);
	return distribution(rng);
}



