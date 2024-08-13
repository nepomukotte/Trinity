//
// test_random.cpp
//
//  Created on: 16 dic 2020
//      Author: cafagna
//
#include <random/flat_random.hpp>
#include <random/normal_random.hpp>
#include <iostream>
#include <vector>



template <typename Vec> void print_ascii(const Vec & v, std::ostream &o=std::cout){
	for(int i=0; i!=v.size(); ++i) {
		std::cout << i << '\t';
		for(int j=0; j!=v[i]; ++j) std::cout << '*';
		std::cout << std::endl;
	}

}

namespace fcr=fcutils::random;

int main(){
	fcr::rand_int rnd{0,4};
	std::vector<int> histo(5);
	for(int i=0;i!=200;++i) ++histo[rnd()];

	std::cout << "The int uniform distribution in range [" << rnd.get_low() <<','<< rnd.get_high() << "] \n\n";
	print_ascii(histo);

	fcr::rand_double rndd{0.,4.5};
	std::vector<double> histod(6);
	for(int i=0;i!=200;++i) ++histod[static_cast<int>(rndd())];

	std::cout << "The double uniform distribution in range [" << rndd.get_low() <<','<< rndd.get_high() << "] \n\n";
	print_ascii(histod);

	fcr::grand grnd{15.,4.};
	std::vector<int> ghisto(30);
	for(int i=0,j=0;i!=500;++i, j=grnd()) if(j<ghisto.size())++ghisto[j]; else std::cout <<"Ooooppsss, out of range :" << j << '\n';
	std::cout << "The normal distribution with mean: " << grnd.get_mean() << ", and sigma: "<< grnd.get_sigma() << " \n\n";
	print_ascii(ghisto);
	return 0;
}




