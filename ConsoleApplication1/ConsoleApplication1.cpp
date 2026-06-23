#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>

#define eval

#pragma omp declare target

long double xi_i(unsigned long long int n2, unsigned int i)
{
    return sqrt(2 * (n2 + i));
}

long double xi_nil(unsigned long long int n2)
{
    return sqrt(2 * n2 + 1);
}

long double frac(long double f, long double s) {
    return 1 / (f + s);
}

long double j(unsigned long long int n2, unsigned int beta1, unsigned int beta3)
{
    return (frac(xi_i(n2, beta1), xi_nil(n2)) * frac(xi_i(n2, beta1), xi_i(n2, beta3)) * frac(xi_nil(n2), xi_i(n2, beta3))) \
        + 2 * std::pow(xi_i(n2, beta3), 2) * (xi_i(n2, beta1) * (2 * xi_nil(n2) + xi_i(n2, beta3)) + xi_nil(n2)*(xi_nil(n2) + 2*xi_i(n2, beta3))) * std::pow(frac(xi_i(n2, beta1), xi_nil(n2)), 2) \
        * std::pow(frac(xi_i(n2, beta1), xi_i(n2, beta3)), 3) * std::pow(frac(xi_nil(n2), xi_i(n2, beta3)), 2) \
        - xi_i(n2, beta3) * (xi_i(n2, beta1) + xi_nil(n2) + 2 * xi_i(n2, beta3)) * frac(xi_i(n2, beta1), xi_nil(n2)) * std::pow(frac(xi_i(n2, beta1), xi_i(n2, beta3)), 2) \
        * std::pow(frac(xi_nil(n2), xi_i(n2, beta3)), 2) \
        - ((xi_i(n2, beta1) * (xi_nil(n2) + xi_i(n2, beta3)) + 2 * xi_nil(n2) * xi_i(n2, beta3)) * std::pow(frac(xi_i(n2, beta1), xi_nil(n2)), 2) * std::pow(frac(xi_i(n2, beta1), xi_i(n2, beta3)), 2) * frac(xi_nil(n2), xi_i(n2, beta3)));
}

#pragma omp end declare target

std::vector<unsigned long long int> maxnum_vec(unsigned int n) {
    std::vector<unsigned long long int> maxnum;

    maxnum.reserve(n);

    unsigned long long int k = 1;

    for (unsigned int i = 1; i <= n; i++) {
        k *= 10;
        maxnum.push_back(k);
    }

    return maxnum;
}

int main()
{
    std::vector<unsigned long long int> max_num = maxnum_vec(10);

#ifdef debug
    for ((unsigned long long int)& f : max_num) std::cout << f << std::endl;
#endif

#ifdef eval
    for (unsigned int i = 0; i < max_num.size(); i++) {
        long double total = j(0, 1, 1);

        long double sum = 0.0;
        unsigned long long int MAX = max_num[i];

#pragma omp target teams distribute parallel for reduction(+:sum)
        for (unsigned long long int n2 = 1; n2 <= MAX; n2++) {
            for (unsigned int beta1 = 0; beta1 <= 1; beta1++) {
                for (unsigned int beta3 = 0; beta3 <= 1; beta3++) {
                    sum += (n2 + 1.0) * ((beta1 + beta3) % 2 == 0 ? 1.0 : -1.0 ) * j(n2, beta1, beta3);
                }
            }
        }

        total += sum;

        std::cout << "Max number:" << MAX << std::endl;
        std::cout.precision(15);
        std::cout << "Result:" << total << std::endl;
    }
#endif
}
