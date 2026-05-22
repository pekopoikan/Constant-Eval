#include <iostream>
#include <cmath>
#include <vector>

#define ui unsigned int

double xi_i(ui n2, ui i)
{
    return std::sqrt(2 * (n2 + i));
}

double xi_nil(ui n2)
{
    return std::sqrt(2 * n2 + 1);
}

auto frac(double f, double s) {
    return 1 / (f + s);
}

auto j(ui n2, ui beta1, ui beta3)
{
    return (frac(xi_i(n2, beta1), xi_nil(n2)) * frac(xi_i(n2, beta1), xi_i(n2, beta3)) * frac(xi_nil(n2), xi_i(n2, beta3))) \
        + 2 * std::pow(xi_i(n2, beta3), 2) * (xi_i(n2, beta1) * (2 * xi_nil(n2) + xi_i(n2, beta3)) + xi_nil(n2)*(xi_nil(n2) + 2*xi_i(n2, beta3))) * std::pow(frac(xi_i(n2, beta1), xi_nil(n2)), 2) \
        * std::pow(frac(xi_i(n2, beta1), xi_i(n2, beta3)), 3) * std::pow(frac(xi_nil(n2), xi_i(n2, beta3)), 2) \
        - xi_i(n2, beta3) * (xi_i(n2, beta1) + xi_nil(n2) + 2 * xi_i(n2, beta3)) * frac(xi_i(n2, beta1), xi_nil(n2)) * std::pow(frac(xi_i(n2, beta1), xi_i(n2, beta3)), 2) \
        * std::pow(frac(xi_nil(n2), xi_i(n2, beta3)), 2) \
        - ((xi_i(n2, beta1) * (xi_nil(n2) + xi_i(n2, beta3)) + 2 * xi_nil(n2) * xi_i(n2, beta3)) * std::pow(frac(xi_i(n2, beta1), xi_nil(n2)), 2) * std::pow(frac(xi_i(n2, beta1), xi_i(n2, beta3)), 2) * frac(xi_nil(n2), xi_i(n2, beta3)));
}

int main()
{
    std::vector<int> max_num = { 100, 1000, 10000, 100000,1000000 };

    for (int i = 0; i < max_num.size(); i++) {
        double total = j(0, 1, 1);

        for (ui n2 = 1; n2 <= max_num[i]; n2++) {
            for (ui beta1 = 0; beta1 <= 1; beta1++) {
                for (ui beta3 = 0; beta3 <= 1; beta3++) {
                    total += (n2 + 1.0) * std::pow(-1, beta1 + beta3) * j(n2, beta1, beta3);
                }
            }
        }
        std::cout << "Max number:" << max_num[i] << std::endl;
        std::cout.precision(15);
        std::cout << "Result:" << total << std::endl;
    }
}