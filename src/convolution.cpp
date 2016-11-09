#include "convolution.h"

std::vector<double> input_signal()
{
  std::vector<double> result(101, 0);

  for(int t = 0; t < 101; ++t) {
    result[t] = skew_gaussian((t - 50)/10.0, 4);
  }

  return result;
}

std::vector<double> output_signal()
{
  auto input(input_signal());
  std::vector<double> result(101, 0);

  const double delay = 8.4;
  auto parts(factor_delay(delay));

  UniformLagrangeSet uls(parts.second);

  for(int t = parts.first + config.interpolation_order; t < 101; ++t) {

    for(int i = 0; i <= config.interpolation_order; ++i) {
      result[t] += uls.weights[1][i]*input[t - parts.first - i];
    }
  }

  return result;
}

std::pair<int, double> factor_delay(const double delay)
{
  std::pair<int, double> parts(std::floor(delay), delay - std::floor(delay));
  return parts;
}
