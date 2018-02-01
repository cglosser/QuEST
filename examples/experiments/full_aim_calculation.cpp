#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "interactions/AIM/aim_interaction.h"
#include "interactions/AIM/nearfield_interaction.h"
#include "interactions/direct_interaction.h"
#include "math_utils.h"

using dbl = std::numeric_limits<double>;

std::vector<Eigen::Vector3d> pos = {
    {0.0954224, 0.0609634, 0.0399523},    {0.0686174, 0.00415103, 0.0582277},
    {0.0195859, 0.0728217, 0.0936565},    {0.0569936, 0.0894418, 0.0167345},
    {0.0971071, 0.0748458, 0.0381879},    {0.0234278, 0.0626966, 0.0936808},
    {0.0314575, 0.08227, 0.0568064},      {0.0298478, 0.0389792, 0.0821379},
    {0.00979624, 0.0586154, 0.0846908},   {0.0488768, 0.0413799, 0.0667558},
    {0.0602131, 0.0518489, 0.0498444},    {0.0331995, 0.0919088, 0.0863656},
    {0.0615184, 0.0300746, 0.0850671},    {0.0779737, 0.0672925, 0.0730399},
    {0.0319344, 0.0812367, 0.00040785},   {0.0178925, 0.0968027, 0.0587978},
    {0.050742, 0.0572987, 0.00303841},    {0.0586277, 0.0935189, 0.0285718},
    {0.0498727, 0.0108634, 0.0357905},    {0.0112937, 0.0234473, 0.0895191},
    {0.0776108, 0.0400216, 0.00125684},   {0.0344503, 0.0281689, 0.0704969},
    {0.0236081, 0.0723685, 0.0935699},    {0.00173399, 0.0353924, 0.0405633},
    {0.0220393, 0.0804712, 0.0380298},    {0.0904788, 0.0399545, 0.0749272},
    {0.0521299, 0.0759129, 0.063545},     {0.0789527, 0.0566495, 0.0776189},
    {0.0625812, 0.0372678, 0.0452288},    {0.0434298, 0.0319982, 0.0104551},
    {0.0491623, 0.0612406, 0.0232817},    {0.0193167, 0.00691747, 0.0797625},
    {0.0948098, 0.0279116, 0.0113815},    {0.0417427, 0.0383013, 0.010508},
    {0.0603137, 0.0428067, 0.028353},     {0.0874788, 0.0353374, 0.0165139},
    {0.0785495, 0.098695, 0.0313298},     {0.00960504, 0.0134366, 0.0955654},
    {0.087883, 0.00781862, 0.0495329},    {0.0602336, 0.090586, 0.0117231},
    {0.0808116, 0.0220835, 0.00957732},   {0.0665134, 0.0093696, 0.0195616},
    {0.0325094, 0.0460713, 0.0917219},    {0.0394414, 0.0420254, 0.0506436},
    {0.048883, 0.0514117, 0.0820448},     {0.0523288, 0.0151846, 0.0756223},
    {0.0719197, 0.051445, 0.0674087},     {0.0239917, 0.0702328, 0.0666112},
    {0.0966959, 0.0363485, 0.0790092},    {0.0769881, 0.00759115, 0.0487968},
    {0.0295569, 0.0200339, 0.0306084},    {0.0810226, 0.0137323, 0.042833},
    {0.000840529, 0.00473207, 0.0754876}, {0.00358216, 0.097329, 0.0754842},
    {0.081163, 0.0387827, 0.0384688},     {0.0130163, 0.0543467, 0.0766495},
    {0.0135389, 0.0428736, 0.0469198},    {0.0360873, 0.0389936, 0.00587089},
    {0.0904368, 0.0259487, 0.0161238},    {0.0695068, 0.0103373, 0.0983447},
    {0.0491897, 0.067746, 0.0469562},     {0.0351955, 0.0383759, 0.05199},
    {0.00789417, 0.0182181, 0.078229},    {0.00417909, 0.0257417, 0.0804767},
    {0.0876679, 0.011037, 0.0976327},     {0.0162453, 0.000480723, 0.0946131},
    {0.0116662, 0.0538591, 0.0913515},    {0.00753958, 0.00648246, 0.0847026},
    {0.0419235, 0.0213808, 0.0181541},    {0.0263381, 0.0317762, 0.0853226},
    {0.0728322, 0.0398005, 0.0661913},    {0.0658767, 0.0676771, 0.0918251},
    {0.0681212, 0.0438428, 0.0809407},    {0.0273335, 0.0603806, 0.0700062},
    {0.0832255, 0.0466491, 0.0380107},    {0.0565619, 0.0880947, 0.0557899},
    {0.0784592, 0.0704105, 0.0307072},    {0.00428413, 0.00575125, 0.0392566},
    {0.0964955, 0.0889505, 0.0591572},    {0.0861959, 0.0477746, 0.0102702},
    {0.0933026, 0.0380676, 0.0479337},    {0.0770583, 0.0311847, 0.0216933},
    {0.0968405, 0.0731884, 0.0245722},    {0.0396713, 0.0225021, 0.0646154},
    {0.0646238, 0.0658822, 0.0255969},    {0.0963668, 0.077285, 0.0616678},
    {0.0331621, 0.0472945, 0.0558715},    {0.025339, 0.0162845, 0.0315776},
    {0.0584159, 0.0225438, 0.0465912},    {0.0241927, 0.00632007, 0.0857199},
    {0.0602682, 0.0697253, 0.0568655},    {0.00310323, 0.0261775, 0.0199437},
    {0.0852613, 0.0444332, 0.0149508},    {0.0580737, 0.0575892, 0.0710267},
    {0.0120281, 0.092466, 0.0402283},     {0.0967731, 0.024107, 0.00221449},
    {0.0853123, 0.0523127, 0.0814275},    {0.0758226, 0.0764629, 0.0528294},
    {0.0337446, 0.0384823, 0.0797621},    {0.0168601, 0.0497013, 0.0558615},
    {0.000760188, 0.00581333, 0.0480155}, {0.0322774, 0.00699445, 0.000347651},
    {0.0448472, 0.0689977, 0.0578067},    {0.0860148, 0.00366898, 0.0472598},
    {0.00811128, 0.0687212, 0.068439},    {0.0200804, 0.0213819, 0.0339454},
    {0.0417784, 0.0511735, 0.0308977},    {0.0984962, 0.0199995, 0.0781608},
    {0.0335732, 0.0974443, 0.0570362},    {0.0318072, 0.0101994, 0.00515029},
    {0.0992792, 0.0920189, 0.0132683},    {0.0825704, 0.0766344, 0.0624423},
    {0.0445048, 0.0428458, 0.073645},     {0.0854136, 0.0224455, 0.0486997},
    {0.0829876, 0.0953205, 0.0541524},    {0.0221395, 0.0795186, 0.0191285},
    {0.0064079, 0.0734483, 0.0854839},    {0.0375544, 0.0795873, 0.0684849},
    {0.0553026, 0.0673325, 0.016394},     {0.0376061, 0.0675043, 0.0464177},
    {0.0777504, 0.0639373, 0.0382685},    {0.00868848, 0.07245, 0.0637925},
    {0.0039547, 0.0874557, 0.0830712},    {0.0097508, 0.0587914, 0.0962246},
    {0.0271717, 0.00197449, 0.041311},    {0.0854734, 0.0129924, 0.0649677},
    {0.0734039, 0.0425005, 0.0138972},    {0.0379465, 0.0820081, 0.00732839}};

struct PARAMETERS {
  int interpolation_order, expansion_order, num_steps, num_dots;
  double c, mu0, dt, total_time, omega;

  Eigen::Array3i num_boxes;
  Eigen::Array3d spacing;

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  AIM::Grid grid;
  AIM::Expansions::ExpansionTable expansions;

  PARAMETERS()
      : interpolation_order(4),
        expansion_order(4),
        num_steps(1024),
        num_dots(pos.size()),

        c(299.792458),
        mu0(2.0133545e-4),  // meV ps^2/e^2 um
        dt(0.25e-2),
        total_time(num_steps * dt),
        omega(2278.9013),

        num_boxes(Eigen::Vector3i(8, 8, 8)),
        spacing(Eigen::Array3d(1, 1, 1) * c * (M_PI / (20 * omega))),

        dots(std::make_shared<DotVector>()),
        history(std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            num_dots, 10, num_steps)),

        grid(spacing, dots, expansion_order),
        expansions(AIM::Expansions::LeastSquaresExpansionSolver::get_expansions(
            expansion_order, grid, *dots))
  {
    history->fill(Eigen::Vector2cd::Zero());
    for(int d = 0; d < num_dots; ++d) {
      for(int i = -10; i < num_steps; ++i) {
        history->array[d][i][0](RHO_01) = gaussian_source(i * dt);
      }
    }
  };

  double gaussian_source(const double t) const
  {
    double arg = (t - total_time / 2.0) / (total_time / 12.0);
    return gaussian(arg);
  }
};

int main()
{
  PARAMETERS params;
  for(const auto &r : pos) {
    params.dots->push_back(
        QuantumDot(r, params.omega, {10, 20}, {0, 0, 5.29177e-4}));
  }
  std::cout << "Spacing: " << params.spacing.transpose() << std::endl;

  // == DIRECT STUFF ==========================================================

  Propagation::RotatingFramePropagator gf(-params.mu0 / (4 * M_PI), params.c,
                                          params.omega);

  DirectInteraction direct(params.dots, params.history, gf,
                           params.interpolation_order, params.c, params.dt);

  std::ofstream direct_file("direct_field.dat");
  direct_file.precision(dbl::max_digits10);
  for(int i = 0; i < params.num_steps; ++i) {
    direct_file << i * params.dt << " " << direct.evaluate(i).transpose()
                << std::endl;
  }

  // == AIM STUFF =============================================================

  AIM::Grid grid(params.spacing, params.dots, params.expansion_order);
  auto expansion_table =
      AIM::Expansions::LeastSquaresExpansionSolver::get_expansions(
          params.expansion_order, grid, *params.dots);

  AIM::NearfieldInteraction nf(params.dots, params.history, gf,
                               params.interpolation_order, params.c, params.dt,
                               grid);

  AIM::AimInteraction aim(
      params.dots, params.history, params.interpolation_order, params.c,
      params.dt, grid, expansion_table,
      AIM::Expansions::RotatingEFIE(
          grid.max_transit_steps(params.c, params.dt) +
              params.interpolation_order,
          params.c, params.dt, params.omega),
      AIM::normalization::Helmholtz(params.omega / params.c,
                                    -params.mu0 / (4 * M_PI)));

  std::ofstream aim_file("aim_field.dat");
  aim_file.precision(dbl::max_digits10);
  aim_file << std::setw(dbl::max_digits10);
  for(int i = 0; i < params.num_steps; ++i) {
    aim_file << i * params.dt << " "
             << (nf.evaluate(i) + aim.evaluate(i)).transpose() << std::endl;
  }

  return 0;
}
