// NangaParbat
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <apfel/apfelxx.h>
/*
 * This code computes S11a0+q as a function of kT at fixed x, t, and Q for
 * a vector of values in xi.
 * The computation ignores (for the moment) the polarization mixing in the matching and
 * focuses only on the diagonal UU channel.
 * The direction of the staple-like Wilson line is assumed to be +\infty, i.e. SIDIS-like (=> s=+1)
 */
int main(int argc, char **argv)
{

  // Bind output stream to either stdout or file stream
  std::ofstream file_stream;
  std::ostream *stream_ptr = nullptr;

  if (argc >= 2 && std::string(argv[1]) == "file")
    {
      file_stream.open("output_CSkernel.dat");
      stream_ptr = &file_stream;
    }
  else
    {
      stream_ptr = &std::cout;
    }

  std::ostream &output_stream = *stream_ptr;

  // Parameters
  const double Qref = 91.1876;
  const double asref = 0.118;
  const double mc = 1.5;
  const double mb = 4.75;
  const double mt = 175;
  const double t = -0.1;
  const double xb = 0.2;
  const double mu = 10;
  const int ifl = 2;
  // const std::vector<double> xiv{0.001, 0.1, 0.3, 0.5, 0.7};
  const std::vector<double> xiv{0.3};

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb, mt};

  // Configure APFEL++
  apfel::SetVerbosityLevel(2);
  apfel::Banner();

  // Running coupling
  apfel::AlphaQCD a{asref, Qref, Thresholds, 1};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 100, 3};
  const auto as = [=](double const &mu) -> double { return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{{100, 1e-7, 3}, {100, 1e-1, 3}, {100, 6e-1, 3}, {100, 8e-1, 3}}};

  // Now compute TMDs in bT space
  // Cache the tmdObj to use it to get CS kernel
  const auto tmdObj = apfel::InitializeTmdObjects(g, Thresholds);
  auto CSkernel = apfel::CollinsSoperKernel(tmdObj, as, 2, 1, 1.0e-7);

  // Tabulation parameters
  const int nbT = 100;
  const double bTmin = 1e-4;
  const double bTmax = 2;
  const double bTstp = (bTmax - bTmin) / (nbT - 1);
  output_stream << std::scientific;
  output_stream << "# bT [GeV^-1]         GTMD (T-even  T-odd) at xi = [" << xiv[0];
  for (size_t i = 1; i < xiv.size(); i++)
    output_stream << ", " << xiv[i];
  output_stream << "]" << std::endl;

  for (double bT = bTmin; bT <= bTmax * (1 + 1e-5); bT += bTstp)
    {
      output_stream << bT << "  ";
      output_stream << M_PI * 0.5 * CSkernel(bT, mu) << "  ";
      output_stream << std::endl;
    }

  return 0;
}
