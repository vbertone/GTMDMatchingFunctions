// PARTONS
#include <ElementaryUtils/logger/CustomException.h>
#include <ElementaryUtils/logger/LoggerManager.h>
#include <partons/Partons.h>
#include <partons/services/automation/AutomationService.h>
#include <partons/ServiceObjectRegistry.h>
#include <partons/services/GPDService.h>
#include <partons/ModuleObjectFactory.h>
#include <partons/modules/gpd/GPDGK16.h>

// NangaParbat
#include <NangaParbat/bstar.h>
#include <NangaParbat/nonpertfunctions.h>

#include <thread>
#include <mutex>

// Global mutex for synchronization
std::mutex mutex_to_print;

void S11a0plusq(std::string filename, int argc, char **argv);
void S11b0plusq(std::string filename, int argc, char **argv);

// Using M^2 = 1 GeV^2:
// For |t| = 4/3, the maximal value of \xi is 1/2
// For |t| =   1, the maximal value of \xi is 1/sqrt(5) ~ 0.4472
// For |t| = 1/2, the maximal value of \xi is 1/3

std::map<double, double> get_t_vector(std::vector<double> const &xiv, const double &t_offset)
{
  (void)t_offset;
  std::map<double, double> res;
  // Assuming that M^2 = 1 GeV^2
  for (const double &xi : xiv)
    res[xi] = -0.5;
  // res[xi] = -std::abs(t_offset) - 4 * pow(xi, 2) / (1 - pow(xi, 2));
  return res;
}

double get_t(double const &xi, const double &t_offset) { return -std::abs(t_offset) - 4 * pow(xi, 2) / (1 - pow(xi, 2)); }

int main(int argc, char **argv)
{

  // S11a0plusq("S11a0plusUp.dat", argc, argv);
  // S11b0plusq("S11b0plusUp.dat", argc, argv);
  S11a0plusq("S11a0plusUp_x_sim_xi.dat", argc, argv);

  return 0;
}

/*
 * This code computes S11a0+q as a function of kT at fixed x, t, and Q for
 * a vector of values in xi.
 * The computation ignores (for the moment) the polarization mixing in the matching and
 * focuses only on the diagonal UU channel.
 * The direction of the staple-like Wilson line is assumed to be +\infty, i.e. SIDIS-like (=> s=+1)
 */
//  // Retrieve GPD service
// PARTONS::GPDService *pGPDService = PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

// // Create GPD module with the BaseModuleFactory
// PARTONS::GPDModule *pGPDModel = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(PARTONS::GPDGK16::classId);

std::vector<std::function<double(double const &)>> S11a0plusq_single_xi(double const &xi, double const &t, PARTONS::GPDService *pGPDService,
                                                                        PARTONS::GPDModule *pGPDModel)
{
  // Parameters
  const double Qref = 91.1876;
  const double asref = 0.118;
  const double mc = 1.5;
  const double mb = 4.75;
  const double mt = 175;
  const double xb = 0.15;
  const double mu = 10;
  const int ifl = 2;

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb, mt};

  // b* prescription
  const std::function<double(double const &, double const &)> bs = NangaParbat::bstarMap.at("bstarmin");

  // Get PV19 parameterisation and set parameters
  NangaParbat::Parameterisation *NPFunc = NangaParbat::GetParametersation("PV19x");
  NPFunc->SetParameters({0.05763262396815123, 2.975948125354748, 1.150602082569575, 1.007523907835578, 0.4106886687573338, 0.02063454913013714,
                         0.09390144417638162, 0.02737272626857176, 0.02266728682675793});

  // Configure APFEL++
  apfel::SetVerbosityLevel(2);
  apfel::Banner();

  // Running coupling
  apfel::AlphaQCD a{asref, Qref, Thresholds, 1};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 100, 3};
  const auto as = [=](double const &mu) -> double { return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{{100, 1e-7, 3}, {100, 1e-1, 3}, {100, 6e-1, 3}, {100, 8e-1, 3}}};

  // Compute the CSkernel
  auto CSkernel = apfel::CollinsSoperKernel(apfel::InitializeTmdObjects(g, Thresholds), as, 2, 1, 1.0e-4);

  // Run over values of xi
  // T-even part UU, T-odd part UU, T-even part UT, T-odd part UT
  std::vector<std::function<double(double const &)>> txGb(4);

  // Input GPDs from PARTONS

  const auto InGPDs = [=](double const &x, double const &Q) -> std::map<int, double>
  {
    const double xi2 = xi * xi;
    const PARTONS::GPDResult gpdResult = pGPDService->computeSingleKinematic(PARTONS::GPDKinematic{x, xi, t, Q * Q, Q * Q}, pGPDModel);
    std::map<int, double> PhysMap{{-6, 0}, {-5, 0}, {-4, 0}, {-3, 0}, {-2, 0}, {-1, 0}, {0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}};
    PhysMap[3] = (1 - xi2) * x *
                     gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H)
                         .getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE)
                         .getQuarkDistribution() -
                 xi2 * x *
                     gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E)
                         .getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE)
                         .getQuarkDistribution();
    PhysMap[2] =
        (1 - xi2) * x *
            gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistribution() -
        xi2 * x *
            gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistribution();
    PhysMap[1] =
        (1 - xi2) * x *
            gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistribution() -
        xi2 * x *
            gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistribution();
    PhysMap[0] = (1 - xi2) * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getGluonDistribution().getGluonDistribution() -
                 xi2 * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getGluonDistribution().getGluonDistribution();
    PhysMap[-1] = (1 - xi2) * x *
                      gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H)
                          .getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN)
                          .getQuarkDistributionPlus() -
                  xi2 * x *
                      gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E)
                          .getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN)
                          .getQuarkDistributionPlus() -
                  PhysMap[1];
    PhysMap[-2] = (1 - xi2) * x *
                      gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H)
                          .getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP)
                          .getQuarkDistributionPlus() -
                  xi2 * x *
                      gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E)
                          .getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP)
                          .getQuarkDistributionPlus() -
                  PhysMap[2];
    PhysMap[-3] = (1 - xi2) * x *
                      gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H)
                          .getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE)
                          .getQuarkDistributionPlus() -
                  xi2 * x *
                      gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E)
                          .getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE)
                          .getQuarkDistributionPlus() -
                  PhysMap[3];
    return apfel::PhysToQCDEv(PhysMap);
  };

  const auto InGPDs_trans = [=](double const &x, double const &Q) -> std::map<int, double>
  {
    const double xi2 = xi * xi;
    const double Ratio = 0.1;
    const PARTONS::GPDResult gpdResult = pGPDService->computeSingleKinematic(PARTONS::GPDKinematic{x, xi, t, Q * Q, Q * Q}, pGPDModel);
    std::map<int, double> PhysMap{{-6, 0}, {-5, 0}, {-4, 0}, {-3, 0}, {-2, 0}, {-1, 0}, {0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}};
    PhysMap[3] = 0;
    PhysMap[2] = 0;
    PhysMap[1] = 0;
    // E_T - \xi \tilde{E}_T + 2 \tilde{H}_T \sim 0.1 * (E - \xi \tilde{E} + 2 \tilde{H})
    PhysMap[0] = Ratio * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getGluonDistribution().getGluonDistribution() -
                 Ratio * xi * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::Et).getGluonDistribution().getGluonDistribution() +
                 Ratio * 2 * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::Ht).getGluonDistribution().getGluonDistribution();
    PhysMap[-1] = 0;
    PhysMap[-2] = 0;
    PhysMap[-3] = 0;
    return apfel::PhysToQCDEv(PhysMap);
  };

  // Evolve and tabulate GPDs
  apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs{
      *(BuildDglap(InitializeGpdObjects(g, Thresholds, xi), InGPDs, sqrt(pGPDModel->getMuF2Ref()), 0, as)), 100, 1, 1000, 3};
  const auto CollGPDs = [=](double const &mu) -> apfel::Set<apfel::Distribution> { return TabulatedGPDs.Evaluate(mu); };

  // Evolve and tabulate GPDs
  apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs_trans{
      *(BuildDglap(InitializeGpdObjects(g, Thresholds, xi), InGPDs_trans, sqrt(pGPDModel->getMuF2Ref()), 0, as)), 100, 1, 1000, 3};
  const auto CollGPDs_trans = [=](double const &mu) -> apfel::Set<apfel::Distribution> { return TabulatedGPDs_trans.Evaluate(mu); };

  // phi mixing angle
  // adding non-perturbative part of CSkernel: - g2 b^2 / 4
  const double g2 = 0.12840; // Value from MAP22b2
  const auto NP_CSK = [&g2](double const &bT) -> double { return -g2 * bT * bT * 0.25; };
  const auto mixing_angle = [=](double const &bT) -> double { return xb > xi ? 0 : M_PI * 0.5 * (CSkernel(bs(bT, mu), mu) + NP_CSK(bT)); };

  // Now compute GTMDs in bT space
  // Tabulate directly with the mixing angle
  const auto EvGTMDs = BuildGtmds(apfel::InitializeGtmdObjectsEvenUU(g, Thresholds, xi), CollGPDs, as, 2);
  const auto EvGTMDs_odd = BuildGtmds(apfel::InitializeGtmdObjectsOddUU(g, Thresholds, xi), CollGPDs, as, 2);

  const auto EvGTMDs_trans = BuildGtmds(apfel::InitializeGtmdObjectsEvenUT(g, Thresholds, xi), CollGPDs_trans, as, 2);
  const auto EvGTMDs_trans_odd = BuildGtmds(apfel::InitializeGtmdObjectsOddUT(g, Thresholds, xi), CollGPDs_trans, as, 2);

  // First the UU channel
  const std::function<double(double const &)> xGb = [=](double const &bT) -> double
  {
    double even = QCDEvToPhys(EvGTMDs(bs(bT, mu), mu, mu * mu).GetObjects()).at(ifl).Evaluate(xb);
    return bT * NPFunc->Evaluate(xb, bT, mu * mu, 0) * (even);
  };
  const apfel::TabulateObject<double> TabxGb{
      xGb, 100, 0.00005, 100, 5, {}, [](double const &x) -> double { return log(x); }, [](double const &x) -> double { return exp(x); }};

  const std::function<double(double const &)> xGb_odd = [=](double const &bT) -> double
  {
    double odd = QCDEvToPhys(EvGTMDs_odd(bs(bT, mu), mu, mu * mu).GetObjects()).at(ifl).Evaluate(xb);
    return bT * NPFunc->Evaluate(xb, bT, mu * mu, 0) * (odd);
  };
  const apfel::TabulateObject<double> TabxGb_odd{
      xGb_odd, 100, 0.00005, 100, 5, {}, [](double const &x) -> double { return log(x); }, [](double const &x) -> double { return exp(x); }};

  // Then the UT channel
  const std::function<double(double const &)> xGb_trans = [=](double const &bT) -> double
  {
    double even = QCDEvToPhys(EvGTMDs_trans(bs(bT, mu), mu, mu * mu).GetObjects()).at(ifl).Evaluate(xb);
    return bT * NPFunc->Evaluate(xb, bT, mu * mu, 0) * (even);
  };
  const apfel::TabulateObject<double> TabxGb_trans{
      xGb_trans, 100, 0.00005, 100, 5, {}, [](double const &x) -> double { return log(x); }, [](double const &x) -> double { return exp(x); }};

  const std::function<double(double const &)> xGb_trans_odd = [=](double const &bT) -> double
  {
    double odd = QCDEvToPhys(EvGTMDs_trans_odd(bs(bT, mu), mu, mu * mu).GetObjects()).at(ifl).Evaluate(xb);
    return bT * NPFunc->Evaluate(xb, bT, mu * mu, 0) * (odd);
  };
  const apfel::TabulateObject<double> TabxGb_trans_odd{
      xGb_trans_odd, 100, 0.00005, 100, 5, {}, [](double const &x) -> double { return log(x); }, [](double const &x) -> double { return exp(x); }};

  // Now store the results as the appropriate linear combinations
  txGb[0] = [=](double const &bT) -> double
  {
    double l1 = cos(mixing_angle(bT));
    double l2 = sin(mixing_angle(bT));
    return TabxGb.Evaluate(bT) * l1 - TabxGb_odd.Evaluate(bT) * l2;
  };
  txGb[1] = [=](double const &bT) -> double
  {
    double l1 = cos(mixing_angle(bT));
    double l2 = sin(mixing_angle(bT));
    return TabxGb.Evaluate(bT) * l2 + TabxGb_odd.Evaluate(bT) * l1;
  };
  txGb[2] = [=](double const &bT) -> double
  {
    double l1 = cos(mixing_angle(bT));
    double l2 = sin(mixing_angle(bT));
    return TabxGb_trans.Evaluate(bT) * l1 - TabxGb_trans_odd.Evaluate(bT) * l2;
  };
  txGb[3] = [=](double const &bT) -> double
  {
    double l1 = cos(mixing_angle(bT));
    double l2 = sin(mixing_angle(bT));
    return TabxGb_trans.Evaluate(bT) * l2 + TabxGb_trans_odd.Evaluate(bT) * l1;
  };

  return txGb;
  //
  //
}

void S11a0plusq(std::string filename, int argc, char **argv)
{

  // Bind output stream to either stdout or file stream
  std::ofstream file_stream;
  std::ostream *stream_ptr = nullptr;

  file_stream.open(filename);
  stream_ptr = &file_stream;

  std::ostream &output_stream = *stream_ptr;

  // const std::vector<double> xiv = {0.001, 0.2, 0.3};
  const std::vector<double> xiv = {0.15001};
  const double t_offset = -0.1;
  const std::map<double, double> tv = get_t_vector(xiv, t_offset);

  // Init PARTONS application

  PARTONS::Partons *pPartons = PARTONS::Partons::getInstance();
  pPartons->init(argc, argv);

  // Retrieve GPD service
  PARTONS::GPDService *pGPDService = PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

  // Create GPD module with the BaseModuleFactory
  PARTONS::GPDModule *pGPDModel = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(PARTONS::GPDGK16::classId);

  // Results container: functions instead of computed values
  std::vector<std::function<double(double const &)>> res_UU_even(xiv.size());
  std::vector<std::function<double(double const &)>> res_UU_odd(xiv.size());
  std::vector<std::function<double(double const &)>> res_UT_even(xiv.size());
  std::vector<std::function<double(double const &)>> res_UT_odd(xiv.size());

  // Spawn threads to compute and store the callable functions
  for (size_t i = 0; i < xiv.size(); ++i)
    {
      auto tmp = S11a0plusq_single_xi(xiv[i], tv.at(xiv[i]), pGPDService, pGPDModel);
      res_UU_even[i] = tmp[0];
      res_UU_odd[i] = tmp[1];
      res_UT_even[i] = tmp[2];
      res_UT_odd[i] = tmp[3];
    }

  // Double exponential quadrature
  apfel::DoubleExponentialQuadrature DEObj{};
  apfel::DoubleExponentialQuadrature DEObj_trans{2};

  // Tabulation parameters
  const int nqT = 1000;
  const double qTmin = 1e-4;
  const double qTmax = 4;
  const double qTstp = (qTmax - qTmin) / (nqT - 1);
  output_stream << std::scientific;
  output_stream << "# kT [GeV]         GTMD (T-even-UU T-even-UT T-odd-UU T-odd-UT) at (xi, -t) = [(" << xiv[0] << ", " << tv.at(xiv[0]) << ")";
  for (size_t i = 1; i < xiv.size(); i++)
    output_stream << ", (" << xiv[i] << ", " << tv.at(xiv[i]) << ")";
  output_stream << "]" << std::endl;

  // List of prefactors, excluding the 2\pi Cos(2\theta_{k\Delta})
  // Assuming that M^2 = 1 GeV^2
  std::vector<double> prefactor;
  for (const double &xi : xiv)
    prefactor.emplace_back((1 - xi * xi) * (-tv.at(xi) / 4.0) - xi * xi);

  for (double qT = qTmin; qT <= qTmax * (1 + 1e-5); qT += qTstp)
    {
      output_stream << qT << "  ";
      for (size_t i = 0; i < xiv.size(); i++)
        {
          const auto &f_even = res_UU_even[i];
          const auto &f_odd = res_UU_odd[i];

          const auto &u_even = res_UT_even[i];
          const auto &u_odd = res_UT_odd[i];
          output_stream << DEObj.transform(f_even, qT) << " " << prefactor[i] * DEObj_trans.transform(u_even, qT) << "  ";
          output_stream << DEObj.transform(f_odd, qT) << " " << prefactor[i] * DEObj_trans.transform(u_odd, qT) << "  ";
        }
      output_stream << std::endl;
    }

  // Remove pointer references
  PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(pGPDModel, 0);
  pGPDModel = 0;

  // Close PARTONS application properly
  pPartons->close();
}

/*
 * This code computes S11b0+q as a function of kT at fixed x, t, and Q for
 * a vector of values in xi.
 * The computation ignores (for the moment) the polarization mixing in the matching and
 * focuses only on the diagonal UU channel.
 * The direction of the staple-like Wilson line is assumed to be +\infty, i.e. SIDIS-like (=> s=+1)
 */
std::vector<std::function<double(double const &)>> S11b0plusq_single_xi(double const &xi, double const &t, PARTONS::GPDService *pGPDService,
                                                                        PARTONS::GPDModule *pGPDModel)
{
  // Parameters
  const double Qref = 91.1876;
  const double asref = 0.118;
  const double mc = 1.5;
  const double mb = 4.75;
  const double mt = 175;
  const double xb = 0.15;
  const double mu = 10;
  const int ifl = 2;

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb, mt};

  // b* prescription
  const std::function<double(double const &, double const &)> bs = NangaParbat::bstarMap.at("bstarmin");

  // Get PV19 parameterisation and set parameters
  NangaParbat::Parameterisation *NPFunc = NangaParbat::GetParametersation("PV19x");
  NPFunc->SetParameters({0.05763262396815123, 2.975948125354748, 1.150602082569575, 1.007523907835578, 0.4106886687573338, 0.02063454913013714,
                         0.09390144417638162, 0.02737272626857176, 0.02266728682675793});

  // Configure APFEL++
  apfel::SetVerbosityLevel(2);
  apfel::Banner();

  // Running coupling
  apfel::AlphaQCD a{asref, Qref, Thresholds, 1};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 100, 3};
  const auto as = [=](double const &mu) -> double { return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{{100, 1e-7, 3}, {100, 1e-1, 3}, {100, 6e-1, 3}, {100, 8e-1, 3}}};

  // Compute the CSkernel
  auto CSkernel = apfel::CollinsSoperKernel(apfel::InitializeTmdObjects(g, Thresholds), as, 2, 1, 1.0e-4);

  // Run over values of xi
  // T-even part, T-odd part
  std::vector<std::function<double(double const &)>> txGb(2);

  // Input GPDs from PARTONS
  // Lack of models for transversely polarized gluon GPDs, assume that they are their 'non-transverse' counterpart
  // times some small conversion factor
  const auto InGPDs = [=](double const &x, double const &Q) -> std::map<int, double>
  {
    const double xi2 = xi * xi;
    const double Ratio = 0.1;
    const PARTONS::GPDResult gpdResult = pGPDService->computeSingleKinematic(PARTONS::GPDKinematic{x, xi, t, Q * Q, Q * Q}, pGPDModel);
    std::map<int, double> PhysMap{{-6, 0}, {-5, 0}, {-4, 0}, {-3, 0}, {-2, 0}, {-1, 0}, {0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}};
    PhysMap[3] = 0;
    PhysMap[2] = 0;
    PhysMap[1] = 0;
    PhysMap[0] = Ratio * xi * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getGluonDistribution().getGluonDistribution() -
                 Ratio * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::Et).getGluonDistribution().getGluonDistribution();
    PhysMap[-1] = 0;
    PhysMap[-2] = 0;
    PhysMap[-3] = 0;
    return apfel::PhysToQCDEv(PhysMap);
  };

  // Evolve and tabulate GPDs
  apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs{
      *(BuildDglap(InitializeGpdObjects(g, Thresholds, xi), InGPDs, sqrt(pGPDModel->getMuF2Ref()), 0, as)), 100, 1, 1000, 3};
  const auto CollGPDs = [=](double const &mu) -> apfel::Set<apfel::Distribution> { return TabulatedGPDs.Evaluate(mu); };

  // phi mixing angle
  // adding non-perturbative part of CSkernel: - g2 b^2 / 4
  const double g2 = 0.12840; // Value from MAP22b2
  const auto NP_CSK = [&g2](double const &bT) -> double { return -g2 * bT * bT * 0.25; };
  const auto mixing_angle = [=](double const &bT) -> double { return xb > xi ? 0 : M_PI * 0.5 * (CSkernel(bs(bT, mu), mu) + NP_CSK(bT)); };

  // Now compute GTMDs in bT space
  // Tabulate directly with the mixing angle
  std::cout << "Start of computation of GTMD in bT space..." << std::endl;
  const auto EvGTMDs = BuildGtmds(apfel::InitializeGtmdObjectsEvenUT(g, Thresholds, xi), CollGPDs, as, 2);
  const auto EvGTMDs_odd = BuildGtmds(apfel::InitializeGtmdObjectsOddUT(g, Thresholds, xi), CollGPDs, as, 2);

  const std::function<double(double const &)> xGb = [=](double const &bT) -> double
  {
    double even = QCDEvToPhys(EvGTMDs(bs(bT, mu), mu, mu * mu).GetObjects()).at(ifl).Evaluate(xb);
    return bT * NPFunc->Evaluate(xb, bT, mu * mu, 0) * (even);
  };
  const apfel::TabulateObject<double> TabxGb{
      xGb, 100, 0.00005, 100, 5, {}, [](double const &x) -> double { return log(x); }, [](double const &x) -> double { return exp(x); }};

  const std::function<double(double const &)> xGb_odd = [=](double const &bT) -> double
  {
    double odd = QCDEvToPhys(EvGTMDs_odd(bs(bT, mu), mu, mu * mu).GetObjects()).at(ifl).Evaluate(xb);
    return bT * NPFunc->Evaluate(xb, bT, mu * mu, 0) * (odd);
  };
  const apfel::TabulateObject<double> TabxGb_odd{
      xGb_odd, 100, 0.00005, 100, 5, {}, [](double const &x) -> double { return log(x); }, [](double const &x) -> double { return exp(x); }};
  std::cout << "...done" << std::endl;

  // Now store the results as the appropriate linear combinations
  txGb[0] = [=](double const &bT) -> double
  {
    double l1 = cos(mixing_angle(bT));
    double l2 = sin(mixing_angle(bT));
    return TabxGb.Evaluate(bT) * l1 - TabxGb_odd.Evaluate(bT) * l2;
  };
  txGb[1] = [=](double const &bT) -> double
  {
    double l1 = cos(mixing_angle(bT));
    double l2 = sin(mixing_angle(bT));
    return TabxGb.Evaluate(bT) * l2 + TabxGb_odd.Evaluate(bT) * l1;
  };

  return txGb;
}

void S11b0plusq(std::string filename, int argc, char **argv)
{

  // Bind output stream to either stdout or file stream
  std::ofstream file_stream;
  std::ostream *stream_ptr = nullptr;

  file_stream.open(filename);
  stream_ptr = &file_stream;

  std::ostream &output_stream = *stream_ptr;

  // Init PARTONS application
  PARTONS::Partons *pPartons = PARTONS::Partons::getInstance();
  pPartons->init(argc, argv);

  // Retrieve GPD service
  PARTONS::GPDService *pGPDService = PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

  // Create GPD module with the BaseModuleFactory
  PARTONS::GPDModule *pGPDModel = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(PARTONS::GPDGK16::classId);

  const std::vector<double> xiv = {0.001, 0.2, 0.3};
  // offset of t from its minimum value (t = t_min - |t_offset|)
  const double t_offset = -0.1;
  const std::map<double, double> tv = get_t_vector(xiv, t_offset);

  // Results container: functions instead of computed values
  std::vector<std::function<double(double const &)>> res_UT_even(xiv.size());
  std::vector<std::function<double(double const &)>> res_UT_odd(xiv.size());

  // Spawn threads to compute and store the callable functions
  for (size_t i = 0; i < xiv.size(); ++i)
    {
      auto tmp = S11b0plusq_single_xi(xiv[i], tv.at(xiv[i]), pGPDService, pGPDModel);
      res_UT_even[i] = tmp[0];
      res_UT_odd[i] = tmp[1];
    }

  // Double exponential quadrature
  apfel::DoubleExponentialQuadrature DEObj{2};

  // Tabulation parameters
  const int nqT = 1000;
  const double qTmin = 1e-4;
  const double qTmax = 5;
  const double qTstp = (qTmax - qTmin) / (nqT - 1);
  output_stream << std::scientific;
  output_stream << "# kT [GeV]         GTMD (T-even  T-odd) at (xi, -t) = [(" << xiv[0] << ", " << tv.at(xiv[0]) << ")";
  for (size_t i = 1; i < xiv.size(); i++)
    output_stream << ", (" << xiv[i] << ", " << tv.at(xiv[i]) << ")";
  output_stream << "]" << std::endl;

  // List of prefactors, excluding the -2\pi i Sin(2\theta_{k\Delta})
  // Assuming that M^2 = 1 GeV^2
  std::vector<double> prefactor;
  for (const double &xi : xiv)
    prefactor.emplace_back((1 - xi * xi) * (-tv.at(xi) / 4.0) - xi * xi);

  for (double qT = qTmin; qT <= qTmax * (1 + 1e-5); qT += qTstp)
    {
      output_stream << qT << "  ";
      for (size_t i = 0; i < xiv.size(); i++)
        {
          const auto &f_even = res_UT_even[i];
          const auto &f_odd = res_UT_odd[i];
          output_stream << prefactor[i] * DEObj.transform(f_even, qT) << "  ";
          output_stream << prefactor[i] * DEObj.transform(f_odd, qT) << "  ";
        }
      output_stream << std::endl;
    }

  // Remove pointer references
  PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(pGPDModel, 0);
  pGPDModel = 0;

  // Close PARTONS application properly
  pPartons->close();
}

/*
H = 2,        //!< Twist-2 GPD \f$H\f$
E = 3,        //!< Twist-2 GPD \f$E\f$
Ht = 4,       //!< Twist-2 GPD \f$\tilde{H}\f$
Et = 5,       //!< Twist-2 GPD \f$\tilde{E}\f$
HTrans = 6,   //!< Twist-2 GPD \f$H_{T}\f$
ETrans = 7,   //!< Twist-2 GPD \f$E_{T}\f$
HtTrans = 8,  //!< Twist-2 GPD \f$\tilde{H}_{T}\f$
EtTrans = 9,  //!< Twist-2 GPD \f$\tilde{E}_{T}\f$
*/