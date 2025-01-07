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
      file_stream.open("output_Main.dat");
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

  // Init PARTONS application
  PARTONS::Partons *pPartons = PARTONS::Partons::getInstance();
  pPartons->init(argc, argv);

  // Retrieve GPD service
  PARTONS::GPDService *pGPDService = PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

  // Create GPD module with the BaseModuleFactory
  PARTONS::GPDModule *pGPDModel = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(PARTONS::GPDGK16::classId);

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

  // Input PDFs (forward limit of GPDs)
  const auto InPDFs = [=](double const &x, double const &Q) -> std::map<int, double>
  {
    PARTONS::GPDResult gpdResult = pGPDService->computeSingleKinematic(PARTONS::GPDKinematic{x, 0.0001, t, Q * Q, Q * Q}, pGPDModel);
    std::map<int, double> PhysMap{{-6, 0}, {-5, 0}, {-4, 0}, {-3, 0}, {-2, 0}, {-1, 0}, {0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}};
    PhysMap[3] =
        x *
        gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE).getQuarkDistribution();
    PhysMap[2] =
        x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistribution();
    PhysMap[1] =
        x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistribution();
    PhysMap[0] = gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getGluonDistribution().getGluonDistribution();
    PhysMap[-1] =
        x *
        gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistributionPlus();
    PhysMap[-2] =
        x *
        gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistributionPlus();
    PhysMap[-3] = x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H)
                          .getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE)
                          .getQuarkDistributionPlus();
    return apfel::PhysToQCDEv(PhysMap);
  };

  // Evolve and tabulate PDFs
  apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{
      *(BuildDglap(InitializeDglapObjectsQCD(g, Thresholds), InPDFs, sqrt(pGPDModel->getMuF2Ref()), 0, as)), 100, 1, 1000, 3};
  const auto CollPDFs = [=](double const &mu) -> apfel::Set<apfel::Distribution> { return TabulatedPDFs.Evaluate(mu); };

  // Now compute TMDs in bT space
  // Cache the tmdObj to use it to get CS kernel
  const auto tmdObj = apfel::InitializeTmdObjects(g, Thresholds);
  auto CSkernel = apfel::CollinsSoperKernel(tmdObj, as, 2);
  const auto EvTMDs = BuildTmdPDFs(apfel::InitializeTmdObjects(g, Thresholds), CollPDFs, as, 2);
  const std::function<double(double const &)> xTb = [=](double const &bT) -> double
  { return bT * QCDEvToPhys(EvTMDs(bs(bT, mu), mu, mu * mu).GetObjects()).at(ifl).Evaluate(xb) * NPFunc->Evaluate(xb, bT, mu * mu, 0); };
  const apfel::TabulateObject<double> TabxTb{
      xTb, 100, 0.00005, 100, 5, {}, [](double const &x) -> double { return log(x); }, [](double const &x) -> double { return exp(x); }};
  const std::function<double(double const &)> txTb = [=](double const &bT) -> double { return TabxTb.Evaluate(bT); };

  // Run over values of xi
  // T-even part
  std::vector<std::function<double(double const &)>> txGb;
  // T-odd part
  std::vector<std::function<double(double const &)>> txGb_odd;
  for (double const &xi : xiv)
    {
      std::cout << "xi = " << xi << std::endl;
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
        PhysMap[2] = (1 - xi2) * x *
                         gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H)
                             .getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP)
                             .getQuarkDistribution() -
                     xi2 * x *
                         gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E)
                             .getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP)
                             .getQuarkDistribution();
        PhysMap[1] = (1 - xi2) * x *
                         gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H)
                             .getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN)
                             .getQuarkDistribution() -
                     xi2 * x *
                         gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E)
                             .getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN)
                             .getQuarkDistribution();
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

      // Evolve and tabulate GPDs
      apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs{
          *(BuildDglap(InitializeGpdObjects(g, Thresholds, xi), InGPDs, sqrt(pGPDModel->getMuF2Ref()), 0, as)), 100, 1, 1000, 3};
      const auto CollGPDs = [=](double const &mu) -> apfel::Set<apfel::Distribution> { return TabulatedGPDs.Evaluate(mu); };

      // phi mixing angle
      const auto mixing_angle = [=](double const &bT) -> double { return xb > xi ? 0 : M_PI * 0.5 * CSkernel(bT, mu); };

      // Now compute GTMDs in bT space
      // Tabulate directly with the mixing angle
      std::cout << "Start of computation of GTMD in bT space..." << std::endl;
      const auto EvGTMDs = BuildGtmds(apfel::InitializeGtmdObjectsEvenUU(g, Thresholds, xi), CollGPDs, as, 2);
      const auto EvGTMDs_odd = BuildGtmds(apfel::InitializeGtmdObjectsOddUU(g, Thresholds, xi), CollGPDs, as, 2);

      const std::function<double(double const &)> xGb = [=](double const &bT) -> double
      {
        double l1 = cos(mixing_angle(bT));
        double l2 = sin(mixing_angle(bT));
        double even = QCDEvToPhys(EvGTMDs(bs(bT, mu), mu, mu * mu).GetObjects()).at(ifl).Evaluate(xb);
        double odd = QCDEvToPhys(EvGTMDs_odd(bs(bT, mu), mu, mu * mu).GetObjects()).at(ifl).Evaluate(xb);
        return bT * NPFunc->Evaluate(xb, bT, mu * mu, 0) * (l1 * even - l2 * odd);
      };
      const apfel::TabulateObject<double> TabxGb{
          xGb, 100, 0.00005, 100, 5, {}, [](double const &x) -> double { return log(x); }, [](double const &x) -> double { return exp(x); }};

      const std::function<double(double const &)> xGb_odd = [=](double const &bT) -> double
      {
        double l1 = cos(mixing_angle(bT));
        double l2 = sin(mixing_angle(bT));
        double even = QCDEvToPhys(EvGTMDs(bs(bT, mu), mu, mu * mu).GetObjects()).at(ifl).Evaluate(xb);
        double odd = QCDEvToPhys(EvGTMDs_odd(bs(bT, mu), mu, mu * mu).GetObjects()).at(ifl).Evaluate(xb);
        return bT * NPFunc->Evaluate(xb, bT, mu * mu, 0) * (l1 * odd + l1 * even);
      };
      const apfel::TabulateObject<double> TabxGb_odd{
          xGb_odd, 100, 0.00005, 100, 5, {}, [](double const &x) -> double { return log(x); }, [](double const &x) -> double { return exp(x); }};
      std::cout << "...done" << std::endl;

      // Now store the results
      // This makes apfel throw and Intrgrator error:
      // [apfel::Integrator::integrate] Error: Too high accuracy required.
      // when I later call the DoubleExponentialQuadrature if xi>xb
      // it is not the T-odd part that causes the issue, but the mixing angle

      txGb.push_back([=](double const &bT) -> double { return TabxGb.Evaluate(bT); });
      txGb_odd.push_back([=](double const &bT) -> double { return TabxGb_odd.Evaluate(bT); });
    }

  bool transform_in_qT = false;
  if (transform_in_qT)
    {
      // Double exponential quadrature
      apfel::DoubleExponentialQuadrature DEObj{};

      // Tabulation parameters
      const int nqT = 100;
      const double qTmin = 1e-4;
      const double qTmax = 2;
      const double qTstp = (qTmax - qTmin) / (nqT - 1);
      output_stream << std::scientific;
      output_stream << "# kT [GeV]         GTMD (T-even  T-odd) at xi = [" << xiv[0];
      for (size_t i = 1; i < xiv.size(); i++)
        output_stream << ", " << xiv[i];
      output_stream << "]" << std::endl;

      for (double qT = qTmin; qT <= qTmax * (1 + 1e-5); qT += qTstp)
        {
          // output_stream << qT << "  " << DEObj.transform(txTb, qT) << "  ";
          output_stream << qT << "  ";
          for (size_t i = 0; i < txGb.size(); i++)
            {
              const auto &f_even = txGb[i];
              const auto &f_odd = txGb_odd[i];
              output_stream << DEObj.transform(f_even, qT) << "  ";
              output_stream << DEObj.transform(f_odd, qT) << "  ";
            }
          output_stream << std::endl;
        }
    }
  else
    {
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
          for (size_t i = 0; i < txGb.size(); i++)
            {
              const auto &f_even = txGb[i];
              const auto &f_odd = txGb_odd[i];
              output_stream << f_even(bT) << "  ";
              output_stream << f_odd(bT) << "  ";
            }
          output_stream << std::endl;
        }
    }

  // Remove pointer references
  PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(pGPDModel, 0);
  pGPDModel = 0;

  // Close PARTONS application properly
  pPartons->close();

  return 0;
}
