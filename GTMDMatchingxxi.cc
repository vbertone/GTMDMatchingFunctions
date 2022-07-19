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
 * This code computes F11 as a function of x at fixed kT, t, and Q for
 * a vector of values in xi.
 */
int main(int argc, char** argv)
{
  // Parameters
  const double Qref  = 91.1876;
  const double asref = 0.118;
  const double mc    = 1.5;
  const double mb    = 4.75;
  const double mt    = 175;
  const double t     = -0.1;
  const double kT    = 1;
  const double mu    = 10;
  const std::vector<double> xiv{0.1, 0.3, 0.5, 0.7};

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb, mt};

  // Init PARTONS application
  PARTONS::Partons* pPartons = PARTONS::Partons::getInstance();
  pPartons->init(argc, argv);

  // Retrieve GPD service
  PARTONS::GPDService* pGPDService = PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

  // Create GPD module with the BaseModuleFactory
  PARTONS::GPDModule* pGPDModel = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(PARTONS::GPDGK16::classId);

  // b* prescription
  const std::function<double(double const&, double const&)> bs = NangaParbat::bstarMap.at("bstarmin");

  // Get PV19 parameterisation and set parameters
  NangaParbat::Parameterisation *NPFunc = NangaParbat::GetParametersation("DWS");
  //NangaParbat::Parameterisation *NPFunc = NangaParbat::GetParametersation("PV19x");
  //NPFunc->SetParameters({0.05763262396815123, 2.975948125354748, 1.150602082569575, 1.007523907835578, 0.4106886687573338, 0.02063454913013714, 0.09390144417638162, 0.02737272626857176, 0.02266728682675793})

  // Configure APFEL++
  apfel::SetVerbosityLevel(2);
  apfel::Banner();

  // Running coupling
  apfel::AlphaQCD a{asref, Qref, Thresholds, 1};
  const  apfel::TabulateObject<double> Alphas{a, 100, 0.9, 100, 3};
  const auto as = [=] (double const& mur) -> double{ return Alphas.Evaluate(mur); };

  // x-space grid
  const apfel::Grid g{{{100, 1e-7, 3}, {100, 1e-2, 3}, {100, 6e-1, 3}, {100, 8e-1, 3}}};

  // Vector of set of distributions
  std::vector<apfel::Set<apfel::Distribution>> Fv;

  // Run over values of xi
  for (double const& xi: xiv)
    {
      std::cout << "xi = " << xi << std::endl;
      // Input GPDs from PARTONS
      const auto InGPDs = [=] (double const& x, double const& Q) -> std::map<int, double>
	{
	  const double xi2 = xi * xi;
	  const PARTONS::GPDResult gpdResult = pGPDService->computeSingleKinematic(PARTONS::GPDKinematic{x, xi, t, Q * Q, Q * Q}, pGPDModel);
	  std::map<int, double> PhysMap{{-6, 0}, {-5, 0}, {-4, 0}, {-3, 0}, {-2, 0}, {-1, 0}, {0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}};
	  PhysMap[3]  = ( 1 - xi2 ) * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE).getQuarkDistribution()
	  - xi2 * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE).getQuarkDistribution();
	  PhysMap[2]  = ( 1 - xi2 ) * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistribution()
	  - xi2 * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistribution();
	  PhysMap[1]  = ( 1 - xi2 ) * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistribution()
	  - xi2 * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistribution();
	  PhysMap[0]  = ( 1 - xi2 ) * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getGluonDistribution().getGluonDistribution()
	  - xi2 * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getGluonDistribution().getGluonDistribution();
	  PhysMap[-1] = ( 1 - xi2 ) * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistributionPlus()
	  - xi2 * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistributionPlus() - PhysMap[1];
	  PhysMap[-2] = ( 1 - xi2 ) * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistributionPlus()
	  - xi2 * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistributionPlus() - PhysMap[2];
	  PhysMap[-3] = ( 1 - xi2 ) * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE).getQuarkDistributionPlus()
	  - xi2 * x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::E).getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE).getQuarkDistributionPlus() - PhysMap[3];
	  return apfel::PhysToQCDEv(PhysMap);
	};

      // Evolve and tabulate GPDs
      apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs{*(BuildDglap(InitializeGpdObjects(g, Thresholds, xi), InGPDs, sqrt(pGPDModel->getMuF2Ref()), 0, as)), 100, 1, 1000, 3};
      const auto CollGPDs = [=] (double const& muf) -> apfel::Set<apfel::Distribution> { return TabulatedGPDs.Evaluate(muf); };

      // Now compute GTMDs in bT space
      const auto EvGTMDs = BuildGtmds(apfel::InitializeGtmdObjects(g, Thresholds, xi), CollGPDs, as, 2);

      // Double exponential quadrature
      apfel::DoubleExponentialQuadrature DEObj{};

      // Define bT-dependent function returning set of distribution
      const std::function<apfel::Set<apfel::Distribution>(double const&)> xGb = [=] (double const& bT) -> apfel::Set<apfel::Distribution>
	{
	  apfel::Set<apfel::Distribution> G = EvGTMDs(bs(bT, mu), mu, mu * mu);
	  G.SetMap(apfel::ConvolutionMap{"Void"});
	  return [=] (double const& x) -> double { return bT * NPFunc->Evaluate(x, bT, mu * mu, 0); } * G;
	};

      // Compute Fourier tranform
      Fv.push_back(DEObj.transform(xGb, kT));
    }

  // Tabulation parameters
  const int nx = 10000;
  const double xmin = 0.001;
  const double xmax = 0.9;
  const double xstp = exp( log( xmax / xmin ) / ( nx - 1 ) );
  std::cout << std::scientific;
  std::cout << "#    x            GTMD ..." << std::endl;
  for (double x = xmin; x <= xmax * ( 1 + 1e-5 ); x *= xstp)
    {
      std::cout << x << "  ";
      for (auto const& F : Fv)
	std::cout << QCDEvToPhys(F.GetObjects()).at(2).Evaluate(x) << "  ";
      std::cout << std::endl;
    }

  std::cout << "\n#    x            GTMD ..." << std::endl;
  for (double x = xmin; x <= xmax * ( 1 + 1e-5 ); x *= xstp)
    {
      std::cout << x << "  ";
      for (auto const& F : Fv)
	std::cout << QCDEvToPhys(F.GetObjects()).at(0).Evaluate(x) << "  ";
      std::cout << std::endl;
    }

  // Remove pointer references
  PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(pGPDModel, 0);
  pGPDModel = 0;

  // Close PARTONS application properly
  pPartons->close();

  return 0;
}
