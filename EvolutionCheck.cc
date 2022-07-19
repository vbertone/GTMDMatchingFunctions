// PARTONS
#include <ElementaryUtils/logger/CustomException.h>
#include <ElementaryUtils/logger/LoggerManager.h>
#include <partons/Partons.h>
#include <partons/services/automation/AutomationService.h>
#include <partons/ServiceObjectRegistry.h>
#include <partons/services/GPDService.h>
#include <partons/ModuleObjectFactory.h>
#include <partons/modules/evolution/gpd/GPDEvolutionApfel.h>
#include <partons/modules/gpd/GPDGK16.h>
#include <partons/modules/running_alpha_strong/RunningAlphaStrongApfel.h>
#include <partons/modules/active_flavors_thresholds/ActiveFlavorsThresholdsVariable.h>

// NangaParbat
#include <NangaParbat/bstar.h>
#include <NangaParbat/nonpertfunctions.h>

/*
 * This code checkes that the evolution as provided by APFEL++ through
 * PARTONS coincides with that obtained by APFEL++ used as a
 * stand-alone code.
 */
int main(int argc, char** argv)
{
  // Evolution parameters
  const double Qref  = 91.1876;
  const double asref = 0.118;
  const double mc    = 1.5;
  const double mb    = 4.75;
  const double mt    = 175;
  const double t     = -0.1;
  const double xi    = 0.2;
  const double xi2   = xi * xi;
  const double mu    = 10;
  const double mu2   = mu * mu;

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, mc, mb, mt};

  // Init PARTONS instance
  PARTONS::Partons* pPartons = 0;

  // Init PARTONS application
  pPartons = PARTONS::Partons::getInstance();
  pPartons->init(argc, argv);

  // Retrieve GPD service
  PARTONS::GPDService* pGPDService = PARTONS::Partons::getInstance()->getServiceObjectRegistry()->getGPDService();

  // Create GPD module with the BaseModuleFactory
  PARTONS::GPDModule* pGPDModel = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDModule(PARTONS::GPDGK16::classId);

  // Create GPD evolution module with the BaseModuleFactory
  PARTONS::GPDEvolutionModule* pGPDEvolution = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newGPDEvolutionModule(PARTONS::GPDEvolutionApfel::classId);

  // Create alphaS module with the BaseModuleFactory
  PARTONS::RunningAlphaStrongModule* pRunningAlphaStrongModule = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newRunningAlphaStrongModule(PARTONS::RunningAlphaStrongApfel::classId);
  static_cast<PARTONS::RunningAlphaStrongApfel*>(pRunningAlphaStrongModule)->setPertOrder(PARTONS::PerturbativeQCDOrderType::NLO);
  static_cast<PARTONS::RunningAlphaStrongApfel*>(pRunningAlphaStrongModule)->setAlphasRef(asref);
  static_cast<PARTONS::RunningAlphaStrongApfel*>(pRunningAlphaStrongModule)->setMuRef(Qref);
  static_cast<PARTONS::RunningAlphaStrongApfel*>(pRunningAlphaStrongModule)->setThresholds(Thresholds);
  static_cast<PARTONS::RunningAlphaStrongApfel*>(pRunningAlphaStrongModule)->configure(ElemUtils::Parameters{});

    // Create active flavors thresholds module with the
  // BaseModuleFactory
  PARTONS::ActiveFlavorsThresholdsModule* pActiveFlavorsThresholdsModule = PARTONS::Partons::getInstance()->getModuleObjectFactory()->newActiveFlavorsThresholdsModule(PARTONS::ActiveFlavorsThresholdsVariable::classId);
  static_cast<PARTONS::ActiveFlavorsThresholdsVariable*>(pActiveFlavorsThresholdsModule)->setThresholds(Thresholds);

  // Create parameters to configure later GPDEvolutionModule
  ElemUtils::Parameters parameters;
  parameters.add(PARTONS::PerturbativeQCDOrderType::PARAMETER_NAME_PERTURBATIVE_QCD_ORDER_TYPE, PARTONS::PerturbativeQCDOrderType::LO);

  // Configure GPDEvolutionModule with previous parameters.
  pGPDEvolution->configure(parameters);

  // Link modules (set physics assumptions of your computation)
  pGPDEvolution->setRunningAlphaStrongModule(pRunningAlphaStrongModule);
  pGPDEvolution->setActiveFlavorsModule(pActiveFlavorsThresholdsModule);
  pGPDModel->setEvolQcdModule(pGPDEvolution);

  // Create a GPDKinematic(x, xi, t, MuF2, MuR2) to compute
  int nx = 10;
  double xmin = 0.001;
  double xmax = 0.8;
  double xstp = exp( log( xmax / xmin ) / ( nx - 1 ) );
  std::vector<PARTONS::GPDKinematic> gpdKinematicVector;
  for (double x = xmin; x <= xmax * 1.000001; x *= xstp)
    gpdKinematicVector.push_back(PARTONS::GPDKinematic{x, xi, t, mu2, mu2});
  PARTONS::List<PARTONS::GPDKinematic> gpdKinematic{gpdKinematicVector};

  // Run computations only for H
  std::vector<PARTONS::GPDType> vtype{PARTONS::GPDType::Type::H};
  PARTONS::List<PARTONS::GPDType> type{vtype};
  PARTONS::List<PARTONS::GPDResult> gpdResult = pGPDService->computeManyKinematic(gpdKinematic, pGPDModel, type);

  // Print results
  std::cout << std::scientific;
  std::cout << "\n#xi = " << xi << std::endl;
  std::cout << "#t = " << t << " GeV^2" << std::endl;
  std::cout << std::endl;
  std::cout << "APFEL++ through PARTONS" << std::endl;
  std::cout << "#    x             gluon            u(+)            u(-)            c(+)" << std::endl;
  for (auto const& r : gpdResult.getData())
    std::cout << r.getKinematic().getX().getValue() << "\t"
	      << r.getPartonDistribution(PARTONS::GPDType::H).getGluonDistribution().getGluonDistribution() << "\t"
	      << r.getKinematic().getX().getValue() * r.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistributionPlus() << "\t"
	      << r.getKinematic().getX().getValue() * r.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistributionMinus() << "\t"
	      << r.getKinematic().getX().getValue() * r.getPartonDistribution(PARTONS::GPDType::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::CHARM).getQuarkDistributionPlus() << "\t"
	      << std::endl;
  std::cout << std::endl;

  // Configure APFEL++
  // Running coupling
  apfel::AlphaQCD a{asref, Qref, Thresholds, 1};
  const  apfel::TabulateObject<double> Alphas{a, 100, 0.9, 100, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{{100, 1e-7, 3}, {60, 1e-1, 3}, {50, 6e-1, 3}, {50, 8e-1, 3}}};

  // Input GPDs
  const auto InGPDs = [=] (double const& x, double const& Q) -> std::map<int, double>
    {
      PARTONS::GPDResult gpdResult = pGPDService->computeSingleKinematic(PARTONS::GPDKinematic{x, xi, t, Q * Q, Q * Q}, pGPDModel);
      std::map<int, double> PhysMap{{-6, 0}, {-5, 0}, {-4, 0}, {-3, 0}, {-2, 0}, {-1, 0}, {0, 0}, {1, 0}, {2, 0}, {3, 0}, {4, 0}, {5, 0}, {6, 0}};
      PhysMap[3]  = x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE).getQuarkDistribution();
      PhysMap[2]  = x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistribution();
      PhysMap[1]  = x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistribution();
      PhysMap[0]  = gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getGluonDistribution().getGluonDistribution();
      PhysMap[-1] = x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::DOWN).getQuarkDistributionPlus() - PhysMap[1];
      PhysMap[-2] = x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::UP).getQuarkDistributionPlus() - PhysMap[2];
      PhysMap[-3] = x * gpdResult.getPartonDistribution(PARTONS::GPDType::Type::H).getQuarkDistribution(PARTONS::QuarkFlavor::Type::STRANGE).getQuarkDistributionPlus() - PhysMap[3];
      return apfel::PhysToQCDEv(PhysMap);
    };

  // Evolve and tabulate GPDs
  apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs{*(BuildDglap(InitializeGpdObjects(g, Thresholds, xi), InGPDs, sqrt(pGPDModel->getMuF2Ref()), 0, as)), 100, 1, 1000, 3};
  const auto CollGPDs = [&] (double const& mu) -> apfel::Set<apfel::Distribution> { return TabulatedGPDs.Evaluate(mu); };

  // Input PDFs (forward limit of GPDs)
  const auto InPDFs = [=] (double const& x, double const& Q) -> std::map<int, double>
    {
      PARTONS::GPDResult gpdResult = pGPDService->computeSingleKinematic(PARTONS::GPDKinematic{x, 0.0001, t, Q * Q, Q * Q}, pGPDModel);
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

  // Evolve and tabulate PDFs
  apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*(BuildDglap(InitializeGpdObjects(g, Thresholds, xi), InPDFs, sqrt(pGPDModel->getMuF2Ref()), 0, as)), 100, 1, 1000, 3};
  const auto CollPDFs = [&] (double const& mu) -> apfel::Set<apfel::Distribution> { return TabulatedPDFs.Evaluate(mu); };

  // Get evolved GPDs in the physical basis
  const std::map<int, apfel::Distribution> xgps = apfel::QCDEvToPhys(TabulatedGPDs.Evaluate(mu).GetObjects());

  std::cout << std::endl;
  std::cout << "APFEL++ stand alone" << std::endl;
  std::cout << "#    x             gluon            u(+)            u(-)            c(+)" << std::endl;
  for (auto const& r : gpdResult.getData())
    std::cout << r.getKinematic().getX().getValue() << "\t"
	      << xgps.at(0).Evaluate(r.getKinematic().getX().getValue()) << "\t"
	      << xgps.at(2).Evaluate(r.getKinematic().getX().getValue()) + xgps.at(-2).Evaluate(r.getKinematic().getX().getValue()) << "\t"
	      << xgps.at(2).Evaluate(r.getKinematic().getX().getValue()) - xgps.at(-2).Evaluate(r.getKinematic().getX().getValue()) << "\t"
	      << xgps.at(4).Evaluate(r.getKinematic().getX().getValue()) + xgps.at(-4).Evaluate(r.getKinematic().getX().getValue()) << "\t"
	      << std::endl;
  std::cout << std::endl;

  // Remove pointer references
  PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(pActiveFlavorsThresholdsModule, 0);
  pActiveFlavorsThresholdsModule = 0;

  PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(pRunningAlphaStrongModule, 0);
  pRunningAlphaStrongModule = 0;

  PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(pGPDModel, 0);
  pGPDModel = 0;

  // Close PARTONS application properly
  pPartons->close();

  return 0;
}
