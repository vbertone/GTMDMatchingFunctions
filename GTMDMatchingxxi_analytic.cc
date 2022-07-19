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
 * Gluon mathing functions
 */
double C1gq(double const& y, double const& kappa)
{
  const double ky  = kappa * y;
  const double ky2 = pow(ky, 2);
  return (y <= 1 ? 2 * apfel::CF * ( 1 - pow(kappa, 2) ) * y / ( 1 - ky2 ) : 0)
    + (kappa > 1 ? - 2 * apfel::CF * ( 1 - pow(kappa, 2) ) / kappa / ( 1 - ky2 ) : 0);
}

double C1ggR(double const& y, double const& kappa)
{
  const double ky  = kappa * y;
  const double ky2 = pow(ky, 2);
  if (kappa < 1)
    return (y <= 1 ? 8 * apfel::CA * pow(kappa, 2) * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
  else
    return (y <= 1 ? apfel::CA * ( 1 + kappa * ( y + kappa * ( - 5 + 3 * ky ) ) ) / 2 / kappa / pow(1 + ky, 2) : apfel::CA * ( 1 - kappa ) * ( 1 + kappa - ( 1 - 7 * kappa ) * ky2 ) / kappa / pow(1 - ky2, 2));
}
double C1ggL()
{
  return - apfel::CA * apfel::zeta2;
}
double C1ggSPV(double const& y, double const& kappa)
{
  if (kappa > 1 && y < 1)
    return apfel::CA * ( 1 + 3 * pow(kappa, 2) ) / ( 1 - kappa * y ) / 2 / kappa;
  else
    return 0;
}
double C1ggLPV(double const& y, double const& kappa)
{
  if (kappa > 1 && y < 1 && kappa * y < 1)
    return apfel::CA * ( 1 + 3 * pow(kappa, 2) ) * log( kappa * ( 1 - kappa * y ) / ( kappa - 1 ) ) / 2 / pow(kappa, 2);
  else
    return 0;
}

/*
 * This code computes F11 of the gluon as a function of x at fixed kT,
 * t, and Q for a vector of values in xi. Differently from
 * GTMDMatchingkTxi.cc, this code also computes the convolution
 * integrals directly without interpolating. The reason ins that the
 * interpolating code GTMDMatchingkTxi.cc does not do a great job
 * around x = xi where GTMDs are singular.
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
  const double xi    = 0.1;
  const int ipt      = 2;

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

  // Configure APFEL++
  apfel::SetVerbosityLevel(2);
  apfel::Banner();

  // Running coupling
  apfel::AlphaQCD a{asref, Qref, Thresholds, 1};
  const  apfel::TabulateObject<double> Alphas{a, 100, 0.9, 100, 3};
  const auto as = [=] (double const& mur) -> double{ return Alphas.Evaluate(mur); };

  // x-space grid
  const apfel::Grid g{{{100, 1e-7, 3}, {100, 1e-2, 3}, {100, 6e-1, 3}, {100, 8e-1, 3}}};

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

  // Get GTMD objects
  const auto GtmdObj = apfel::InitializeGtmdObjects(g, Thresholds, xi);

  // Now compute GTMDs in bT space
  const auto EvGTMDs = BuildGtmds(GtmdObj, CollGPDs, as, ipt);

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
  const apfel::Set<apfel::Distribution> F = DEObj.transform(xGb, kT);

  // Sudakov form factor
  const std::function<double(double const&, double const&, double const&)> Sf = GluonEvolutionFactor(GtmdObj, as, ipt);

  // Tabulation parameters
  const int nx = 1000;
  const double xmin = 0.01;
  const double xmax = 0.9;
  const double xstp = exp( log( xmax / xmin ) / ( nx - 1 ) );
  std::cout << std::scientific;
  std::cout << "#    x         Numerical      Analytic" << std::endl;
  for (double x = xmin; x <= xmax * ( 1 + 1e-5 ); x *= xstp)
    {
      // kappa parameter
      const double kappa = xi / x;

      // Analytic computation
      const std::function<double(double const&)> xGbAn = [=] (double const& bT) -> double
	{
	  const double b = bs(bT, mu);
	  const double mub = 2 * exp(- apfel::emc) / b;
	  const auto gpsd = TabulatedGPDs.Evaluate(mub);
	  const double gluxi = gpsd.at(0).Evaluate(xi);
	  const std::function<double(double const&)> f = [=] (double const& y) -> double
	  {
	    const double ky  = kappa * y;
	    const double qrk = gpsd.at(1).Evaluate(x/y);
	    const double glu = gpsd.at(0).Evaluate(x/y);
	    return C1gq(y, kappa) * qrk
	    + C1ggR(y, kappa) * glu
	    + C1ggSPV(y, kappa) * ( glu - gluxi * ( 1 + (ky > 1 ? ( 1 - ky ) / ky : 0) ) );
	  };
	  const apfel::Integrator I{f};
	  const double coup = (ipt > 1 ? as(mub) / apfel::FourPi : 0);
	  const double mGPD = gpsd.at(0).Evaluate(x)
	  + coup * ( I.integrate(x, 1.1, apfel::eps7) + I.integrate(1.1, x * 1e7, apfel::eps7)
		     + C1ggL() * gpsd.at(0).Evaluate(x)
		     + C1ggLPV(x, kappa) * gluxi );
	  return bT * Sf(b, mu, mu * mu) * mGPD * NPFunc->Evaluate(x, bT, mu * mu, 0);
	};

      const double nm = QCDEvToPhys(F.GetObjects()).at(0).Evaluate(x);
      const double an = DEObj.transform(xGbAn, kT);
      std::cout << x << "  "
		<<  nm << "  " << an << "  " << nm / an << "  "
		<< std::endl;
    }

  // Remove pointer references
  PARTONS::Partons::getInstance()->getModuleObjectFactory()->updateModulePointerReference(pGPDModel, 0);
  pGPDModel = 0;

  // Close PARTONS application properly
  pPartons->close();

  return 0;
}
