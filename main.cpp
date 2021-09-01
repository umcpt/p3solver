
#include "material.hpp"
#include "slab_geom.hpp"
#include "p3_solver.hpp"
#include "eigen_solver.hpp"
#include "depletion_solver.hpp"

int main() {

  // create isotopes with two group xs
  const auto Egroups = 2;

  // ----> h1
  IsotopeXS h1;
  h1.name    = "h1";
  h1.sigmaT  = { 20.0, 38.0 };
  h1.sigmaF  = { 0.0, 0.0 };
  h1.nubar   = { 0.0, 0.0 };
  h1.chi     = { 1.0, 0.0 };
  h1.sigmaS0 = { 20.0 * 0.8, 0.0, 20.0 * 0.2, 37.7 };
  h1.sigmaS1 = { 20.0 * 0.667, 0.0, 0.0, 37.7 * 0.667 };
  h1.sigmaS2 = { 0.0, 0.0, 0.0, 0.0 };
  h1.sigmaS3 = { 0.0, 0.0, 0.0, 0.0 };

  // ----> o16
  IsotopeXS o16;
  o16.name    = "o16";
  o16.sigmaT  = { 4.2, 4.2 };
  o16.sigmaF  = { 0.0, 0.0 };
  o16.nubar   = { 0.0, 0.0 };
  o16.chi     = { 1.0, 0.0 };
  o16.sigmaS0 = { 4.2 * 0.98, 0.0, 4.2 * 0.02, 4.2 };
  o16.sigmaS1 = { 4.2 * 0.04, 0.0, 0.0, 4.2 * 0.04 };
  o16.sigmaS2 = { 0.0, 0.0, 0.0, 0.0 };
  o16.sigmaS3 = { 0.0, 0.0, 0.0, 0.0 };

  // ----> u235
  IsotopeXS u235;
  u235.name    = "u235";
  u235.sigmaT  = { 16.0, 700.0 };
  u235.sigmaF  = { 0.0,  577.0 };
  u235.nubar   = { 2.80,   2.45 };
  u235.chi     = { 1.0, 0.0 };
  u235.sigmaS0 = { 15.0, 0.0, 0.0, 15.0 };
  u235.sigmaS1 = { 0.0, 0.0, 0.0, 0.0 };
  u235.sigmaS2 = { 0.0, 0.0, 0.0, 0.0 };
  u235.sigmaS3 = { 0.0, 0.0, 0.0, 0.0 };

  // ----> u238
  IsotopeXS u238;
  u238.name    = "u238";
  u238.sigmaT  = { 18.0, 12.0 };
  u238.sigmaF  = { 0.5,   0.0 };
  u238.nubar   = { 2.80,  0.0 };
  u238.chi     = { 1.0, 0.0 };
  u238.sigmaS0 = { 15.0, 0.0, 0.0, 10.0 };
  u238.sigmaS1 = { 0.0, 0.0, 0.0, 0.0 };
  u238.sigmaS2 = { 0.0, 0.0, 0.0, 0.0 };
  u238.sigmaS3 = { 0.0, 0.0, 0.0, 0.0 };

  // ----> i135
  IsotopeXS i135;
  i135.name    = "i135";
  i135.sigmaT  = { 10.0, 150.0 };
  i135.sigmaF  = { 0.0,  0.0 };
  i135.nubar   = { 0.0,  0.0 };
  i135.chi     = { 1.0, 0.0 };
  i135.sigmaS0 = { 9.0, 0.0, 0.0, 15.0 };
  i135.sigmaS1 = { 0.0, 0.0, 0.0, 0.0 };
  i135.sigmaS2 = { 0.0, 0.0, 0.0, 0.0 };
  i135.sigmaS3 = { 0.0, 0.0, 0.0, 0.0 };

  i135.fyield = 0.06386;
  i135.decay  = std::make_pair( 2.875e-5, "xe135" );
  i135.sigmaGamma = std::make_pair( xs1d_table( { 1.0, 135.0 } ), "i136" );

  // ----> xe135
  IsotopeXS xe135;
  xe135.name    = "xe135";
  xe135.sigmaT  = { 10.0, 3.0e6 };
  xe135.sigmaF  = { 0.0,  0.0 };
  xe135.nubar   = { 0.0,  0.0 };
  xe135.chi     = { 1.0, 0.0 };
  xe135.sigmaS0 = { 9.0, 0.0, 0.0, 3.0e5 };
  xe135.sigmaS1 = { 0.04, 0.0, 0.0, 1.0e3 };
  xe135.sigmaS2 = { 0.0, 0.0, 0.0, 0.0 };
  xe135.sigmaS3 = { 0.0, 0.0, 0.0, 0.0 };

  xe135.fyield = 0.00228;
  xe135.decay  = std::make_pair( 2.092e-5, "cs135" );
  xe135.sigmaGamma = std::make_pair( xs1d_table( { 1.0, 2.7e6 } ), "xe136" );

  // water
  Material water;
  water.atom_den = { 0.0223, 0.0112 };
  water.isotope  = { std::make_shared<IsotopeXS>( h1 ),
                     std::make_shared<IsotopeXS>( o16 ) };

  // water(75%) + uo2(25%) w/ u at 4% 235u fuel mix
  Material fuel_mix;
  fuel_mix.atom_den = { 0.012639, 0.010532, 0.000084, 0.002022, 0.0, 0.0 };
  fuel_mix.isotope  = { std::make_shared<IsotopeXS>( h1 ),
                        std::make_shared<IsotopeXS>( o16 ),
                        std::make_shared<IsotopeXS>( u235 ),
                        std::make_shared<IsotopeXS>( u238 ),
                        std::make_shared<IsotopeXS>( i135 ),
                        std::make_shared<IsotopeXS>( xe135 ) };
  fuel_mix.deplete = true;

  // construct geometry (20 cm - 60 cm - 20 cm, water - fuel_mix - water )
  // fuel is in six 10-cm zones for depletion
  auto fuel_mixes = std::vector<Material>( 6, fuel_mix );

  SlabGeometry geom;
  geom.ints = 2000;
  geom.xcoords = { 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0,
                   60.0, 65.0, 70.0, 75.0, 80.0, 100.0 };
  geom.mat = { std::make_shared<Material>( water ),
               std::make_shared<Material>( fuel_mixes[0] ),
               std::make_shared<Material>( water ),
               std::make_shared<Material>( fuel_mixes[1] ),
               std::make_shared<Material>( water ),
               std::make_shared<Material>( fuel_mixes[2] ),
               std::make_shared<Material>( water ),
               std::make_shared<Material>( water ),
               std::make_shared<Material>( fuel_mixes[3] ),
               std::make_shared<Material>( water ),
               std::make_shared<Material>( fuel_mixes[4] ),
               std::make_shared<Material>( water ),
               std::make_shared<Material>( fuel_mixes[5] ),
               std::make_shared<Material>( water ) };
/*
  geom.xcoords = { 20.0, 80.0, 100.0 };
  geom.mat = { std::make_shared<Material>( water ),
               std::make_shared<Material>( fuel_mix ),
               std::make_shared<Material>( water ) };
*/

  std::vector<std::vector<double>> source( geom.ints+1, std::vector<double>( 2, 0.0 ));
  for ( auto i = round(0.2*geom.ints) ; i <= round(0.8*geom.ints) ; ++i ) {
    source[i][0] = 1.0/(80.0-20.0);
  }

  auto geom_ptr = std::make_shared<SlabGeometry>( geom );

  auto p3sol   = std::make_shared<P3Solver>( Egroups, geom_ptr );
  auto eigsol  = std::make_shared<EigenSolver>( Egroups, geom_ptr, p3sol );
  auto deplsol = std::make_shared<DepletionSolver>( Egroups, geom_ptr, eigsol );

  eigsol->solve( source );

/*
  for ( auto n = 0 ; n < 1 ; ++n ) {
    deplsol->solve( source, 4.0 * 3600.0, 1.0e3 );
  }
*/
  return 0;
}
