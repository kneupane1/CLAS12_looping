
#ifndef DT_H_GUARD
#define DT_H_GUARD
static const double MASS_P = 0.93827203;
static const double MASS_N = 0.93956556;
static const double MASS_E = 0.000511;
static const double MASS_PIP = 0.13957018;
static const double MASS_PIM = 0.13957018;
static const double MASS_PI0 = 0.1349766;
static const double MASS_KP = 0.493677;
static const double MASS_KM = 0.493677;
static const double MASS_G = 0.0;
static const double MASS_OMEGA = 0.78265;
static const double CLAS12_E = 10.7;
const double c_special_units = 29.9792458;

double vertex_time(double sc_time, double sc_pathlength,
                   double relatavistic_beta) {
  return (sc_time - (sc_pathlength / (relatavistic_beta * c_special_units)));
}

double delta_t(double electron_vertex_time, double mass, double momentum,
               double sc_t, double sc_r) {
  double relatavistic_beta =
      1.0 / sqrt(1.0 + ((mass / momentum) * (mass / momentum)));
  return (electron_vertex_time - vertex_time(sc_t, sc_r, relatavistic_beta));
}
#endif
