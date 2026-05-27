#include "allometry.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>

#include <cmath>
#include <algorithm>

// --- Griese2025Allometry Implementation ---
double Griese2025Allometry::H_vs_DBH(double dbh) const
{
  // dbh input is likely in meters or cm depending on your framework,
  // assuming matching units with your original code.
  return 36.03 * std::pow(1.0 - std::exp(-0.05 * dbh), 1.1);
}

double Griese2025Allometry::DBH_vs_H(double H) const
{
  // Transition height chosen just below the 36.03 asymptote
  // to ensure numerical stability and a clean slope.
  const double H_tr = 35.5;

  double DBH_cm;

  if (H < H_tr)
  {
    double ratio = H / 36.03;
    DBH_cm = -1.0 / 0.05 * std::log(1.0 - std::pow(ratio, 1.0 / 1.1));
  }
  else
  {
    const double ratio_tr = H_tr / 36.03;
    const double dbh_tr = -1.0 / 0.05 * std::log(1.0 - std::pow(ratio_tr, 1.0 / 1.1));

    const double exp_term = std::exp(-0.05 * dbh_tr);
    const double dH_ddbh = 36.03 * 1.1 * std::pow(1.0 - exp_term, 0.1) * (0.05 * exp_term);

    const double m = 1.0 / dH_ddbh;
    const double b = dbh_tr - m * H_tr;

    DBH_cm = m * H + b;
  }

  return DBH_cm / 100.0;
}
// --- AllometryDataBase Implementation ---
std::unique_ptr<Allometry> AllometryDataBase::getAllometry(std::string_view name)
{
  if (name == "Griese2025")
  {
    return std::make_unique<Griese2025Allometry>();
  }

  throw std::invalid_argument("Unknown allometry model name: " + std::string(name));
}
