#include "allometry.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <functional>
#include <unordered_map>

// --- Griese2025Allometry Implementation ---
double Griese2025Allometry::H_vs_DBH(double dbh) const
{
  return 36.03 * std::pow(1.0 - std::exp(-0.05 * dbh), 1.1);
}

double Griese2025Allometry::DBH_vs_H(double H) const
{
  // Transition height chosen just below the 36.03 asymptote
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

double CostaCysneiros2020::H_vs_DBH(double dbh) const
{
  constexpr double b0 = 1.029;
  constexpr double b1 = 0.567;

  return std::exp(b0) * std::pow(dbh*100, b1);
}

double CostaCysneiros2020::DBH_vs_H(double H) const
{
  constexpr double b0 = 1.029;
  constexpr double b1 = 0.567;

  return std::exp((std::log(H) - b0) / b1)/100;
}

double Chenge2020::H_vs_DBH(double dbh) const
{
  constexpr double a = 5.135;
  constexpr double b = 0.451;

  return a * std::pow(dbh*100, b);
}

double Chenge2020::DBH_vs_H(double H) const
{
  constexpr double a = 5.135;
  constexpr double b = 0.451;

  return std::pow(H / a, 1.0 / b)/100;
}


double CacaoAllometry::H_vs_DBH(double dbh) const
{
  return 36.03 * std::pow(1.0 - std::exp(-0.05 * dbh), 1.1);
}

double CacaoAllometry::DBH_vs_H(double H) const
{
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

  return 3*(DBH_cm / 100.0);
}

std::unique_ptr<Allometry> AllometryDataBase::getAllometry(const std::string& name)
{
  using Factory = std::function<std::unique_ptr<Allometry>()>;

  static const std::unordered_map<std::string_view, Factory> registry = {
    { "Griese2025",         []() { return std::make_unique<Griese2025Allometry>(); } },
    { "Chenge2020",         []() { return std::make_unique<Chenge2020>(); } },
    { "CostaCysneiros2020", []() { return std::make_unique<CostaCysneiros2020>(); } },
    { "Cacao",              []() { return std::make_unique<CacaoAllometry>(); } }
  };

  auto it = registry.find(name);

  if (it != registry.end())
  {
    return it->second();
  }

  throw std::invalid_argument("Unknown allometry model name: " + std::string(name));
}
