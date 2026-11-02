#include "allometry.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>

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

std::string Griese2025Allometry::authors() const { return "Griese, N., Ritzert, M. & N\u00f6lke, N."; }
std::string Griese2025Allometry::reference() const { return "https://doi.org/10.1038/s41597-025-06421-7"; }
std::string Griese2025Allometry::country() const { return "Germany"; }

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

std::string CostaCysneiros2020::authors() const { return "Costa Cysneiros, V., Pelissari, A. L., Gaui, T. D., Fiorentin, L. D., de Carvalho, D. C., Silveira Filho, T. B. & do Amaral Machado, S."; }
std::string CostaCysneiros2020::reference() const { return "https://doi.org/10.1139/cjfr-2020-0060"; }
std::string CostaCysneiros2020::country() const { return "Brazil"; }

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

std::string Chenge2020::authors() const { return "Chenge, I. B."; }
std::string Chenge2020::reference() const { return "https://doi.org/10.1016/j.tfp.2020.100051"; }
std::string Chenge2020::country() const { return "Nigeria"; }

double FSMiombo::H_vs_DBH(double dbh) const
{
  constexpr double a = 2.759;
  constexpr double b = 0.465;

  return a * std::pow(dbh*100, b);
}

double FSMiombo::DBH_vs_H(double H) const
{
  constexpr double a = 2.759;
  constexpr double b = 0.465;

  return std::pow(H / a, 1.0 / b)/100;
}


std::string FSMiombo::authors() const { return "r-lidar, inc"; }
std::string FSMiombo::reference() const { return "Unpublished"; }
std::string FSMiombo::country() const { return "Zambia"; }

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

std::string CacaoAllometry::authors() const { return "r-lidar, inc"; }
std::string CacaoAllometry::reference() const { return "Unpublished"; }
std::string CacaoAllometry::country() const { return "Unknown"; }

namespace
{
  using Factory = std::unique_ptr<Allometry>(*)();

  template <typename T>
  std::unique_ptr<Allometry> make() { return std::make_unique<T>(); }

  // Ordered list of available models. This is the single source of truth for
  // both lookup by name and listing (e.g. available_allometries() in R)
  const std::vector<std::pair<std::string, Factory>>& registry()
  {
    static const std::vector<std::pair<std::string, Factory>> reg = {
      { "CostaCysneiros2020", &make<CostaCysneiros2020> },
      { "Griese2025",         &make<Griese2025Allometry> },
      { "Chenge2020",         &make<Chenge2020> },
      { "FSMiombo",           &make<FSMiombo> }
      //{ "Cacao",              &make<CacaoAllometry> }
    };
    return reg;
  }
}

std::unique_ptr<Allometry> AllometryDataBase::getAllometry(const std::string& name)
{
  for (const auto& [key, factory] : registry())
  {
    if (key == name) return factory();
  }

  throw std::invalid_argument("Unknown allometry model name: " + name);
}

std::vector<std::string> AllometryDataBase::names()
{
  std::vector<std::string> out;
  for (const auto& entry : registry()) out.push_back(entry.first);
  return out;
}

std::vector<AllometryDataBase::Entry> AllometryDataBase::entries()
{
  std::vector<Entry> out;
  for (const auto& [key, factory] : registry())
  {
    auto model = factory();
    out.push_back({ key, model->authors(), model->reference(), model->country() });
  }
  return out;
}
