#ifndef ALLOMETRY_H
#define ALLOMETRY_H

#include <memory>

class Allometry
{
public:
  virtual ~Allometry() = default;
  virtual double H_vs_DBH(double dbh) const = 0;
  virtual double DBH_vs_H(double H) const = 0;
};

// Griese, N., Ritzert, M. & Nölke, N. A large dataset of labelled single tree point clouds, QSMs and
// tree graphs. Sci Data 12, 1953 (2025). https://doi.org/10.1038/s41597-025-06421-7
class Griese2025Allometry : public Allometry
{
public:
  double H_vs_DBH(double dbh) const override;
  double DBH_vs_H(double H) const override;
};

// Vinicius Costa Cysneiros, Allan Libanio Pelissari, Tatiana Dias Gaui, Luan Demarco Fiorentin,
// Daniel Costa de Carvalho, Telmo Borges Silveira Filho, and Sebastião do Amaral Machado
// Canadian Journal of Forest Research 2020 50:12, 1289-1298 10.1139/cjfr-2020-0060
class CostaCysneiros2020 : public Allometry
{
public:
  double H_vs_DBH(double dbh) const override;
  double DBH_vs_H(double H) const override;
};


// Chenge, I. B. (2021). Height–diameter relationship of trees in Omo strict nature forest reserve, Nigeria.
// Trees, Forests and People, 3, 100051. https://doi.org/10.1016/j.tfp.2020.100051
class Chenge2020 : public Allometry
{
public:
  double H_vs_DBH(double dbh) const override;
  double DBH_vs_H(double H) const override;
};

class CacaoAllometry : public Allometry
{
public:
  double H_vs_DBH(double dbh) const override;
  double DBH_vs_H(double H) const override;
};

class AllometryDataBase
{
public:
  static std::unique_ptr<Allometry> getAllometry(const std::string& name);
};

#endif // ALLOMETRY_H
