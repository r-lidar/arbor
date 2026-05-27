#ifndef ALLOMETRY_H
#define ALLOMETRY_H

#include <memory>
#include <string_view>

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

class AllometryDataBase
{
public:
  static std::unique_ptr<Allometry> getAllometry(std::string_view name);
};

#endif // ALLOMETRY_H
