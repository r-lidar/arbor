#include <string>

#include "QSM.h"


class QSF
{
public:
  QSF() = default;
  void add_qsm(const std::string& name, const QSM& q);
  void write(const std::string& dir, const std::string& format, bool binary = true) const;

private:
  std::unordered_map<std::string, QSM> qsm_;
};
