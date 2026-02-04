#include "QSF.h"

#include <string>
#include <filesystem>
#include <stdexcept>

namespace fs = std::filesystem;

void QSF::add_qsm(const std::string& name, const QSM& q)
{
  qsm_[name] = q;
}

void QSF::write(const std::string& dir, const std::string& format, bool binary) const
{
  if (dir.empty())
    throw std::invalid_argument("QSF::write: output directory is empty");

  if (format.empty())
    throw std::invalid_argument("QSF::write: format is empty");

  fs::path base_dir(dir);

  // Create base directory if needed
  if (!fs::exists(base_dir))
  {
    fs::create_directories(base_dir);
  }
  else if (!fs::is_directory(base_dir))
  {
    throw std::runtime_error("QSF::write: path exists but is not a directory");
  }

  // Optional subfolder per format
  fs::path out_dir = base_dir / format;

  if (!fs::exists(out_dir))
    fs::create_directories(out_dir);

  // Write each QSM
  for (const auto& [key, qsm] : qsm_)
  {
    if (key.empty()) continue;

    fs::path filename = out_dir / (key + "." + format);
    qsm.write(filename.string(), binary);
  }
}
