/**
 * @file QSF.h
 * Project: Arbor
 *
 * Copyright (C) 2026 Jean-Romain Roussel (r-lidar) <info @ r-lidar.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef QSF_H
#define QSF_H

#include "QSM.h"
#include <string>
#include <unordered_map>
#include <filesystem>

namespace arbor::qsm {

class QSF
{
public:
  QSF() = default;
  void add_qsm(const QSM& q);

  // Writes every QSM held by this QSF.
  //
  // - If 'path' is (or would be) a directory, one file per QSM is written
  //   to 'path/format/<name-or-id>.<format>' (legacy behaviour).
  // - If 'path' names a file whose extension matches 'format', e.g.
  //   write("trees.obj", "obj"), write("trees.ply", "ply") or
  //   write("trees.stl", "stl"), all QSMs are written into that single
  //   combined file instead. Combined single-file output is only supported
  //   for the 'obj', 'ply' and 'stl' formats. How well each QSM stays
  //   identifiable inside that combined file depends on what the format
  //   itself can express:
  //     - OBJ:          each QSM becomes its own named 'o' object. Fully
  //                      supported and preserved by essentially all OBJ
  //                      readers (e.g. Blender).
  //     - PLY:          all QSMs share one vertex/face list, but every face
  //                      carries a 'qsm_id' property and the header records
  //                      one comment per QSM (id + name), so QSMs can be
  //                      recovered by filtering on that property. Generic
  //                      PLY viewers will show a single unnamed mesh.
  //     - STL (ASCII):  each QSM is written as its own
  //                      'solid <name> ... endsolid <name>' block. This is
  //                      only a de facto convention, not part of the STL
  //                      standard: many tools import multiple solids as
  //                      separate named objects, but some only read the
  //                      first block or silently merge all of them, so
  //                      verify with the target software before relying on
  //                      it.
  //     - STL (binary): the binary format has no concept of an object name
  //                      or grouping at all - it is a flat, unlabelled
  //                      triangle list. All QSMs are therefore merged into
  //                      a single anonymous mesh with no way to recover
  //                      which triangle belongs to which tree. Use ASCII
  //                      STL, OBJ, or PLY instead if that distinction
  //                      matters.
  void write(const std::string& path, const std::string& format, bool binary = true) const;
  static QSF read(const std::string& path); // Reads a .qsf manifest file and returns the QSF it describes.

  const std::unordered_map<int, QSM>& get_qsm_map() const { return qsm_; }

private:
  std::unordered_map<int, QSM> qsm_;

  void write_obj(const std::filesystem::path& file) const;
  void write_ply(const std::filesystem::path& file, bool binary) const;
  void write_stl(const std::filesystem::path& file, bool binary) const;
  void write_qsf(const std::filesystem::path& qsm_dir, const std::unordered_map<int, std::filesystem::path>& written_files) const;
};

}

#endif