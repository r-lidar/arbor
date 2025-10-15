#' Transfer attributes between point clouds
#'
#' Transfers one or more attributes from a decimated (sparser) point cloud
#' to a denser point cloud using a 1-nearest neighbor (1-NN) approach.
#' Each point in the target point cloud receives the attribute value of
#' its nearest neighbor in the source point cloud.
#'
#' @param from A `LAS` object representing the source (typically decimated) point cloud
#'   containing the attributes to transfer.
#' @param to A `LAS` object representing the target (typically denser) point cloud
#'   that will receive the transferred attributes.
#' @param attributes A character vector specifying the names of the attributes
#'   to transfer from `from` to `to`.
#'
#' @details
#' The function performs a 1-NN search between the `from` and `to` point clouds
#' using `lidR::knnx()`. For each attribute listed in `attributes`, the values are
#' assigned from the nearest point in `from`. If available, the function also
#' propagates the corresponding attribute description from the Extra Bytes
#' metadata (`VLR$Extra_Bytes`).
#'
#' @return A `LAS` object identical to `to`, but with the requested attributes added.
#'
#' @seealso [lidR::knnx()], [lidR::add_lasattribute()]
#'
#' @examples
#' \dontrun{
#' las <- transfer_attributes(las_decimated, las, attributes = "treeID", "foliage")
#' }
#'
#' @export
transfer_attributes = function(from, to, attributes)
{
  # Sanity check
  stopifnot(all(attributes %in% names(from@data)))

  # Compute 1-NN mapping (index + distance)
  nn   <- lidR::knnx(from, to, k = 1)
  idx  <- nn[[1]]
  dist <- nn[[2]]

  # For each requested attribute, transfer it
  for (attr in attributes)
  {
    eba <- from@header@VLR$Extra_Bytes[["Extra Bytes Description"]][[attr]]
    values <- from@data[[attr]][idx]
    to@data[[attr]] <- values

    if (!is.null(eba))
      to <- lidR::add_lasattribute(to,  name = attr,  desc = eba$description)
  }

  return(to)
}
