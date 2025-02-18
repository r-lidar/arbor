cc_export_labels = function(heights, circles, file)
{
  labels1 = heights[,1:3]
  labels1$label = paste0(round(heights$Z, 1), " m")

  labels2 = circles[,1:3]
  labels2$label = paste0(round(circles$radius*2*100, 1), " cm")
  names(labels2) = names(labels1)

  labels = rbind(labels1, labels2)

  circles = sf::st_as_sf(circles, coords = c("center_x", "center_y"))
  circles = sf::st_buffer(circles, circles$radius)
  circles = sf::st_cast(circles, "POINT")

  path = tools::file_path_sans_ext(file)

  flabel = paste0(path, "_labels.txt")
  frings = paste0(path, "_rings.shp")

  data.table::fwrite(labels, flabel)
  sf::st_write(circles["center_z"], frings, append = FALSE)
}
