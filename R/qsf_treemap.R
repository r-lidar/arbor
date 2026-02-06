qsf_treemap = function(qsf)
{
  ans = Filter(function(x) inherits(x, "qsm"), qsf)
  ans = lapply(qsf, function(x) qsm_stats(x)$stats_global)
  ans = data.table::rbindlist(ans)
  ans = ans[,-c("X_bh", "Y_bh", "Z_bh")]
  ans = sf::st_as_sf(ans, coords = c("X_root", "Y_root", "Z_root"))
  ans
}
