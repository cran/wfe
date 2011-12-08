
Index <- function(index_name, uniq_index_name, len_u_index, len_data) {

  return(.C("Index", as.integer(index_name), as.integer(uniq_index_name), as.integer(len_u_index),
            as.integer(len_data),
            result = integer(len_data), package = "wfe")$result)

}



Vectorize <- function(W, time.index, dyad.index, n.row) {

  return(.C("Vectorize", as.double(W), as.integer(nrow(W)), as.integer(ncol(W)),
            as.integer(time.index), as.integer(dyad.index), as.integer(n.row),
            results = double(n.row), package = "wfe")$results)

}



Transform <- function(y, treat, pscore) {

  return(.C("Transform", as.double(y), as.integer(length(y)), as.integer(treat),
            as.double(pscore), ytrans = double(length(y)), package = "wfe")$ytrans)

}



GenTime <- function(unit_index, len_data, len_u_index) {
  return(.C("GenTime", as.integer(unit_index), as.integer(len_data),
            as.integer(len_u_index), time_index = double(len_data), package = "wfe")$time_index)
}



GenWeights <- function(unit_index, time_index, tr, C_it, len_data, len_u_index, len_t_index, ate, att, size) {
  return(.C("GenWeights", as.integer(unit_index), as.integer(time_index), as.integer(tr), as.integer(C_it),
            as.integer(len_data), as.integer(len_u_index), as.integer(len_t_index),
            as.integer(ate), as.integer(att),
            weight = double(len_u_index*len_t_index), package = "wfe")$weight)
}

GenWeightsFD <- function(unit_index, time_index, tr, C_it, len_data, len_u_index, len_t_index, ate, att) {
  return(.C("GenWeightsFD", as.integer(unit_index), as.integer(time_index), as.integer(tr), as.integer(C_it),
            as.integer(len_data), as.integer(len_u_index), as.integer(len_t_index),
            as.integer(ate), as.integer(att),
            weightfd = double(len_u_index*len_t_index), package = "wfe")$weightfd)
}


Demean <- function(var_name, index, len_index, len_data) {
  return(.C("Demean", as.double(var_name), as.integer(index), as.integer(len_index),
            as.integer(len_data), demean = double(len_data), package = "wfe")$demean)
}




WDemean <- function(var_name, weight, index, len_index, len_data) {
  return(.C("WDemean", as.double(var_name), as.double(weight), as.integer(index), as.integer(len_index),
            as.integer(len_data), Wdemean = double(len_data), package = "wfe")$Wdemean)
}




WWDemean <- function(var_name, weight, index, len_index, len_data) {
  return(.C("WWDemean", as.double(var_name), as.double(weight), as.integer(index), as.integer(len_index),
            as.integer(len_data), WWdemean = double(len_data), package = "wfe")$WWdemean)
}





Wdemean <- function (x, w, unit.index, time.index, data) {
  data$u.index <- Index(data[,unit.index])
  data$t.index <- Index(data[,time.index])
  uniq.unit <- unique(data$u.index)
  uniq.time <- unique(data$t.index)
  wdemean.x <- c()
  wdemean.x <- new.tr.x <- matrix(NA, nrow = length(uniq.time), ncol = length(uniq.unit))
  for (i in 1:length(uniq.unit)){
    sub.i <- data[data$u.index == uniq.unit[i], ]
    nr <- nrow(sub.i)
    denom <- as.numeric(sum(sub.i[,w]))
    x.star <- as.numeric((sub.i[,x]%*%sub.i[,w])/(denom))
    tr.x <- as.vector(sub.i[,x] - rep(x.star, nr))
    wdemean.x[,i] <- tr.x
      # sqrt(w)* (weighted demean) for SE
      for (j in 1:length(uniq.time)){
       new.tr.x[j,i] <- sqrt(sub.i[,w][j])*(sub.i[,x][j] - x.star)
      }
  }
  w.demean <- as.vector(wdemean.x)
  w.tr.x <- as.vector(new.tr.x)
  list(w.demeaned = wdemean.x, new.tr.x = w.tr.x)
}
