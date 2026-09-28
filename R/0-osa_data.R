# Fold a raw composition proportion row for tail accumulation, mirroring the
# per-sex-block young/old fold in ceattle.cpp (Slot 2). `prop` is the raw
# proportion vector of length `nblk * nbins_blk` (block-major: block 0's bins,
# then block 1's for joint-sex). Bins below `yng` fold into the `yng` bin and
# above `old` into the `old` bin, within each sex block, yielding an
# `nblk * (old - yng + 1)` vector. yng/old are 1-based, already clamped by the
# caller (1 <= yng <= old <= nbins_blk). Returns the folded raw proportions; the
# caller adds the offset once per folded bin via append_composition(), matching
# the cpp order (fold raw props, then += comp_prop_offset).
.fold_comp_bins <- function(prop, nbins_blk, nblk, yng, old) {
  nkeep  <- old - yng + 1L
  folded <- numeric(nblk * nkeep)
  for (b in seq_len(nblk) - 1L) {
    for (j in seq_len(nbins_blk) - 1L) {   # 0-based source bin within the block
      tgt <- j
      if (j < yng - 1L) tgt <- yng - 1L    # fold young tail into the yng bin
      if (j > old - 1L) tgt <- old - 1L    # fold old tail into the old bin
      k <- b * nkeep + (tgt - (yng - 1L)) + 1L
      folded[k] <- folded[k] + prop[b * nbins_blk + j + 1L]
    }
  }
  folded
}

#' Build the flat observation vector and metadata for OSA residuals
#'
#' @description
#' Internal helper used by [rearrange_data()] to assemble the inputs the TMB
#' template needs to compute one-step-ahead (OSA) residuals via
#' [TMB::oneStepPredict()].
#'
#' OSA residuals require every observation that enters the likelihood to live in
#' a single flat vector (`obsvec`) with a companion `keep` indicator. This
#' function walks the already-built control/observation matrices and produces:
#'
#' \itemize{
#'   \item `obsvec`, the flat vector of observations. Aggregate (catch and
#'     index) observations are stored as log(observation) because their
#'     likelihood is lognormal; composition (comp) and conditional-age-at-length
#'     (caal) observations are stored as bin counts, `(proportion + 1e-5) * N`,
#'     matching exactly what the TMB likelihood forms during fitting.
#'   \item `obs_ctl`, a data frame with one row per element of `obsvec`,
#'     mapping each position back to its data type, fleet, species, year,
#'     age/length bin, etc., so residuals stay interpretable. R-side metadata
#'     only; removed before the data list is passed to TMB.
#'   \item `*_obsvec_idx`, per observation-row index vectors. For aggregate
#'     series this is the 0-based `obsvec` position of each row; for composition
#'     and caal it is the 0-based `obsvec` position of the row's FIRST bin (the
#'     template reads the row's bins as `obsvec.segment(start, n_bins)`). `-1`
#'     marks rows excluded from the likelihood.
#' }
#'
#' The inclusion rules and per-row bin counts below must match the guards in the
#' TMB model exactly, so that every observation the model evaluates has a
#' valid `obsvec` position. All fitted observation types are handled: aggregate
#' catch and index, comp and caal composition, and predator diet composition.
#'
#' @param data_list A `data_list` whose `*_ctl`/`*_obs`/`*_n` matrices have
#'   already been built (and coerced to matrices) by [rearrange_data()].
#' @param build_osa Logical. When `TRUE`, build the full OSA observation data
#'   for every type (aggregate, composition, caal, and diet) so
#'   [osa_residuals()] can be computed. When `FALSE` (the default), only the
#'   aggregate index/catch entries the TMB template always reads are built and
#'   the (much larger) composition/caal/diet metadata is skipped, this is the
#'   fast path for simulation testing (e.g. [run_mse()]), where the fitted
#'   objective is identical but OSA composition residuals are not produced.
#'
#' @details The composition proportion offset is read from `data_list$comp_offset`
#'   (defaulting to `1e-5`), so the comp/caal bin counts are `(proportion +
#'   comp_offset) * Neff`, the same offset the TMB likelihood applies when
#'   fitting. Set it via `fit_control(comp_offset = )` or on `data_list` directly.
#'
#' @return The input `data_list` with `obsvec`, `obs_ctl`, `osa_mode`,
#'   `comp_offset`, and the per-type `*_obsvec_idx` vectors added.
#'
#' @keywords internal
build_osa_data <- function(data_list, build_osa = FALSE) {

  # Proportion offset added to comp/CAAL bins before the likelihood, so the fit
  # and the OSA obsvec use one value. Filled by switch_check() and by
  # fit_control(comp_offset =) on the fit path, and here (default 1e-5) on the
  # exported rearrange_data() path, which runs neither.
  comp_offset <- data_list$comp_offset
  if (is.null(comp_offset)) comp_offset <- 1e-5
  comp_offset <- as.numeric(comp_offset)[1]
  if (!is.finite(comp_offset) || comp_offset < 0) {
    stop("`comp_offset` must be a single non-negative number.")
  }

  # Lognormal bias-correction toggles (DATA_SCALARs read by the TMB template).
  # Default to 1 (correction on) when absent so a data_list produced here is
  # usable by MakeADFun directly -- the same guarantee comp_offset gives; fit_mod()
  # overrides both from fit_control() before fitting.
  bias_adjust_obs  <- data_list$bias_adjust_obs
  if (is.null(bias_adjust_obs))  bias_adjust_obs  <- 1
  bias_adjust_obs  <- as.numeric(bias_adjust_obs)[1]
  bias_adjust_proc <- data_list$bias_adjust_proc
  if (is.null(bias_adjust_proc)) bias_adjust_proc <- 1
  bias_adjust_proc <- as.numeric(bias_adjust_proc)[1]

  endyr    <- data_list$endyr
  flt_type <- data_list$flt_type   # indexed by Fleet_code, matching the TMB cpp
  nages    <- data_list$nages
  nlengths <- data_list$nlengths

  obsvec_parts <- list()           # numeric blocks; concatenated once at the end
  ctl_parts    <- list()           # obs_ctl row blocks; rbind'd once at the end
  next_pos     <- 0L               # 0-based position counter (TMB convention)

  # Append a block of observations to `obsvec`/`obs_ctl` and return their
  # 0-based positions. `value` is what the likelihood reads. When
  # `one_group = TRUE` all elements share a single group id (one decomposition
  # unit, e.g. the bins of one composition); otherwise each element is its own
  # group (e.g. one aggregate observation). Blocks are accumulated via `<<-` and
  # joined once after the loops (see the assembly step at the end).
  append_obs <- function(value, source, data_row, fleet_code, species, year,
                         sex = NA_integer_, age_length_bin = NA_integer_,
                         length = NA_integer_, bin_index = NA_integer_,
                         comp_type = NA_integer_, is_last_bin = FALSE,
                         stomach_id = NA_integer_, one_group = FALSE,
                         accumulated = FALSE) {
    n <- length(value)
    if (n == 0) return(integer(0))
    pos      <- next_pos + seq_len(n) - 1L
    group_id <- if (one_group) rep(pos[1], n) else pos
    k <- length(obsvec_parts) + 1L
    obsvec_parts[[k]] <<- as.numeric(value)
    ctl_parts[[k]]    <<- data.frame(
      obs_pos        = as.integer(pos),
      source         = source,
      data_row       = as.integer(data_row),
      fleet_code     = as.integer(fleet_code),
      species        = as.integer(species),
      sex            = as.integer(sex),
      year           = as.integer(year),
      age_length_bin = as.integer(age_length_bin),
      length         = as.integer(length),
      bin_index     = as.integer(bin_index),
      comp_type     = as.integer(comp_type),
      is_last_bin   = as.logical(is_last_bin),
      accumulated   = as.logical(accumulated),
      stomach_id    = as.integer(stomach_id),
      group_id      = as.integer(group_id),
      stringsAsFactors = FALSE)
    next_pos <<- next_pos + n
    pos
  }

  # Append one composition row (comp or caal): its bin counts =
  # (proportion + comp_offset) * Neff, one decomposition group, final bin
  # flagged. Returns the obsvec start position of the row's first bin.
  # `bin_label` is the age/length ordinal each element of `obs_row` stands for.
  # It matters once tail accumulation has folded the row: the folded vector is
  # shorter than the fleet's bin dimension, so numbering it 1..n_bins labels
  # each residual with a bin it does not represent. Defaults to 1..n_bins,
  # which is what an unfolded row means.
  append_composition <- function(source, obs_row, n_bins, Neff, fleet, sp, sex,
                                 yr, data_row, comp_type = NA_integer_,
                                 length = NA_integer_, bin_label = NULL,
                                 accumulated = NULL) {
    counts <- (as.numeric(obs_row[seq_len(n_bins)]) + comp_offset) * Neff
    if (is.null(bin_label))   bin_label   <- seq_len(n_bins)
    if (is.null(accumulated)) accumulated <- rep(FALSE, n_bins)
    append_obs(value = counts, source = source, data_row = data_row,
               fleet_code = fleet, species = sp, sex = sex, year = yr,
               age_length_bin = bin_label, length = length,
               bin_index = seq_len(n_bins) - 1L, comp_type = comp_type,
               is_last_bin = seq_len(n_bins) == n_bins, one_group = TRUE,
               accumulated = accumulated)[1]
  }

  # ---- Index (survey) observations ----
  # TMB guard: Year in (0, endyr], fleet on (flt_type > 0), observation > 0.
  # The value laid into obsvec depends on the index likelihood family
  # (index_ll_type), matching what the cpp reads in OSA mode so oneStepPredict()
  # residualizes the same model that was fit:
  #   0 lognormal IID  -> log(obs)                     (dnorm on the log scale)
  #   3 Normal, 4 TruncatedNormal -> obs                (dnorm on the natural
  #       scale; family 4 adds the truncation constant in the cpp, which shifts
  #       the density but not the value stored here)
  #   1/2 MVN / MVNORM -> z = L^-1 obs, whitened by the
  #       lower Cholesky of the fleet's covariance Sigma = L L', so the correlated
  #       block becomes independent standard normals (the cpp whitens the mean with
  #       the same L). The innovation z - L^-1(q*pred) is the multivariate-Gaussian
  #       one-step-ahead residual (Thygesen et al. 2017; the SAM/TMB construction).
  # The per-family scale is only needed for the OSA build; the ordinary fit reads
  # obsvec only for family 0 (and reads index_obs directly for 1/2/3/4), so the
  # fast path (build_osa == FALSE) keeps the original log(obs) layout for every
  # fleet.
  index_ctl <- data_list$index_ctl
  index_obs <- data_list$index_obs
  index_obsvec_idx <- rep(-1L, nrow(index_obs))
  if (nrow(index_obs) > 0) {
    inc <- which(index_ctl[, 3] > 0 & index_ctl[, 3] <= endyr &
                   flt_type[index_ctl[, 1]] > 0 & index_obs[, 1] > 0)
    ill <- data_list$index_ll_type
    if (build_osa && !is.null(ill) && length(inc) > 0) {
      # MVN / MVNORM (families 1, 2): whiten each fleet's fitted observation block
      # with the lower Cholesky of its covariance. Rows are taken in ascending
      # index_obs order, matching Sigma's row order (.align_index_cov() builds
      # Sigma in this order) and the cpp residual assembly. The lower-triangular
      # whitening conditions observation k on rows 1..k, so this row order is also
      # the one-step-ahead conditioning order; CEATTLE index_obs is chronological,
      # making it the chronological order. A fleet whose fitted years are NOT in
      # ascending order is excluded (rather than given a non-chronological, and
      # potentially Sigma-misaligned, decomposition).
      for (f in sort(unique(index_ctl[inc, 1][ill[index_ctl[inc, 1]] %in% c(1L, 2L)]))) {
        rows <- inc[index_ctl[inc, 1] == f]
        Sigma <- tryCatch(as.matrix(data_list$index_cov_mat[[f]]), error = function(e) NULL)
        chrono <- !is.unsorted(index_ctl[rows, 3], strictly = FALSE)
        L <- if (!is.null(Sigma) && nrow(Sigma) == length(rows) && chrono)
          tryCatch(t(chol(Sigma)), error = function(e) NULL) else NULL
        if (is.null(L)) {
          # Malformed / non-PD / mis-dimensioned covariance, or non-chronological
          # rows: fall back to excluding this fleet from the OSA residuals rather
          # than emit a wrong or ambiguously-ordered residual.
          # A fleet at obsvec position -1 adds nothing to the residualization
          # model, so under random effects excluding it moves every other residual.
          warning(sprintf(paste0(
            "OSA residuals: index fleet %d has survey rows out of year order or a ",
            "missing / non-positive-definite / non-%dx%d covariance matrix. Its ",
            "survey likelihood is dropped from the residual model, so it gets no ",
            "residuals and, under random effects, the other fleets' residuals ",
            "shift."), f, length(rows), length(rows)))
        } else {
          z <- as.numeric(forwardsolve(L, index_obs[rows, 1]))   # L^-1 obs (whitened)
          index_obsvec_idx[rows] <- append_obs(
            value = z, source = "index", data_row = rows,
            fleet_code = index_ctl[rows, 1], species = index_ctl[rows, 2],
            year = index_ctl[rows, 3])
        }
        inc <- setdiff(inc, rows)
      }
      # Natural-scale families -- "Normal" (3) and "TruncatedNormal" (4) -- store
      # the UNTRANSFORMED observation, because that is what the cpp's OSA branch
      # scores it against. Both belong here: family 4 differs only by the
      # truncation constant, which is a function of the predicted index and not
      # of the observation, so it changes the density but not the value stored.
      # A natural-scale family that misses this list falls through to the
      # lognormal default below and is residualized as log(obs) against a
      # natural-scale mean.
      rows_nat <- inc[ill[index_ctl[inc, 1]] %in% c(3L, 4L)]
      if (length(rows_nat) > 0) {
        index_obsvec_idx[rows_nat] <- append_obs(
          value = index_obs[rows_nat, 1], source = "index", data_row = rows_nat,
          fleet_code = index_ctl[rows_nat, 1], species = index_ctl[rows_nat, 2],
          year = index_ctl[rows_nat, 3])
        inc <- setdiff(inc, rows_nat)
      }
    }
    # Lognormal IID (family 0) -- and every fleet on the fast fitting path -- as
    # log(obs).
    if (length(inc) > 0) {
      index_obsvec_idx[inc] <- append_obs(
        value = log(index_obs[inc, 1]), source = "index", data_row = inc,
        fleet_code = index_ctl[inc, 1], species = index_ctl[inc, 2],
        year = index_ctl[inc, 3])
    }
  }

  # ---- Catch (fishery) observations: lognormal, stored as log(obs) ----
  # TMB guard: Year in (0, endyr], fishery fleet (flt_type == 1), catch > 0.
  catch_ctl <- data_list$catch_ctl
  catch_obs <- data_list$catch_obs
  catch_obsvec_idx <- rep(-1L, nrow(catch_obs))
  if (nrow(catch_obs) > 0) {
    inc <- which(catch_ctl[, 3] > 0 & catch_ctl[, 3] <= endyr &
                   flt_type[catch_ctl[, 1]] == 1 & catch_obs[, 1] > 0)
    if (length(inc) > 0) {
      catch_obsvec_idx[inc] <- append_obs(
        value = log(catch_obs[inc, 1]), source = "catch", data_row = inc,
        fleet_code = catch_ctl[inc, 1], species = catch_ctl[inc, 2],
        year = catch_ctl[inc, 3])
    }
  }

  # ---- Environmental-covariate (Rogers QAR1) observations ----
  # An observed covariate (observe=) is a Gaussian measurement of the latent
  # state, dnorm(obs, latent_re, obs_sd) per observed year (ceattle.cpp row
  # JNLL_LINKAGE_RE). Lay each observed slot into obsvec so oneStepPredict() can
  # residualize it. Slots are one per year over the model span, so slot g maps to
  # year styr + g - 1; unobserved years (obs_mask == 0) get no obsvec entry.
  # Indexing keys off the absolute slot, so interior gaps do not shift later
  # years. slot_year is a display label only (observe= grouping is annual).
  n_re <- length(data_list$linkage_re_obs_value)
  linkage_re_obsvec_idx <- rep(-1L, n_re)
  if (n_re > 0 && !is.null(data_list$linkage_re_obs_mask)) {
    obs_mask  <- data_list$linkage_re_obs_mask
    slot_year <- data_list$styr + seq_len(n_re) - 1L
    for (g in which(obs_mask == 1L)) {
      linkage_re_obsvec_idx[g] <- append_obs(
        value = data_list$linkage_re_obs_value[g], source = "ecov",
        data_row = g, fleet_code = NA_integer_, species = 1L,
        year = slot_year[g])
    }
  }

  # ---- Age/length composition: bin counts (proportion + 1e-5) * N ----
  # TMB guard: Year in (0, endyr], fleet on, Sample_size > 0. Bins per row:
  # nages (age comp) or nlengths (length comp), doubled for joint-sex (Sex==3).
  # One decomposition unit per row; the last bin is fixed by the sum-to-N
  # constraint (is_last_bin -> residual NA).
  comp_ctl <- data_list$comp_ctl
  comp_obs <- data_list$comp_obs
  comp_n   <- data_list$comp_n
  comp_obsvec_idx <- rep(-1L, nrow(comp_obs))
  if (build_osa && nrow(comp_obs) > 0) {
    accum_yng <- data_list$comp_accum_young; accum_old <- data_list$comp_accum_old
    for (r in seq_len(nrow(comp_obs))) {
      fleet     <- comp_ctl[r, 1]; sp <- comp_ctl[r, 2]; sex <- comp_ctl[r, 3]
      comp_type <- comp_ctl[r, 4]; yr <- comp_ctl[r, 5]; Neff <- comp_n[r, 2]
      if (!(yr > 0 && yr <= endyr && flt_type[fleet] > 0 && Neff > 0)) next
      joint_adjust <- if (sex == 3) 2L else 1L          # joint-sex doubles the bins
      nbins_blk <- if (comp_type == 0) nages[sp] else nlengths[sp]
      # Composition tail accumulation (Comp_accum_young/old): fold the RAW
      # proportion row exactly as the fitting path does (ceattle.cpp Slot 2)
      # BEFORE append_composition() adds the offset, so the OSA obsvec matches the
      # folded comp_hat_tmp / n_comp the un-gated cpp fold now produces on the OSA
      # path too. yng/old clamp mirrors the cpp; the default (yng == 1, old >=
      # nbins_blk) is a no-op, leaving non-accumulating fleets bit-identical.
      obs_row <- comp_obs[r, ]
      n_comp  <- nbins_blk * joint_adjust
      yng <- if (!is.null(accum_yng)) accum_yng[fleet] else 1L
      old <- if (!is.null(accum_old)) accum_old[fleet] else 0L
      if (is.na(yng) || yng < 1L) yng <- 1L                 # NA/0 -> no young accum
      if (is.na(old) || old < 1L || old > nbins_blk) old <- nbins_blk  # NA/0/>=nbins -> no old accum
      if (yng > old) yng <- old
      # Label by the fleet's own bin dimension, offset by sex block. Unfolded
      # rows give `block*nbins + 1..nbins`, the previous numbering.
      blk_off <- rep((seq_len(joint_adjust) - 1L) * nbins_blk,
                     each = old - yng + 1L)
      bin_label <- blk_off + rep(seq.int(yng, old), times = joint_adjust)
      # The boundary bins carry the folded tails, so they are not comparable with
      # an unaccumulated bin of the same ordinal; flag them rather than let a
      # reader assume `age_length_bin` means that age alone.
      acc_one <- c(yng > 1L, rep(FALSE, max(0L, old - yng - 1L)),
                   if (old > yng) old < nbins_blk else logical(0))
      accumulated <- rep(acc_one, times = joint_adjust)
      if (yng > 1L || old < nbins_blk) {
        obs_row <- .fold_comp_bins(as.numeric(obs_row[seq_len(n_comp)]),
                                   nbins_blk, joint_adjust, yng, old)
        n_comp  <- joint_adjust * (old - yng + 1L)
      }
      comp_obsvec_idx[r] <- append_composition(
        "comp", obs_row, n_comp, Neff, fleet, sp, sex, yr, r,
        comp_type = comp_type, bin_label = bin_label,
        accumulated = accumulated)
    }
  }

  # ---- Conditional age-at-length (caal): age-bin counts per length bin ----
  # TMB guard: Year in (0, endyr], fleet on, Sample_size > 0. n_caal = nages.
  caal_ctl <- data_list$caal_ctl
  caal_obs <- data_list$caal_obs
  caal_n   <- data_list$caal_n
  caal_obsvec_idx <- rep(-1L, nrow(caal_obs))
  if (build_osa && nrow(caal_obs) > 0) {
    for (r in seq_len(nrow(caal_obs))) {
      fleet <- caal_ctl[r, 1]; sp <- caal_ctl[r, 2]; sex <- caal_ctl[r, 3]
      yr    <- caal_ctl[r, 4]; len_bin <- caal_ctl[r, 5]; Neff <- caal_n[r, 1]
      if (!(yr > 0 && yr <= endyr && flt_type[fleet] > 0 && Neff > 0)) next
      caal_obsvec_idx[r] <- append_composition(
        "caal", caal_obs[r, ], nages[sp], Neff, fleet, sp, sex, yr, r,
        comp_type = 0L, length = len_bin)            # age bins conditioned on length
    }
  }

  # ---- Diet composition: per-stomach prey counts + "other prey" ----
  # TMB guard: the predator's suitability is estimated (suitMode[pred] > 0) and
  # the stomach has >= 1 prey item. The composition for a stomach is its prey
  # items plus an "other prey" category (the last bin, fixed by the sum-to-1
  # constraint and dropped). One decomposition unit per stomach; diet_obsvec_idx
  # is indexed by 0-based stomach id, so it has length n_stomach_obs.
  n_stomach <- data_list$n_stomach_obs
  diet_obsvec_idx <- if (is.null(n_stomach)) integer(0) else rep(-1L, n_stomach)
  if (build_osa && !is.null(n_stomach) && n_stomach > 0 &&
      !is.null(data_list$diet_ctl) && nrow(data_list$diet_ctl) > 0) {
    diet_ctl   <- data_list$diet_ctl
    diet_obs   <- data_list$diet_obs
    stomach_id <- data_list$stomach_id
    suitMode   <- data_list$suitMode
    for (i in seq_len(n_stomach) - 1L) {        # stomach ids are 0-based
      # which() is order-independent where the C++ diet likelihood instead takes
      # a contiguous run of stomach_id (ceattle.cpp, section 13.2). The two agree
      # because stomach_id is numbered 0, 1, ... with each stomach's rows
      # consecutive -- established by clean_data() and enforced by data_check().
      # Were that to lapse, this would still find every prey row while the C++
      # scan would not, and the diet OSA "other prey" bin would misalign.
      rows   <- which(stomach_id == i)
      n_prey <- length(rows)
      if (n_prey == 0) next
      rsp <- diet_ctl[rows[1], 1]               # predator species (1-based)
      if (is.null(suitMode) || suitMode[rsp] <= 0) next
      N_s  <- diet_obs[rows[1], 1]              # sample size (number of stomachs)
      # observed prey proportions plus the residual "other prey" category, then
      # the same offset/normalize the TMB likelihood applies, scaled to counts.
      obs_p  <- as.numeric(diet_obs[rows, 2])
      vec    <- c(obs_p, 1 - min(sum(obs_p), 1)) + comp_offset
      counts <- vec / sum(vec) * N_s
      diet_obsvec_idx[i + 1L] <- append_obs(
        value          = counts, source = "diet",
        data_row       = c(rows, NA_integer_),
        fleet_code     = rsp, species = rsp,
        year           = diet_ctl[rows[1], 7],
        age_length_bin = c(diet_ctl[rows, 2], NA_integer_),   # prey id; NA = other
        bin_index      = seq_len(n_prey + 1L) - 1L,
        is_last_bin   = seq_len(n_prey + 1L) == (n_prey + 1L),
        stomach_id    = i, one_group = TRUE)[1]
    }
  }

  # Assemble the flat vector and metadata in a single pass from the blocks
  # accumulated above, instead of growing them with c() / rbind() inside the
  # loops (which is quadratic and runs on every fit_mod() call).
  obsvec  <- if (length(obsvec_parts)) unlist(obsvec_parts, use.names = FALSE) else numeric(0)
  obs_ctl <- if (length(ctl_parts)) do.call(rbind, ctl_parts) else .new_obs_ctl()
  rownames(obs_ctl) <- NULL

  # TMB needs a non-empty DATA_VECTOR; use a length-1 sentinel when no
  # observations are included (the *_obsvec_idx vectors are then all -1 and the
  # sentinel value is never read by the model).
  if (length(obsvec) == 0) obsvec <- 0

  data_list$obsvec           <- as.numeric(obsvec)
  data_list$obs_ctl          <- obs_ctl
  data_list$index_obsvec_idx <- as.integer(index_obsvec_idx)
  data_list$catch_obsvec_idx <- as.integer(catch_obsvec_idx)
  data_list$comp_obsvec_idx  <- as.integer(comp_obsvec_idx)
  data_list$caal_obsvec_idx  <- as.integer(caal_obsvec_idx)
  data_list$diet_obsvec_idx  <- as.integer(diet_obsvec_idx)
  data_list$linkage_re_obsvec_idx <- as.integer(linkage_re_obsvec_idx)
  data_list$osa_mode         <- 0L
  # Read by the TMB DATA_INTEGER. Defaulted here rather than in fit_mod() so a
  # data list built before this option existed, or rebuilt by .osa_build_obj(),
  # still carries the field MakeADFun() requires.
  if (is.null(data_list$adreport_sel)) data_list$adreport_sel <- 0L
  data_list$comp_offset      <- comp_offset       # read by the TMB DATA_SCALAR
  data_list$bias_adjust_obs  <- bias_adjust_obs   # read by the TMB DATA_SCALAR
  data_list$bias_adjust_proc <- bias_adjust_proc  # read by the TMB DATA_SCALAR

  # Simulation switches, read by the TMB DATA_IVECTORs. Set here because this is
  # the single funnel every data list passes through on its way to MakeADFun
  # (rearrange_data() for a fit, .osa_build_obj() for a rebuild). Process error
  # defaults OFF: redrawing a process changes what a self-test measures, so it
  # is asked for rather than assumed. sim_mod() overrides them per call.
  #   simulate_state : 0 recruitment (annual and initial), 1 M, 2 growth,
  #                    3 catchability, 4 selectivity -- the linkage process
  #                    codes, so the model can index by them directly
  #   simulate_period: 0 the fitted window, 1 outside it. Only slot 0 is read
  #                    by the model, and only by the PROCESS draws; the
  #                    observation draws are ungated by period on purpose (see
  #                    ceattle.cpp 2.4.2c). Slot 1 is reserved and inert.
  if (is.null(data_list$simulate_state))  data_list$simulate_state  <- rep(0L, 5)
  if (is.null(data_list$simulate_period)) data_list$simulate_period <- rep(1L, 2)
  data_list$simulate_state  <- as.integer(data_list$simulate_state)
  data_list$simulate_period <- as.integer(data_list$simulate_period)

  data_list
}


#' Empty `obs_ctl` metadata frame with the correct column types
#'
#' @return A 0-row data frame with the columns used to map each `obsvec` element
#'   back to its observation. See [build_osa_data()].
#' @keywords internal
.new_obs_ctl <- function() {
  data.frame(
    obs_pos        = integer(0),
    source         = character(0),
    data_row       = integer(0),
    fleet_code     = integer(0),
    species        = integer(0),
    sex            = integer(0),
    year           = integer(0),
    age_length_bin = integer(0),
    length         = integer(0),
    bin_index     = integer(0),
    comp_type     = integer(0),
    is_last_bin   = logical(0),
    stomach_id    = integer(0),
    group_id      = integer(0),
    stringsAsFactors = FALSE)
}
