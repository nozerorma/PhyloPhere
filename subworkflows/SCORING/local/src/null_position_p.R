# Position-level empirical p and BH over the fixed tested-position family.
# p.emp = (k_emp + 1)/(N + 1), inclusive upper tail. Null-only positions enter BH at p=1.
# Legacy fit/c_k inputs are accepted but ignored; no pooled conditional-tail estimator is produced.

# obs_max   data frame (Gene, Position, .obs): the observed positive positions and their best-side score
# counters  data frame (Gene, Position, k_emp, nd): one row per observed position
# fit       ignored legacy argument, retained for cached callers
# n_u       size of U among the null-positive positions (distinct positions)
# N         cycles of the null (the roster)
# Returns list(pos = data frame (Gene, Position, p.emp, p.adj_bh, nd, null_unseen),
#              n_matched, n_obs, n_unseen, n_family)
null_position_p <- function(obs_max, counters, fit, n_u, N, TIE_TOL=1e-12, P_EMP_GUARD_MIN_POS=10L,
                            n_family=NULL) {
  keys <- c('Gene','Position')
  if(!all(c(keys,'.obs') %in% names(obs_max)) || !all(c(keys,'k_emp','nd') %in% names(counters)) ||
     anyDuplicated(obs_max[keys]) || anyDuplicated(counters[keys])) stop('Invalid/duplicate position counter keys')
  if(length(N)!=1 || !is.finite(N) || N<1 || N!=floor(N)) stop('Invalid null cycle count')
  key <- function(d) paste(d$Gene,d$Position,sep='\r')
  at <- match(key(obs_max),key(counters))
  if(anyNA(at) || nrow(counters)!=nrow(obs_max)) stop('Null counters and observed positions disagree')
  ctr <- counters[at,,drop=FALSE]
  for(k in c('k_emp','nd')) ctr[[k]] <- suppressWarnings(as.numeric(ctr[[k]]))
  for(k in c('k_emp','nd')) if(anyNA(ctr[[k]]) || any(!is.finite(ctr[[k]]) | ctr[[k]]<0 | ctr[[k]]!=floor(ctr[[k]]))) stop('Invalid null counters')
  if(any(ctr$k_emp>ctr$nd | ctr$nd>N)) stop('Null counters violate 0 <= k_emp <= nd <= N')
  if(any(!is.finite(obs_max$.obs) | obs_max$.obs<0 | obs_max$.obs>1)) stop('Invalid observed scores')
  nd <- ctr$nd
  p <- ifelse(obs_max$.obs<=0,1,(ctr$k_emp+1)/(N+1))
  # Compatibility for explicit helper callers. Production always supplies the fixed tested family.
  if(is.null(n_family)) n_family <- n_u+sum(nd==0)
  if(length(n_family)!=1 || !is.finite(n_family) || n_family<nrow(obs_max) || n_family!=floor(n_family)) stop('Invalid fixed BH family size')
  q <- p.adjust(p,'BH',n=n_family)
  list(pos=data.frame(Gene=obs_max$Gene,Position=obs_max$Position,p.emp=p,nd=nd,
                      null_unseen=nd==0,p.adj_bh=q,stringsAsFactors=FALSE),
       n_matched=sum(nd>0),n_obs=nrow(obs_max),n_unseen=sum(nd==0),n_family=n_family)
}

# Rows of null_universe_summary.tsv (metric, cycles, value, note): N, |U|, observed positive positions, the ones the null
# never scored, and the saturation curve of the null. `base` holds what the null itself fixes (null_universe_base.tsv of
# null_position_stats.py); an unseen fraction is a diagnostic of the adequacy of the null (does N reach the positions that
# can be scored?), not a reason to demote a position.
null_universe_summary <- function(base, n_obs, n_unseen) {
  val <- function(m) base$value[base$metric == m & is.na(base$cycles)][1]
  note <- function(m) { v <- base$note[base$metric == m & is.na(base$cycles)][1]; if (is.na(v)) "" else v }
  n_u <- val("null_positive_positions_in_universe")
  metrics <- c("N_cycles", "N_cycles_present_in_null", "null_positive_positions_all", "null_positive_positions_in_universe",
               "observed_positive_positions", "observed_unseen_positions", "fraction_unseen", "universe_size")
  curve <- base[base$metric == "null_positive_positions_cumulative", c("metric", "cycles", "value", "note")]
  # the gene list the universe received: value = number of genes, note = digest (null_position_stats.py)
  genes <- base[base$metric == "gene_list_sha256" & is.na(base$cycles), c("metric", "cycles", "value", "note")]
  rbind(
    data.frame(metric = metrics, cycles = NA_integer_,
               value = c(val("N_cycles"), val("N_cycles_present_in_null"), val("null_positive_positions_all"), n_u, n_obs,
                         n_unseen, if (n_obs > 0) n_unseen / n_obs else NA_real_, if ('tested_positions_family' %in% base$metric) val('tested_positions_family') else n_u + n_unseen),
               note = c(note("N_cycles"), "", "", note("null_positive_positions_in_universe"), "",
                        "positions whose null detections are 0 (p = 1/(N+1), kept)", "",
                        "fixed tested-position BH family (legacy fallback explicitly identified by missing tested_positions_family)"),
               stringsAsFactors = FALSE),
    genes, curve)
}
