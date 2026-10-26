# Encode an Rceattle linkage table into TMB-friendly inputs

The template reads the linkage table as parallel integer and numeric
vectors plus a design matrix. Text columns (`process`, `param`, `link`,
`prior_family`) become 0-based integer codes, an `NA` stratum id
(`species`, `sex`, `age_bin`, `fleet`) becomes `0` ("every level"), and
`X_col` becomes 0-based.
