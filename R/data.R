#' A dataset from a landmark colon cancer trial
#'
#' @description A landmark colon cancer trial on the efficacy of levamisole and fluorouracil
#' was reported by Moertel et al. (1990). The trial recruited 929 patients with stage C disease
#' and randomly assigned them to levamisole treatment alone, levamisole combined with fluorouracil,
#' and the control. The dataset here is restricted to the comparison between the combined treatment
#' and control groups, consisting of 304 and 314 patients, respectively.
#'
#'
#'
#' @format A data frame with 915 rows and 6 variables:
#' \describe{
#'   \item{id}{Unique patient ID.}
#'   \item{time}{Event time (years).}
#'   \item{status}{Event type; 1 = cancer relapse, 2 = death.}
#'   \item{rx}{"Lev+5FU" = combined treatment, "Control" = control.}
#'   \item{age}{Patient age (years) at randomization.}
#'   \item{sex}{0 = female, 1 = male.}
#'   }
#' @references MOERTEL, C. G., FLEMING, T. R., MACDONALD, J. S., HALLER, D. G., LAURIE, J. A., GOODMAN, P. J.,
#' UNGERLEIDER, J. S., EMERSON, W. A., TORMEY, D. C., GLICK, J. H. et al. (1990). Levamisole and fluorouracil
#'  for adjuvant therapy of resected colon carcinoma. New Engl. J. Med. 322, 352--358.
"colon_lev"


#' A dataset from the HF-ACTION trial
#'
#' @description Over two thousand heart failure patients across the USA, Canada, and France
#' participated in the Heart Failure: A Controlled Trial Investigating Outcomes of Exercise
#' Training (HF-ACTION) between 2003--2007 (O'Connor et al., 2009).
#' The primary objective of the trial was to evaluate the effect of adding exercise
#' training to the usual patient care on the composite endpoint of all-cause hospitalization
#'  and death.
#' The dataset here contains a subgroup of 426 non-ischemic patients
#' with baseline cardio-pulmonary exercise test less than or equal to nine minutes.
#'
#'
#'
#' @format A data frame with 1,448 rows and 5 variables:
#' \describe{
#'   \item{patid}{Unique patient ID.}
#'   \item{time}{Event time (years).}
#'   \item{status}{Event type; 1 = hospitalization, 2 = death.}
#'   \item{trt_ab}{1 = exercise training, 0 = usual care.}
#'   \item{age60}{1 = 60 years or older, 0 = otherwise.}
#'   }
#' @references O'CONNOR, C. M., WHELLAN, D. J., LEE, K. L., KETEYIAN, S. J.,
#'  COOPER, L. S., ELLIS, S. J., LEIFER, E. S.,
#' KRAUS, W. E., KITZMAN, D. W., BLUMENTHAL, J. A. et al. (2009).
#' Efficacy and safety of exercise training in
#'  patients with chronic heart failure: Hf-action randomized controlled trial.
#'  J. Am. Med. Assoc. 301, 1439--1450.
"hfaction"
