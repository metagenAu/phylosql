#' Convert master sheet metadata to CMS format
#'
#' @param data A data frame containing at least the CMS fields listed below.
#'
#' @return A tibble ready for upload to the `cmsdata` table.
#' @export
ms_as_cmsData <- function(data) {
  required <- c(
    "MetagenNumber", "PropertyName", "BlockName", "CropName",
    "SurveyDate", "Location", "AgronomistName"
  )
  cms <- ensure_dataframe(data, required, "master sheet metadata")

  if (!"TargetCropTypeName" %in% names(cms)) {
    cms$TargetCropTypeName <- NA_character_
  }

  tibble::as_tibble(cms) %>%
    dplyr::mutate(
      SurveyID = NA_character_,
      TargetCropTypeName = dplyr::coalesce(rlang::.data$TargetCropTypeName, NA_character_)
    ) %>%
    dplyr::select(
      "MetagenNumber", "PropertyName", "SurveyID", "BlockName",
      "CropName", "SurveyDate", "Location", "AgronomistName",
      "TargetCropTypeName"
    )
}

#' Normalise CMS exports for upload
#'
#' @param data CMS data where `SurveyId` and `BarcodeId` may require
#'   normalisation.
#'
#' @return A tibble ready for upload to the `cmsdata` table.
#' @export
cms_as_cmsData <- function(data) {
  cms <- ensure_dataframe(data, NULL, "CMS metadata")

  if (!"TargetCropTypeName" %in% names(cms)) {
    cms$TargetCropTypeName <- NA_character_
  }

  cms %>%
    tibble::as_tibble() %>%
    dplyr::rename(
      SurveyID = "SurveyId",
      MetagenNumber = "BarcodeId"
    ) %>%
    dplyr::mutate(SurveyID = as.character(rlang::.data$SurveyID)) %>%
    dplyr::select(
      "MetagenNumber", "PropertyName", "SurveyID", "BlockName",
      "CropName", "SurveyDate", "Location", "AgronomistName",
      "TargetCropTypeName"
    )
}

#' Tidy laboratory measurements for upload
#'
#' @param labdata Data frame of laboratory metrics.
#'
#' @return A tibble in long format with columns `MetagenNumber`, `variable`, and
#'   `value`.
#' @export
as_labData <- function(labdata) {
  required <- c(
    "MetagenNumber", "pH", "ActiveCarbon", "AggregateStability",
    "Phosphatase", "B_Glucosidase", "DNA.conc...ng.ul.", "SoilMoisture"
  )

  lab <- ensure_dataframe(labdata, required, "laboratory data") %>%
    tibble::as_tibble()

  renamed <- dplyr::rename(lab, DNAConc = "DNA.conc...ng.ul.")

  cleaned <- renamed %>%
    dplyr::mutate(
      SoilMoisture = suppressWarnings(as.numeric(rlang::.data$SoilMoisture)),
      dplyr::across(
        c("ActiveCarbon", "pH", "Phosphatase", "B_Glucosidase", "DNAConc", "SoilMoisture"),
        ~ replace(.x, .x == 0, NA_real_)
      ),
      SoilMoisture = replace(rlang::.data$SoilMoisture, rlang::.data$SoilMoisture < 0, NA_real_)
    )

  cleaned %>%
    tidyr::pivot_longer(
      cols = -"MetagenNumber",
      names_to = "variable",
      values_to = "value",
      values_drop_na = TRUE
    )
}

#' Derive soil moisture ratios from wet/dry masses
#'
#' @param labdata Laboratory measurements including `SoilDW` and `SoilWW`.
#'
#' @return The input data with an updated `SoilMoisture` column.
#' @export
fix_soil_moisture <- function(labdata) {
  data <- ensure_dataframe(labdata, c("SoilDW", "SoilWW", "SoilMoisture"), "laboratory data")

  ratio <- with(data, SoilDW / SoilWW)
  data$SoilMoisture[!is.na(ratio)] <- ratio[!is.na(ratio)]
  data
}

#' Convert CMS data to long format
#'
#' @param data CMS data in wide format.
#'
#' @return A tibble with columns `MetagenNumber`, `Factor`, and `Level`.
#' @export
cms_as_longCms <- function(data) {
  cms <- ensure_dataframe(data, c("MetagenNumber", "PropertyName"), "CMS metadata")

  cms %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      cols = -"MetagenNumber",
      names_to = "Factor",
      values_to = "Level",
      values_drop_na = TRUE
    )
}

#' Convert miscellaneous wide data to CMS long format
#'
#' @param data A data frame whose first column is `MetagenNumber`.
#'
#' @return A tibble with columns `MetagenNumber`, `Factor`, and `Level`.
#' @export
misc_as_longCms <- function(data) {
  if (!ncol(data)) {
    stop("`data` must contain at least one column.", call. = FALSE)
  }

  cms <- ensure_dataframe(data, colnames(data)[1], "miscellaneous CMS metadata")
  names(cms)[1] <- "MetagenNumber"

  cms %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(
      cols = -"MetagenNumber",
      names_to = "Factor",
      values_to = "Level",
      values_drop_na = TRUE
    )
}

#' Categorise qualitative DNA yield measurements
#'
#' @param dnas Data frame containing a `DNA` column with qualitative values.
#'
#' @return The input data with an added `range` column categorising the result
#'   as "above", "below", or "inrange".
#' @export
clean_dna_yields <- function(dnas) {
  data <- ensure_dataframe(dnas, c("DNA"), "DNA yield records")
  dna_values <- as.character(data$DNA)

  above <- grepl(">", dna_values, fixed = TRUE)
  below <- grepl("<", dna_values, fixed = TRUE)

  numeric <- gsub("[<>]", "", dna_values)

  data$DNA <- numeric
  data$range <- "inrange"
  data$range[above] <- "above"
  data$range[below] <- "below"

  data
}

#' Convert DNA yield to kilograms
#'
#' @param dna Numeric vector of DNA concentrations (ng/µL).
#' @param dna_elution_volume Elution volume in microlitres. Defaults to 2000.
#'
#' @return Numeric vector of yields in kilograms.
#' @export
transform_dna_yield_to_kg <- function(dna, dna_elution_volume = 2000) {
  dna * dna_elution_volume * 3.33 * 1.55 * 10 * 1e-6
}
