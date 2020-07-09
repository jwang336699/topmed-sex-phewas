#' Sets up a workspace for PIC-SURE analyses. The working directory must contain a working PIC-SURE
#' token and should also have the basic PIC-SURE processing libraries (available on this GitHub
#' repository: https://github.com/jwang336699/topmed-sex-phewas)
#'
#' @param wd a path to a working directory containing a PIC-SURE token in a file named 'token.txt'
#' and the R_lib directory from the above repository.
#'
#' @return A PIC-SURE resource connection
#' @import picsure hpds
#' @importFrom picsure connect
#' @importFrom hpds get.resource
#' @export
setup <- function(wd) {
  setwd(wd)
  source("R_lib/requirements.R")
  source("R_lib/utils.R")
  PICSURE_network_URL <- "https://picsure.biodatacatalyst.nhlbi.nih.gov/picsure"
  resource_id <- "02e23f52-f354-4e8b-992c-d37c8b9ba140"
  token_file <- "token.txt"
  token <- scan(token_file, what = "character")
  connection <- picsure::connect(url = PICSURE_network_URL,
                                 token = token)
  resource <- hpds::get.resource(connection,
                                 resourceUUID = resource_id)
  return(resource)
}

#' Retrieves the variable dictionaries for a study in the TOPMed database. 
#'
#' @param resource a PIC-SURE resource connection
#' @param study_name a character singleton containing the name of a TOPMed study; see
#' (https://docs.google.com/document/d/1oVmdBSETxHNpB2DIWAh05TH_uUMPeJQKgyFsAGXpQVU/edit#)
#'
#' @return a list containing cdict - a dictionary of consent variables, and vdict - a 
#' dictionary of the variables in the given study
#' @import hpds dplyr
#' @importFrom hpds find.in.dictionary, extract.dataframe
#' @importFrom dplyr bind_rows
#' @export
select_data <- function(resource, study_name) {
  # obtain the study variables
  variables <- hpds::find.in.dictionary(resource, study_name)
  # and the consent variable
  consent_variable <- hpds::find.in.dictionary(resource, "Study Accession with Consent Code")
  
  # find the descriptive dataframes
  study_dict <- hpds::extract.dataframe(variables)
  consent_dict <- hpds::extract.dataframe(consent_variable)
  
  # bind them together
  plain_variablesDict <- dplyr::bind_rows(study_dict, consent_dict) %>% dplyr::arrange(name)
  variablesDict <- get_multiIndex_variablesDict(plain_variablesDict)
  return(list(cdict = consent_dict, vdict = variablesDict))
}

#' Retrieves the study cohort data for a set of variables in the TOPMed database. 
#'
#' @param resource a PIC-SURE resource connection
#' @param vdict a variable dictionary as outputted by \function{select_data} above.
#' @param cdict a specialized consent variable dictionary as outputted by \function{select_data}.
#' @param lvls a list containing vectors of accepted categories for each variable level
#' @param study_id a character singleton containing the study ID as listed here:
#' (https://docs.google.com/document/d/1oVmdBSETxHNpB2DIWAh05TH_uUMPeJQKgyFsAGXpQVU/edit#)
#' @param toggle a binary toggle for matching variables
#'
#' @return a list containing: selected_vars - a binary vector of selected variables, vdict - an
#' enriched dataframe representing a variable dictionary, and dframe - the actual study data
#' @import hpds
#' @importFrom hpds new.query, query.filter.add, query.select.add, query.run
#' @export
retrieve_data <- function(resource, vdict, cdict, lvls = list(NA, NA, NA, NA), study_id, toggle = F) {
  # first select for variables matching the study-level category
  base <- vdict[['level_0']] %in% lvls[[1]]
  
  # match variables based on lower levels
  keep_vars <- rep(toggle, nrow(vdict))
  for (i in 2:4) {
    lvl <- lvls[[i]]
    if (!is.na(lvl)) {
      keep_vars <- keep_vars | (vdict[[paste0('level_', i-1)]] %in% lvl)
    }
  }
  
  # drop some unused variables
  to_drop <- vdict[["simplified_name"]] %in% list("Dbgap_id", "De-identified site code", "A1AD: phenotype/genotype")
  
  # finalize selected variables
  keep_vars <- keep_vars & base
  keep_vars <- keep_vars & !to_drop
  selected_vars <- vdict[keep_vars, ]$name %>% unname() %>% as.list()
  
  # create filter for variables in the desired study
  phs_values = strsplit(cdict[["categoryValues"]], ",") %>% unlist()
  phs_selected = grep(paste0(study_id, "\\.*"), phs_values, perl=T, value=T) %>% as.list()
  
  # build query and execute
  my_query = hpds::new.query(resource = resource)
  hpds::query.filter.add(query = my_query,
                         keys = cdict[["name"]],
                         phs_selected)
  
  hpds::query.select.add(query = my_query,
                         keys = selected_vars)
  final = hpds::query.run(query = my_query, result.type = "dataframe")
  
  # parse names for R formatting (no spaces)
  vdict[["df_name"]] <- parsing_varNames(vdict[["name"]])
  
  return(list(selected_vars = keep_vars, vdict = vdict, dframe = final))
}

summarize_variable <- function(df, vdict, vname) {
  summary(as.factor(df[[vdict[vdict[['simplified_name']] == vname,][['df_name']]]]))
}
