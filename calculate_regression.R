# Title     : Calculates regression
# Objective : Calculates regression and create new columms with forumlas.
# Created by: Wondim
# Created on: 2/11/2019
#install if necessary
#install.packages("stringr")
list.of.packages <- c("mgsub", "gtools", "caret", "stringr")
#install.packages("rJava")
#install.packages('xml2')
new.packages <- list.of.packages[
! (list.of.packages %in% installed.packages()[, "Package"])
]
if (length(new.packages)) install.packages(new.packages)

#load library
library(stringr)
#library(qdap)
library(mgsub)
library(gtools)
library(caret)

permutate_equation_by_band <- function(input_csv, input_list, function_string, minimum_rsquared){
    #' Uses an equation and tries different combinations of columns in the input_list to
    #' come up with list of equation above the R-squared provided.
    #' @param input_csv_path is the file path containing the first column as
    #' a Y axis and the subsequent columns as X axis.
    #' @param input_list are headers of columns that will be used in the equations
    #' @param function_string is an equation to be tested.
    #' @param minimum_rsquared the minimum r-squared used to save the equation.
    #' If it is above the minimum r-squared, the equation is saved in a result file for review.

    place_holders <- regmatches(function_string, gregexpr("[[:alpha:]]+", function_string))
    place <- unlist(place_holders) # convert list to vector

    permutation <- permutations(n = length(input_list), r = length(place), v = input_list, repeats.allowed = T)
   # print (dim(permutation))

    file_conn <- gsub('csv', 'txt', input_csv)
    if (! file_test("-f", file_conn)) {
        file.create(file_conn)
    }
    for (i in 1 : dim(permutation)[1]) {
        permutation_vec <- as.vector(permutation[i,])

        equation <- mgsub::mgsub(function_string, place, permutation_vec) # replace the placeholders with permetuations

        csv <- read.csv(file = input_csv, header = TRUE, sep = ",", row.names = NULL)

        calculated <- (eval(parse(text = equation)))

        r_squared <- postResample(pred = calculated, obs = csv$SPM)['Rsquared']
        # if (! is.na(r_squared)) {
        #     if (r_squared > 0.5) {
        #         print (paste(r_squared, ' - ', equation))
        #     }
        # }


        if (! is.na(r_squared) && r_squared > minimum_rsquared) {
            correl <- cor(calculated,  csv$SPM)

            write(paste(r_squared, ' - ', correl, ' - ', equation), file_conn, append = TRUE)
            #write(paste(r_squared, ' - ', equation), file_conn, append = TRUE)
        }
    }
}

extract_best_equation = function(filepath) {
  con = file(filepath, "r")
  options(digits = 14)
  equation_list = c()
  r_squared_list = c()
  while (TRUE) {
    line = readLines(con, n = 1)
    if (length(line) == 0) {
      break
    }
    line_list = strsplit(line, '  -  ')
    r_squared = as.numeric(line_list[[1]][1])

    equation = line_list[[1]][2]
    equation_list <- append(equation_list, equation)
    r_squared_list <- append(r_squared_list, r_squared)
  }
  close(con)
  # r_squared_list
  return (equation_list[which(r_squared_list == max(r_squared_list))])
}

extract_numbers = function(eq, position = 0, inc = 0) {
  #matches = unlist(regmatches(eq, gregexpr("[^a-zA-Z0-9_$+(^][0-9]+|[0-9]+([.][0-9]+)", eq)))
  pattern = "(?:[^\\d,]|^)(\\d+(?:(?:,\\d+)*,\\d)?\\.\\d+)"
  matches = unlist(str_extract_all(eq, pattern))
  final_matches = ifelse(is.na(matches),
                         gsub("([\\^])|[:ascii:]]|[[:space:]]", "", matches),
                         matches) # strip if found

  return (final_matches)
}
increment_equation_by_position = function(equation0, position) {
  sequence = seq(-10, 10, by = 0.001)
  #print (paste(1, equation1) )
  #return()

  numbers = extract_numbers(equation0)

  equations = c()
  number_count = (length(numbers))
  for (inc in sequence) {
    #print (paste(1, equation))
    # print (paste(inc, numbers, position))
    numbers = extract_numbers(equation0)
    equation1 = equation0
    to <- sprintf("[VAR%d]", 1:number_count)
    for(i in seq_along(numbers)) equation1 = sub(numbers[i], to[i], equation1, fixed = TRUE)
    #print (paste(1, equation1) )
    numbers[position] = toString(as.numeric(numbers[position])+inc)
    #print (numbers)
    for(i in seq_along(numbers)) equation1 = sub(to[i], numbers[i], equation1, fixed = TRUE)
    #print (paste(3, equation1))
    #print (equation)
    equations <- append(equations, equation1)
    # break()
    # if (position + 1 <= number_count) {
    #   increment_equation_by_position(equation1, position + 1)
    # }
  }

  return (equations)
}

handle_equation_placeholders = function(equation_result_file) {
  equation_list = extract_best_equation(equation_result_file)
  place_holdered_equations = c()

  for (eq in unique(equation_list)) {
    final_matches = extract_numbers(eq)
    updated_eq <-
      mgsub::mgsub (eq, final_matches, replicate(length(final_matches), "%f"))
    place_holdered_equations <-
      append(place_holdered_equations, updated_eq)
  }
  return (place_holdered_equations)
}

increment_equations_calc <- function(input_csv, equation, min_r_squared_value = NULL) {
    str_n <- str_count(equation, "%f")
    sequence = seq(-10, 10, by = 0.001)
    # Get number of placeholders
    # Construct arguments list
    arg_list <- lapply(c(1:str_n), function(x)
      sequence)

    arg_list$fmt <- equation # Named argument
    # Call sprintf with constructed arguments
    input_equations <- do.call("sprintf", arg_list)

    equation_list = c()
    r_squared_list = c()
    file_conn = gsub('.csv', '_optimized.txt', input_csv)

    csv <- read.csv(
        file = input_csv,
        header = TRUE,
        sep = ",",
        row.names = NULL
    )
    for (eq in input_equations) {
      for (i in 1:str_n) {
        equations = increment_equation_by_position(eq, i)
        #print (equations)
        for (eq in equations) {
          calculated <- (eval(parse(text = eq)))

          r_squared <- postResample(pred = calculated, obs = csv$SPM)['Rsquared']
          if (!missing(min_r_squared_value)) {
            # Use range - above value
            if (!is.na(r_squared) &&
                r_squared > min_r_squared_value) {

              correl <- cor(calculated,  csv$SPM)

              write(paste(r_squared, ' - ',correl, ' - ', eq), file_conn, append = TRUE)
            }
            equation_list <- append(equation_list, eq)
            r_squared_list <- append(r_squared_list, r_squared)
          }
          else {
            # Max
            if (!is.na(r_squared)) {
              equation_list <- append(equation_list, eq)

              r_squared_list <- append(r_squared_list, r_squared)
            }
          }
        }
      }
    }

    # print (r_squared_list)
    #print (length(r_squared_list))

    if (missing(min_r_squared_value)) {
      # Use max value
      best_equations = equation_list[which(r_squared_list == max(r_squared_list))]

      for (eq in best_equations) {
        r_squared = r_squared_list[which(equation_list == eq)]
        print (paste(r_squared, ' - ', eq))
        write(paste(r_squared, ' - ', eq), file_conn, append = TRUE)
      }
    }
  }

optimize_algorithm = function(input_csv_path, headers, equations, min_r_squared) {
    #' Uses equations and tries different combinations of columns to
    #' come up with list above the R-squared provided.
    #' @param input_csv_path is the file path containing the first column as
    #' a Y axis and the subsequent columns as X axis.
    #' @param headers are headers of columns that will be used in the equations
    #' @param equations are vector of equations that are already tested to bring a good result.
    #' @param min_r_squared the minimum r-squared used to save the equation.
    #' If above the minimum r-squared, the equation is saved in a result file for review.
    result_path = gsub('csv', 'txt', input_csv_path)
    # loop through all equations
    for (equation in equations) {
      permutate_equation_by_band(input_csv_path, headers, equation, min_r_squared)
    }
    # equations with coefficents are extracted.
    equations_with_placeholder = handle_equation_placeholders(result_path)
    # The coefficents of the equations are tested to come up with the best coefficients.
    for (eq in equations_with_placeholder) {
        increment_equations_calc(input_csv_path, eq)
    }
}
# res = stringr::str_extract_all("((csv$RedEdge_717*csv$Red_667^3.400000)/(csv$NIR_831*csv$Green_557^(344)))+(csv$Blue_482/csv$RedEdge_717)", "[^a-zA-Z0-9_$*^][0-9]+|[0-9]+([.][0-9]+)")
# #gsub("[:ascii:]]|[[:punct:]]|[[:space:]]", "", res)
# matches <- regmatches("((csv$RedEdge_717*csv$Red_667^3.400000)/(csv$NIR_831*csv$Green_557^(344)))+(csv$Blue_482/csv$RedEdge_717)", gregexpr("[^a-zA-Z0-9_$*^][0-9]+|[0-9]+([.][0-9]+)", "((csv$RedEdge_717*csv$Red_667^3.400000)/(csv$NIR_831*csv$Green_557^(344)))+(csv$Blue_482/csv$RedEdge_717)"))
# final_matches = ifelse(is.na(as.double(unlist(matches))), gsub("[:ascii:]]|[[:punct:]]|[[:space:]]", "", unlist(matches)), unlist(matches))
# as.double(final_matches)
# # qu = '/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/hypers_total_rrs.txt'
# # increment_equations_calc('/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/mic_data.csv', eq)
#optimize_algorithm('/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/hypers_total_rrs.csv')


#mic_headers <- c('csv$Blue_482', 'csv$Green_557', 'csv$Red_667', 'csv$RedEdge_717', 'csv$NIR_831')
mic_headers <- c('csv$Green_557', 'csv$Red_667', 'csv$RedEdge_717', 'csv$NIR_831')
modis_headers <- c('csv$band8_412', 'csv$band9_443', 'csv$band10_488', 'csv$band11_531','csv$band12_555', 'csv$band13_667', 'csv$band14_678', 'csv$band15_750', 'csv$band16_867')

hypers_headers <- c("csv$Rrs400", "csv$Rrs401", "csv$Rrs402", "csv$Rrs403", "csv$Rrs404", "csv$Rrs405", "csv$Rrs406",
"csv$Rrs407", "csv$Rrs408", "csv$Rrs409", "csv$Rrs410", "csv$Rrs411", "csv$Rrs412", "csv$Rrs413",
"csv$Rrs414", "csv$Rrs415", "csv$Rrs416", "csv$Rrs417", "csv$Rrs418", "csv$Rrs419", "csv$Rrs420",
"csv$Rrs421", "csv$Rrs422", "csv$Rrs423", "csv$Rrs424", "csv$Rrs425", "csv$Rrs426", "csv$Rrs427",
"csv$Rrs428", "csv$Rrs429", "csv$Rrs430", "csv$Rrs431", "csv$Rrs432", "csv$Rrs433", "csv$Rrs434",
"csv$Rrs435", "csv$Rrs436", "csv$Rrs437", "csv$Rrs438", "csv$Rrs439", "csv$Rrs440", "csv$Rrs441",
"csv$Rrs442", "csv$Rrs443", "csv$Rrs444", "csv$Rrs445", "csv$Rrs446", "csv$Rrs447", "csv$Rrs448",
"csv$Rrs449", "csv$Rrs450", "csv$Rrs451", "csv$Rrs452", "csv$Rrs453", "csv$Rrs454", "csv$Rrs455",
"csv$Rrs456", "csv$Rrs457", "csv$Rrs458", "csv$Rrs459", "csv$Rrs460", "csv$Rrs461", "csv$Rrs462",
"csv$Rrs463", "csv$Rrs464", "csv$Rrs465", "csv$Rrs466", "csv$Rrs467", "csv$Rrs468", "csv$Rrs469",
"csv$Rrs470", "csv$Rrs471", "csv$Rrs472", "csv$Rrs473", "csv$Rrs474", "csv$Rrs475", "csv$Rrs476",
"csv$Rrs477", "csv$Rrs478", "csv$Rrs479", "csv$Rrs480", "csv$Rrs481", "csv$Rrs482", "csv$Rrs483",
"csv$Rrs484", "csv$Rrs485", "csv$Rrs486", "csv$Rrs487", "csv$Rrs488", "csv$Rrs489", "csv$Rrs490",
"csv$Rrs491", "csv$Rrs492", "csv$Rrs493", "csv$Rrs494", "csv$Rrs495", "csv$Rrs496", "csv$Rrs497",
"csv$Rrs498", "csv$Rrs499", "csv$Rrs500", "csv$Rrs501", "csv$Rrs502", "csv$Rrs503", "csv$Rrs504",
"csv$Rrs505", "csv$Rrs506", "csv$Rrs507", "csv$Rrs508", "csv$Rrs509", "csv$Rrs510", "csv$Rrs511",
"csv$Rrs512", "csv$Rrs513", "csv$Rrs514", "csv$Rrs515", "csv$Rrs516", "csv$Rrs517", "csv$Rrs518",
"csv$Rrs519", "csv$Rrs520", "csv$Rrs521", "csv$Rrs522", "csv$Rrs523", "csv$Rrs524", "csv$Rrs525",
"csv$Rrs526", "csv$Rrs527", "csv$Rrs528", "csv$Rrs529", "csv$Rrs530", "csv$Rrs531", "csv$Rrs532",
"csv$Rrs533", "csv$Rrs534", "csv$Rrs535", "csv$Rrs536", "csv$Rrs537", "csv$Rrs538", "csv$Rrs539",
"csv$Rrs540", "csv$Rrs541", "csv$Rrs542", "csv$Rrs543", "csv$Rrs544", "csv$Rrs545", "csv$Rrs546",
"csv$Rrs547", "csv$Rrs548", "csv$Rrs549", "csv$Rrs550", "csv$Rrs551", "csv$Rrs552", "csv$Rrs553",
"csv$Rrs554", "csv$Rrs555", "csv$Rrs556", "csv$Rrs557", "csv$Rrs558", "csv$Rrs559", "csv$Rrs560",
"csv$Rrs561", "csv$Rrs562", "csv$Rrs563", "csv$Rrs564", "csv$Rrs565", "csv$Rrs566", "csv$Rrs567",
"csv$Rrs568", "csv$Rrs569", "csv$Rrs570", "csv$Rrs571", "csv$Rrs572", "csv$Rrs573", "csv$Rrs574",
"csv$Rrs575", "csv$Rrs576", "csv$Rrs577", "csv$Rrs578", "csv$Rrs579", "csv$Rrs580", "csv$Rrs581",
"csv$Rrs582", "csv$Rrs583", "csv$Rrs584", "csv$Rrs585", "csv$Rrs586", "csv$Rrs587", "csv$Rrs588",
"csv$Rrs589", "csv$Rrs590", "csv$Rrs591", "csv$Rrs592", "csv$Rrs593", "csv$Rrs594", "csv$Rrs595",
"csv$Rrs596", "csv$Rrs597", "csv$Rrs598", "csv$Rrs599", "csv$Rrs600", "csv$Rrs601", "csv$Rrs602",
"csv$Rrs603", "csv$Rrs604", "csv$Rrs605", "csv$Rrs606", "csv$Rrs607", "csv$Rrs608", "csv$Rrs609",
"csv$Rrs610", "csv$Rrs611", "csv$Rrs612", "csv$Rrs613", "csv$Rrs614", "csv$Rrs615", "csv$Rrs616",
"csv$Rrs617", "csv$Rrs618", "csv$Rrs619", "csv$Rrs620", "csv$Rrs621", "csv$Rrs622", "csv$Rrs623",
"csv$Rrs624", "csv$Rrs625", "csv$Rrs626", "csv$Rrs627", "csv$Rrs628", "csv$Rrs629", "csv$Rrs630",
"csv$Rrs631", "csv$Rrs632", "csv$Rrs633", "csv$Rrs634", "csv$Rrs635", "csv$Rrs636", "csv$Rrs637",
"csv$Rrs638", "csv$Rrs639", "csv$Rrs640", "csv$Rrs641", "csv$Rrs642", "csv$Rrs643", "csv$Rrs644",
"csv$Rrs645", "csv$Rrs646", "csv$Rrs647", "csv$Rrs648", "csv$Rrs649", "csv$Rrs650", "csv$Rrs651",
"csv$Rrs652", "csv$Rrs653", "csv$Rrs654", "csv$Rrs655", "csv$Rrs656", "csv$Rrs657", "csv$Rrs658",
"csv$Rrs659", "csv$Rrs660", "csv$Rrs661", "csv$Rrs662", "csv$Rrs663", "csv$Rrs664", "csv$Rrs665",
"csv$Rrs666", "csv$Rrs667", "csv$Rrs668", "csv$Rrs669", "csv$Rrs670", "csv$Rrs671", "csv$Rrs672",
"csv$Rrs673", "csv$Rrs674", "csv$Rrs675", "csv$Rrs676", "csv$Rrs677", "csv$Rrs678", "csv$Rrs679",
"csv$Rrs680", "csv$Rrs681", "csv$Rrs682", "csv$Rrs683", "csv$Rrs684", "csv$Rrs685", "csv$Rrs686",
"csv$Rrs687", "csv$Rrs688", "csv$Rrs689", "csv$Rrs690", "csv$Rrs691", "csv$Rrs692", "csv$Rrs693",
"csv$Rrs694", "csv$Rrs695", "csv$Rrs696", "csv$Rrs697", "csv$Rrs698", "csv$Rrs699", "csv$Rrs700",
"csv$Rrs701", "csv$Rrs702", "csv$Rrs703", "csv$Rrs704", "csv$Rrs705", "csv$Rrs706", "csv$Rrs707",
"csv$Rrs708", "csv$Rrs709", "csv$Rrs710", "csv$Rrs711", "csv$Rrs712", "csv$Rrs713", "csv$Rrs714",
"csv$Rrs715", "csv$Rrs716", "csv$Rrs717", "csv$Rrs718", "csv$Rrs719", "csv$Rrs720", "csv$Rrs721",
"csv$Rrs722", "csv$Rrs723", "csv$Rrs724", "csv$Rrs725", "csv$Rrs726", "csv$Rrs727", "csv$Rrs728",
"csv$Rrs729", "csv$Rrs730", "csv$Rrs731", "csv$Rrs732", "csv$Rrs733", "csv$Rrs734", "csv$Rrs735",
"csv$Rrs736", "csv$Rrs737", "csv$Rrs738", "csv$Rrs739", "csv$Rrs740", "csv$Rrs741", "csv$Rrs742",
"csv$Rrs743", "csv$Rrs744", "csv$Rrs745", "csv$Rrs746", "csv$Rrs747", "csv$Rrs748", "csv$Rrs749",
"csv$Rrs750", "csv$Rrs751", "csv$Rrs752", "csv$Rrs753", "csv$Rrs754", "csv$Rrs755", "csv$Rrs756",
"csv$Rrs757", "csv$Rrs758", "csv$Rrs759", "csv$Rrs760", "csv$Rrs761", "csv$Rrs762", "csv$Rrs763",
"csv$Rrs764", "csv$Rrs765", "csv$Rrs766", "csv$Rrs767", "csv$Rrs768", "csv$Rrs769", "csv$Rrs770",
"csv$Rrs771", "csv$Rrs772", "csv$Rrs773", "csv$Rrs774", "csv$Rrs775", "csv$Rrs776", "csv$Rrs777",
"csv$Rrs778", "csv$Rrs779", "csv$Rrs780", "csv$Rrs781", "csv$Rrs782", "csv$Rrs783", "csv$Rrs784",
"csv$Rrs785", "csv$Rrs786", "csv$Rrs787", "csv$Rrs788", "csv$Rrs789", "csv$Rrs790", "csv$Rrs791",
"csv$Rrs792", "csv$Rrs793", "csv$Rrs794", "csv$Rrs795", "csv$Rrs796", "csv$Rrs797", "csv$Rrs798",
"csv$Rrs799", "csv$Rrs800", "csv$Rrs801", "csv$Rrs802", "csv$Rrs803", "csv$Rrs804", "csv$Rrs805",
"csv$Rrs806", "csv$Rrs807", "csv$Rrs808", "csv$Rrs809", "csv$Rrs810", "csv$Rrs811", "csv$Rrs812",
"csv$Rrs813", "csv$Rrs814", "csv$Rrs815", "csv$Rrs816", "csv$Rrs817", "csv$Rrs818", "csv$Rrs819",
"csv$Rrs820", "csv$Rrs821", "csv$Rrs822", "csv$Rrs823", "csv$Rrs824", "csv$Rrs825", "csv$Rrs826",
"csv$Rrs827", "csv$Rrs828", "csv$Rrs829", "csv$Rrs830", "csv$Rrs831", "csv$Rrs832", "csv$Rrs833",
"csv$Rrs834", "csv$Rrs835", "csv$Rrs836", "csv$Rrs837", "csv$Rrs838", "csv$Rrs839", "csv$Rrs840",
"csv$Rrs841", "csv$Rrs842", "csv$Rrs843", "csv$Rrs844", "csv$Rrs845", "csv$Rrs846", "csv$Rrs847",
"csv$Rrs848", "csv$Rrs849", "csv$Rrs850", "csv$Rrs851", "csv$Rrs852", "csv$Rrs853", "csv$Rrs854",
"csv$Rrs855", "csv$Rrs856", "csv$Rrs857", "csv$Rrs858", "csv$Rrs859", "csv$Rrs860", "csv$Rrs861",
"csv$Rrs862", "csv$Rrs863", "csv$Rrs864", "csv$Rrs865", "csv$Rrs866", "csv$Rrs867", "csv$Rrs868",
"csv$Rrs869", "csv$Rrs870", "csv$Rrs871", "csv$Rrs872", "csv$Rrs873", "csv$Rrs874", "csv$Rrs875",
"csv$Rrs876", "csv$Rrs877", "csv$Rrs878", "csv$Rrs879", "csv$Rrs880", "csv$Rrs881", "csv$Rrs882",
"csv$Rrs883", "csv$Rrs884", "csv$Rrs885", "csv$Rrs886", "csv$Rrs887", "csv$Rrs888", "csv$Rrs889",
"csv$Rrs890", "csv$Rrs891", "csv$Rrs892", "csv$Rrs893", "csv$Rrs894", "csv$Rrs895", "csv$Rrs896",
"csv$Rrs897", "csv$Rrs898", "csv$Rrs899", "csv$Rrs900", "csv$Rrs901", "csv$Rrs902", "csv$Rrs903",
"csv$Rrs904", "csv$Rrs905", "csv$Rrs906", "csv$Rrs907", "csv$Rrs908", "csv$Rrs909", "csv$Rrs910",
"csv$Rrs911", "csv$Rrs912", "csv$Rrs913", "csv$Rrs914", "csv$Rrs915", "csv$Rrs916", "csv$Rrs917",
"csv$Rrs918", "csv$Rrs919", "csv$Rrs920", "csv$Rrs921", "csv$Rrs922", "csv$Rrs923", "csv$Rrs924",
"csv$Rrs925", "csv$Rrs926", "csv$Rrs927", "csv$Rrs928", "csv$Rrs929", "csv$Rrs930", "csv$Rrs931",
"csv$Rrs932", "csv$Rrs933", "csv$Rrs934", "csv$Rrs935", "csv$Rrs936", "csv$Rrs937", "csv$Rrs938",
"csv$Rrs939", "csv$Rrs940", "csv$Rrs941", "csv$Rrs942", "csv$Rrs943", "csv$Rrs944", "csv$Rrs945",
"csv$Rrs946", "csv$Rrs947", "csv$Rrs948", "csv$Rrs949", "csv$Rrs950", "csv$Rrs951", "csv$Rrs952",
"csv$Rrs953", "csv$Rrs954", "csv$Rrs955", "csv$Rrs956", "csv$Rrs957", "csv$Rrs958", "csv$Rrs959",
"csv$Rrs960", "csv$Rrs961", "csv$Rrs962", "csv$Rrs963", "csv$Rrs964", "csv$Rrs965", "csv$Rrs966",
"csv$Rrs967", "csv$Rrs968", "csv$Rrs969", "csv$Rrs970", "csv$Rrs971", "csv$Rrs972", "csv$Rrs973",
"csv$Rrs974", "csv$Rrs975", "csv$Rrs976", "csv$Rrs977", "csv$Rrs978", "csv$Rrs979", "csv$Rrs980",
"csv$Rrs981", "csv$Rrs982", "csv$Rrs983", "csv$Rrs984", "csv$Rrs985", "csv$Rrs986", "csv$Rrs987",
"csv$Rrs988", "csv$Rrs989", "csv$Rrs990", "csv$Rrs991", "csv$Rrs992", "csv$Rrs993", "csv$Rrs994",
"csv$Rrs995", "csv$Rrs996", "csv$Rrs997", "csv$Rrs998", "csv$Rrs999", "csv$Rrs1000")

hypers_headers2 <- c(
"csv$Rrs550", "csv$Rrs670", "csv$Rrs701",
"csv$Rrs713", "csv$Rrs850"
)

hypers_headers3 <- c('csv$Rrs550', 'csv$Rrs670', 'csv$Rrs631', 'csv$Rrs704',
'csv$Rrs717', 'csv$Rrs701', 'csv$Rrs713', 'csv$Rrs850')

hypers_headers4 <- c(
'csv$Rrs975',
'csv$Rrs907',
'csv$Rrs926',
'csv$Rrs969',
'csv$Rrs629',
'csv$Rrs712',
'csv$Rrs814',
'csv$Rrs830',
'csv$Rrs832',
'csv$Rrs978',
'csv$Rrs920',
'csv$Rrs429',
'csv$Rrs766',
'csv$Rrs739',
'csv$Rrs655',
'csv$Rrs711')


equations <- c(
'(a^2)^5.478/(b^1.133*c)^5.049',
'((a/b)^3.1/(c/d)^3.1)^0.3*(e/f)^3.1/(g/h)^3.1',
'a/(b-c)',
'a/(b+c)',
'a/(b*c)',
'(a-b)/(c-d)',
'(a*b)/(c*d)',
'(a+b)/(c+d)',
'(a-b)/(c+d)',
'(a+b)/(c-d)',
'a-b','(((a*b)^3.1/(c*d)^3.1)^3.1)+(e/f)^3.1',
  '(((a*b^3.1)/(c*d^3.2)))+(e/f)^3.1',
'(((a*b^3.1)/(c*d^3.2)))-(e/f)^3.1',
'(((a/b)^3.1/(c/d)^3.1) * ((e/f)^3.1/(g/h)^3.1))'
)
equations2 <- c(
  '(a^b)/2'
)
#equations = c('(a^2)^5.478/(b^1.133*c)^5.049')

# for (equation in equations) {
#   function_permutation('/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/hypers_total_rrs.csv', hypers_headers2, equation, 0.8)
# }


#eq = '((csv$RedEdge_717*csv$Red_667^%f)/(csv$NIR_831*csv$Green_557^(%f+0.1)))+(csv$Blue_482/csv$RedEdge_717)'
eqs = c('csv$Red_667^%f/(csv$Green_557+csv$NIR_831)^(%f+0.1)',
'((csv$Rrs633/csv$Rrs816)/(csv$Rrs701/csv$Rrs633))^0.3*(csv$Rrs713/csv$Rrs701)/(csv$Rrs563/csv$Rrs633)',
'((csv$Rrs633/csv$Rrs701)/(csv$Rrs816/csv$Rrs633))^0.3*(csv$Rrs633/csv$Rrs563)/(csv$Rrs701/csv$Rrs713)',
'((csv$Rrs633/csv$Rrs701)/(csv$Rrs816/csv$Rrs633))^0.3*(csv$Rrs633/csv$Rrs701)/(csv$Rrs563/csv$Rrs713)',
'((csv$Rrs633/csv$Rrs701)/(csv$Rrs816/csv$Rrs633))^0.3*(csv$Rrs713/csv$Rrs563)/(csv$Rrs701/csv$Rrs633)',
'((csv$Rrs633/csv$Rrs701)/(csv$Rrs816/csv$Rrs633))^0.3*(csv$Rrs713/csv$Rrs701)/(csv$Rrs563/csv$Rrs633)',
'((csv$Rrs633/csv$Rrs816)/(csv$Rrs701/csv$Rrs633))^0.3*(csv$Rrs633/csv$Rrs563)/(csv$Rrs563/csv$Rrs633)',
'((csv$Rrs633/csv$Rrs816)/(csv$Rrs701/csv$Rrs633))^0.3*(csv$Rrs633/csv$Rrs563)/(csv$Rrs701/csv$Rrs713)',
'((csv$Rrs633/csv$Rrs816)/(csv$Rrs701/csv$Rrs633))^0.3*(csv$Rrs633/csv$Rrs701)/(csv$Rrs563/csv$Rrs713)',
'((csv$Rrs633/csv$Rrs816)/(csv$Rrs701/csv$Rrs633))^0.3*(csv$Rrs713/csv$Rrs563)/(csv$Rrs701/csv$Rrs633)',
'((csv$Rrs633/csv$Rrs816)/(csv$Rrs701/csv$Rrs633))^0.3*(csv$Rrs713/csv$Rrs701)/(csv$Rrs563/csv$Rrs633)',
'((csv$Rrs633/csv$Rrs701)/(csv$Rrs816/csv$Rrs633))^0.3*(csv$Rrs633/csv$Rrs563)/(csv$Rrs701/csv$Rrs713)',
'((csv$Rrs633/csv$Rrs701)/(csv$Rrs816/csv$Rrs633))^0.3*(csv$Rrs633/csv$Rrs701)/(csv$Rrs563/csv$Rrs713)',
'((csv$Rrs633/csv$Rrs701)/(csv$Rrs816/csv$Rrs633))^0.3*(csv$Rrs713/csv$Rrs563)/(csv$Rrs701/csv$Rrs633)',
'((csv$Rrs633/csv$Rrs701)/(csv$Rrs816/csv$Rrs633))^0.3*(csv$Rrs713/csv$Rrs701)/(csv$Rrs563/csv$Rrs633)',
'((csv$Rrs633/csv$Rrs816)/(csv$Rrs701/csv$Rrs633))^0.3*(csv$Rrs633/csv$Rrs563)/(csv$Rrs701/csv$Rrs713)',
'((csv$Rrs633/csv$Rrs816)/(csv$Rrs701/csv$Rrs633))^0.3*(csv$Rrs633/csv$Rrs701)/(csv$Rrs563/csv$Rrs713)',
'((csv$Rrs633/csv$Rrs816)/(csv$Rrs701/csv$Rrs633))^0.3*(csv$Rrs713/csv$Rrs563)/(csv$Rrs701/csv$Rrs633)'
)



# optimize_algorithm('/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/hypers_total_rrs.csv', hypers_header2, 0.6)
# optimize_algorithm('/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/hypers_total_rrs.csv', hypers_header3, 0.6)
#optimize_algorithm('/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/hypers_total_rrs.csv', hypers_headers3, equations, 0.7)
#optimize_algorithm('/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/WMS_Rrs_micasense_SPM_alg_31_51_out.csv', mic_headers,  equations, 0.6)

# optimize_algorithm('/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/Hyperspectral_SPM_2020.csv', hypers_headers, equations2, 0.4)

optimize_algorithm('/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/modis_algorithm.csv', modis_headers,  equations, 0.68)

# choose_best_cols('/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/hypers_total_rrs.csv', 2:602)

# choose_best_cols('/mnt/13aa104a-192c-43e5-95af-68aba6ac57a9/MSU/SPM/mic_data.csv', 2:6)