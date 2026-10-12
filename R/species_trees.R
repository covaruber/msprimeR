# Species tree parsing (port of msprime/species_trees.py), including a small
# Newick parser supporting quoted labels, branch lengths and [comments].

.newick_node <- function() {
  e <- new.env(parent = emptyenv())
  e$name <- NULL
  e$length <- 0
  e$comment <- NULL
  e$descendants <- list()
  e
}

.parse_newick_string <- function(s) {
  chars <- strsplit(s, "")[[1]]
  n <- length(chars)
  pos <- 1
  peek <- function() if (pos <= n) chars[pos] else ""
  skip_ws <- function() while (pos <= n && grepl("[[:space:]]", chars[pos])) pos <<- pos + 1
  read_comment <- function() {
    # Assumes chars[pos] == "["
    depth <- 0
    start <- pos
    repeat {
      if (pos > n) stop("Unterminated comment in newick string", call. = FALSE)
      ch <- chars[pos]
      if (ch == "[") depth <- depth + 1
      if (ch == "]") depth <- depth - 1
      pos <<- pos + 1
      if (depth == 0) break
    }
    paste(chars[(start + 1):(pos - 2)], collapse = "")
  }
  add_comment <- function(node, txt) {
    node$comment <- if (is.null(node$comment)) txt else paste0(node$comment, txt)
  }
  read_label <- function() {
    skip_ws()
    if (peek() == "'") {
      pos <<- pos + 1
      out <- character(0)
      repeat {
        if (pos > n) stop("Unterminated quoted label in newick string", call. = FALSE)
        if (chars[pos] == "'") {
          if (pos + 1 <= n && chars[pos + 1] == "'") {
            out <- c(out, "'")
            pos <<- pos + 2
            next
          }
          pos <<- pos + 1
          break
        }
        out <- c(out, chars[pos])
        pos <<- pos + 1
      }
      return(paste(out, collapse = ""))
    }
    start <- pos
    while (pos <= n && !chars[pos] %in% c("(", ")", ",", ":", ";", "[")) pos <<- pos + 1
    if (pos == start) return(NULL)
    lab <- paste(chars[start:(pos - 1)], collapse = "")
    lab <- trimws(lab)
    if (nzchar(lab)) lab else NULL
  }
  parse_node <- function() {
    node <- .newick_node()
    skip_ws()
    if (peek() == "(") {
      pos <<- pos + 1
      repeat {
        child <- parse_node()
        node$descendants[[length(node$descendants) + 1]] <- child
        skip_ws()
        ch <- peek()
        if (ch == ",") {
          pos <<- pos + 1
        } else if (ch == ")") {
          pos <<- pos + 1
          break
        } else {
          stop("Malformed newick string", call. = FALSE)
        }
      }
    }
    node$name <- read_label()
    skip_ws()
    while (peek() == "[") {
      add_comment(node, read_comment())
      skip_ws()
    }
    if (peek() == ":") {
      pos <<- pos + 1
      skip_ws()
      while (peek() == "[") {
        add_comment(node, read_comment())
        skip_ws()
      }
      start <- pos
      while (pos <= n && !chars[pos] %in% c("(", ")", ",", ";", "[") &&
             !grepl("[[:space:]]", chars[pos])) pos <<- pos + 1
      len <- suppressWarnings(as.double(paste(chars[start:(pos - 1)], collapse = "")))
      if (pos == start || is.na(len)) stop("Invalid branch length in newick string",
                                           call. = FALSE)
      node$length <- len
      skip_ws()
      while (peek() == "[") {
        add_comment(node, read_comment())
        skip_ws()
      }
    }
    node
  }
  s <- trimws(s)
  if (!nzchar(s)) return(NULL)
  root <- parse_node()
  skip_ws()
  if (peek() == ";") pos <- pos + 1
  skip_ws()
  if (pos <= n) stop(sprintf("Not a valid newick tree: '%s'", s), call. = FALSE)
  root
}

.newick_walk <- function(root, order = c("preorder", "postorder")) {
  order <- match.arg(order)
  out <- list()
  visit <- function(node) {
    if (order == "preorder") out[[length(out) + 1]] <<- node
    for (child in node$descendants) visit(child)
    if (order == "postorder") out[[length(out) + 1]] <<- node
  }
  visit(root)
  out
}

.parse_newick <- function(tree, branch_length_multiplier) {
  root <- tryCatch(.parse_newick_string(tree), error = function(e) {
    stop(sprintf("Not a valid newick tree: '%s'", tree), call. = FALSE)
  })
  if (is.null(root)) stop(sprintf("Not a valid newick tree: '%s'", tree), call. = FALSE)
  stack <- list(list(root, 0))
  num_nodes <- 0
  max_depth <- 0
  while (length(stack) > 0) {
    item <- stack[[length(stack)]]
    stack[[length(stack)]] <- NULL
    node <- item[[1]]
    depth <- item[[2]]
    if (depth > max_depth) max_depth <- depth
    num_nodes <- num_nodes + 1
    node$depth <- depth
    for (child in node$descendants) stack[[length(stack) + 1]] <- list(child, depth + child$length)
  }
  if (num_nodes < 3) stop("Newick tree must have at least three nodes", call. = FALSE)
  for (node in .newick_walk(root)) {
    node$time <- (max_depth - node$depth) * branch_length_multiplier
    if (length(node$descendants) == 0 && abs(node$time) > 1e-8) {
      .stopf("All leaf populations must be at time 0: time=%s", node$time)
    }
  }
  root
}

.check_generation_time <- function(generation_time) {
  g <- suppressWarnings(as.double(generation_time))
  if (length(g) != 1 || is.na(g)) stop("Generation time must be numeric.", call. = FALSE)
  if (g <= 0) stop("Generation time must be > 0.", call. = FALSE)
  g
}

.generations_per_branch_length_unit <- function(time_units, generation_time) {
  if (time_units == "gen") 1 else if (time_units == "myr") 10^6 / generation_time else
    1 / generation_time
}

.number_or_mapping <- function(value, message) {
  if (is.numeric(value) && length(value) == 1 && is.null(names(value))) {
    x <- as.double(value)
    return(function(name) x)
  }
  if ((is.list(value) || is.numeric(value)) && !is.null(names(value))) {
    return(function(name) {
      key <- name %||% ""
      if (!key %in% names(value)) {
        if (".default" %in% names(value)) return(as.double(value[[".default"]]))
        .stopf("No value specified for population '%s'", key)
      }
      as.double(value[[key]])
    })
  }
  stop(message, call. = FALSE)
}

#' Parse a species tree into a demographic model
#'
#' Parses a Newick species tree (rooted, ultrametric, with branch lengths) and
#' returns the corresponding [Demography]: one population per node, with
#' population splits at internal nodes. Used by `Demography$from_species_tree()`.
#'
#' @param tree The Newick string.
#' @param initial_size A single size, or a named vector/list mapping population
#'   names to sizes (an element named `.default` is used for unnamed nodes).
#' @param time_units Branch length units: "gen", "yr" or "myr".
#' @param generation_time Years per generation (required unless "gen").
#' @param growth_rate A single growth rate or a named mapping (default 0).
#' @return A [Demography].
#' @export
parse_species_tree <- function(tree, initial_size, time_units = "gen",
                               generation_time = NULL, growth_rate = NULL) {
  if (!time_units %in% c("myr", "yr", "gen")) {
    .stopf(paste0("The specified units for branch lengths (\"%s\") are not accepted. ",
                  "Accepted units are \"myr\" (millions of years), \"yr\" (years), and ",
                  "\"gen\" (generations)."), time_units)
  }
  size_of <- .number_or_mapping(initial_size, paste0(
    "initial_size argument must be a single number or a mapping from species names to ",
    "their population sizes."))
  growth_of <- .number_or_mapping(growth_rate %||% 0, paste0(
    "growth_rate argument must be a single number or a mapping from species names to ",
    "their exponential growth rates."))
  if (!is.null(generation_time)) generation_time <- .check_generation_time(generation_time)
  if (time_units == "gen") {
    if (!is.null(generation_time)) {
      stop("With branch lengths in units of generations (\"gen\"), a generation time ",
           "should not be specified additionally.", call. = FALSE)
    }
  } else if (is.null(generation_time)) {
    .stopf(paste0("With branch lengths in units of \"%s\", a generation time must be ",
                  "specified additionally."), time_units)
  }
  mult <- .generations_per_branch_length_unit(time_units, generation_time)
  root <- .parse_newick(tree, mult)
  demography <- Demography$new()
  add_population <- function(node) {
    name <- node$name
    if (!is.null(name)) {
      name <- trimws(name)
      if (!nzchar(name)) name <- NULL
    }
    pop <- demography$add_population(initial_size = size_of(name),
                                     growth_rate = growth_of(name), name = name)
    node$pop_name <- pop$name
    pop$name
  }
  for (node in .newick_walk(root)) {
    if (length(node$descendants) == 0) add_population(node)
  }
  for (node in .newick_walk(root, "postorder")) {
    if (length(node$descendants) > 0) {
      pop_name <- add_population(node)
      demography$add_population_split(
        time = node$time, derived = lapply(node$descendants, function(c) c$pop_name),
        ancestral = pop_name)
    }
  }
  demography$sort_events()
  demography$validate()
  demography
}

#' Parse a StarBEAST species tree
#'
#' Parses a species tree written by TreeAnnotator from a StarBEAST posterior
#' (Nexus format with `dmv` population size annotations) into a [Demography].
#' Used by `Demography$from_starbeast()`.
#'
#' @param tree The Nexus string (typically the whole file content).
#' @param generation_time Years per generation.
#' @param time_units "myr" (default) or "yr".
#' @return A [Demography].
#' @export
parse_starbeast <- function(tree, generation_time, time_units = "myr") {
  if (!time_units %in% c("myr", "yr")) {
    .stopf(paste0("The specified units for branch lengths (\"%s\") are not accepted. ",
                  "Accepted units are \"myr\" (millions of years) or \"yr\" (years)."),
           time_units)
  }
  generation_time <- .check_generation_time(generation_time)
  mult <- .generations_per_branch_length_unit(time_units, generation_time)
  parsed <- .parse_nexus(tree)
  name_map <- .parse_translate_command(parsed$translate)
  root <- .parse_newick(parsed$tree, mult)
  demography <- Demography$new()
  for (node in .newick_walk(root)) {
    if (is.null(node$comment)) stop("Annotation missing for one or more nodes.", call. = FALSE)
    m <- regmatches(node$comment, regexec("&dmv=\\{([0-9.]+?)\\}", node$comment))[[1]]
    if (length(m) < 2) stop("No dmv annotation for node", call. = FALSE)
    node$pop_size <- as.double(m[2]) * mult
  }
  add_population <- function(node) {
    name <- if (!is.null(node$name)) name_map[[node$name]] else NULL
    pop <- demography$add_population(initial_size = node$pop_size, name = name)
    node$pop_name <- pop$name
    pop$name
  }
  for (node in .newick_walk(root)) {
    if (length(node$descendants) == 0) add_population(node)
  }
  for (node in .newick_walk(root, "postorder")) {
    if (length(node$descendants) > 0) {
      pop_name <- add_population(node)
      demography$add_population_split(
        time = node$time, ancestral = pop_name,
        derived = lapply(node$descendants, function(c) c$pop_name))
    }
  }
  demography$sort_events()
  demography$validate()
  demography
}

.parse_translate_command <- function(translate_command) {
  stopifnot(substr(translate_command, 1, 10) == "translate ")
  items <- strsplit(substring(translate_command, 11), ",")[[1]]
  mapping <- list()
  for (item in items) {
    parts <- strsplit(trimws(item), "[[:space:]]+")[[1]]
    parts <- parts[nzchar(parts)]
    if (length(parts) <= 1) stop("Missing translation in translation block.", call. = FALSE)
    if (length(parts) != 2) {
      stop("Species IDs in the translation block appear to include whitespace. This is ",
           "not supported.", call. = FALSE)
    }
    if (parts[1] %in% names(mapping)) {
      .stopf("Newick ID %s defined multiple times in translation", parts[1])
    }
    mapping[[parts[1]]] <- parts[2]
  }
  if (anyDuplicated(unlist(mapping))) {
    stop("Duplicate species names in translation", call. = FALSE)
  }
  mapping
}

.parse_nexus <- function(nexus) {
  s <- gsub("\n", " ", nexus, fixed = TRUE)
  s <- paste(strsplit(trimws(s), "[[:space:]]+")[[1]], collapse = " ")
  s <- tolower(s)
  if (substr(s, 1, 6) != "#nexus") {
    stop("The species tree does not appear to be in Nexus format.", call. = FALSE)
  }
  s <- gsub(" ;", ";", s, fixed = TRUE)
  start <- regexpr("begin trees;", s, fixed = TRUE)
  if (start < 0) stop("The Nexus string does not include a complete trees block.", call. = FALSE)
  rest <- substring(s, start)
  end <- regexpr("end;", rest, fixed = TRUE)
  if (end < 0) stop("The Nexus string does not include a complete trees block.", call. = FALSE)
  block <- substr(rest, 1, end + 3)
  commands <- trimws(strsplit(block, ";", fixed = TRUE)[[1]])
  if (commands[1] != "begin trees" || commands[length(commands)] != "end" ||
      length(commands) < 3) {
    stop("Tree block malformed", call. = FALSE)
  }
  commands <- commands[-c(1, length(commands))]
  first_word <- vapply(strsplit(commands, " "), function(x) x[1], character(1))
  translate <- commands[first_word == "translate"]
  trees <- commands[first_word == "tree"]
  if (length(translate) != 1) {
    stop("The Nexus string does not contain exactly one translate command.", call. = FALSE)
  }
  if (length(trees) != 1) {
    stop("The Nexus string does not contain exactly one tree command.", call. = FALSE)
  }
  paren <- regexpr("(", trees, fixed = TRUE)
  if (paren < 0) stop("No parentheses in tree string", call. = FALSE)
  list(translate = translate, tree = substring(trees, paren))
}
