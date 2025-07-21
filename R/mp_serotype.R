#' This function creates maps of serotypes by municipality.
#'
#' @param path_sinave is the directory where the files are located.
#' @param cve_edo is the id of state.
#' @param palette is is the color palette.
#' @author Felipe Antonio Dzul Manzanilla \email{felipe.dzul.m@gmail.com}
#' @return a mapview objetc
#' @export
#'
#' @examples

usethis::use_pipe(export = TRUE)


require(dplyr)
require(tidyr)
require(readr)
require(viridis)
require(stringr)
require(mapview)
require(rgeomex)
require(magrittr)

#' @export

mp_serotype <- function(path_sinave = 'DENGUE2_.txt', cve_edo = paste(str_pad(1:32, pad = "0", width = 2, side = "left")), palette = viridis::viridis){

    # Step 1. load the dengue dataset ####
    x_serotipo <- readr::read_delim(path_sinave, delim = "|",
                                    locale = reader::locale(encoding = "latin1"),
                                    quote = "",
                                    comment = "") %>%
        dplyr::filter(ESTATUS_CASO == 2) %>%
        dplyr::mutate(CVE_EDO_REP = stringr::str_pad(CVE_EDO_REP,
                                                     width = 2,
                                                     side = "left",
                                                     pad = "0")) %>%
        dplyr::mutate(CVE_MPO_REP = stringr::str_pad(CVE_MPO_REP,
                                                     width = 3,
                                                     side = "left",
                                                     pad = "0")) %>%
        dplyr::filter(DENGUE_SER_TRIPLEX %in% 1:4) %>%
        dplyr::count(CVE_EDO_REP, 
                     CVE_MPO_REP,
                     DENGUE_SER_TRIPLEX) %>%
        tidyr::pivot_wider(names_from = DENGUE_SER_TRIPLEX,
                           values_from = n,
                           names_prefix = "D",
                           values_fill = 0) %>%
        dplyr::full_join(
            dplyr::tibble(CVE_EDO_REP = NA, 
                         CVE_MPO_REP = NA,
                         D1 = NA,
                         D2 = NA,
                         D3 = NA,
                         D4 = NA)
        ) %>%
        dplyr::filter(!is.na(CVE_EDO_REP)) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(D1_binary = ifelse(D1 >0, 1, 0),
                      D2_binary = ifelse(D2 >0, 1, 0),
                      D3_binary = ifelse(D3 >0, 1, 0),
                      D4_binary = ifelse(D4 >0, 1, 0)) %>%
        dplyr::mutate(n_serotipo = D1_binary + D2_binary + D3_binary + D4_binary) %>%
        dplyr::mutate(D1_text = ifelse(D1 >0, "D1", ""),
                      D2_text = ifelse(D2 >0, "D2", ""),
                      D3_text = ifelse(D3 >0, "D3", ""),
                      D4_text = ifelse(D4 >0, "D4", "")) %>%
        dplyr::mutate(serotype_combination= paste(D1_text,
                                                  D2_text,
                                                  D3_text,
                                                  D4_text,
                                                  sep = "")) %>%
        dplyr::select(CVE_EDO_REP, CVE_MPO_REP,
                      D1, D2, D3, D4,
                      n_serotipo, serotype_combination) %>%
        dplyr::filter(CVE_EDO_REP %in% c(cve_edo))

    # Step 1. load the dengue dataset ####
    x_casos <- readr::read_delim(path_sinave, delim = "|",
                                    locale = reader::locale(encoding = "latin1"),
                                    quote = "",
                                    comment = "") %>%
        dplyr::filter(ESTATUS_CASO == 2) %>%
        dplyr::mutate(CVE_EDO_REP = stringr::str_pad(CVE_EDO_REP,
                                                     width = 2,
                                                     side = "left",
                                                     pad = "0")) %>%
        dplyr::mutate(CVE_MPO_REP = stringr::str_pad(CVE_MPO_REP,
                                                     width = 3,
                                                     side = "left",
                                                     pad = "0")) %>%
        dplyr::count(CVE_EDO_REP, 
                     CVE_MPO_REP) %>%
        dplyr::filter(CVE_EDO_REP %in% c(cve_edo))


    # Step 2. load boundary of municipality ####
    y <- rgeomex::AGEM_inegi19_mx %>%
        dplyr::filter(CVE_ENT  %in% c(cve_edo))


    # Step 3. left joint between munipality and serotype ####
    xy_serotipo <- dplyr::left_join(x = y,
                                    y = x_serotipo,
                                    by = c("CVE_ENT" = "CVE_EDO_REP",
                                           "CVE_MUN" = "CVE_MPO_REP")) %>%
        dplyr::filter(!is.na(n_serotipo))

    xy_casos <- dplyr::left_join(x = y,
                                 y = x_casos,
                                 by = c("CVE_ENT" = "CVE_EDO_REP",
                                        "CVE_MUN" = "CVE_MPO_REP")) %>%
        dplyr::filter(!is.na(n))



    map <- mapview::mapview(y,
                            col.regions = "gray90",
                            layer.name = "Municipios") +
        mapview::mapview(x = xy_casos,
                         col.regions = "#2EB67D",
                         layer.name = "Municipios Positivos") +
        mapview::mapview(x = xy_casos,
                         zcol = "n",
                         col.regions = palette,
                         layer.name = "Número de Casos")


    # Step 4. maps ####
    if(length(unique(xy_serotipo$n_serotipo)) == 4){
        if(sum(unique(xy_serotipo$n_serotipo)) == 10){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 1),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "1 Serotipos") +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 2),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "2 Serotipos") +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 3),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "3 Serotipos") +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 4),
                                 zcol = "serotype_combination",
                                 col.regions = "#E01E5A",
                                 layer.name = "4 Serotipos") + map
        } else if(sum(unique(xy_serotipo$n_serotipo)) == 9){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 2),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "2 Serotipo")  +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 3),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "3 Serotipos") +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 4),
                                 zcol = "serotype_combination",
                                 col.regions = "#E01E5A",
                                 layer.name = "4 Serotipos") + map
        } else if(sum(unique(xy_serotipo$n_serotipo)) == 8){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 1),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "1 Serotipo")  +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 3),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "3 Serotipos") +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 4),
                                 zcol = "serotype_combination",
                                 col.regions = "#E01E5A",
                                 layer.name = "4 Serotipos") + map

        } else if(sum(unique(xy_serotipo$n_serotipo)) == 7){
            if(unique(x_serotipo$n_serotype) %in% c(4, 2, 1)){
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 1),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "1 Serotipo")  +
                    mapview::mapview(x = xy_serotipo |>
                                         dplyr::filter(n_serotipo == 2),
                                     zcol = "serotype_combination",
                                     col.regions = palette,
                                     layer.name = "2 Serotipos") +
                    mapview::mapview(x = xy_serotipo |>
                                         dplyr::filter(n_serotipo == 4),
                                     zcol = "serotype_combination",
                                     col.regions = "#E01E5A",
                                     layer.name = "4 Serotipos") + map
            } else{

                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 3),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "3 Serotipos") +
                    mapview::mapview(x = xy_serotipo |>
                                         dplyr::filter(n_serotipo == 4),
                                     zcol = "serotype_combination",
                                     col.regions = "#E01E5A",
                                     layer.name = "4 Serotipos") + map

            }

        } else if(sum(unique(xy_serotipo$n_serotipo)) == 6){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 2),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "2 Serotipos") +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 4),
                                 zcol = "serotype_combination",
                                 col.regions = "#E01E5A",
                                 layer.name = "4 Serotipos") + map
        } else if(sum(unique(xy_serotipo$n_serotipo)) == 5){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 1),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "1 Serotipo")  +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 4),
                                 zcol = "serotype_combination",
                                 col.regions = "#E01E5A",
                                 layer.name = "4 Serotipos") + map
        } else if(sum(unique(xy_serotipo$n_serotipo)) == 4) {
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 4),
                             zcol = "serotype_combination",
                             col.regions = "#E01E5A",
                             layer.name = "4 Serotipos") + map
        } else{

        }

    } else if(length(unique(xy_serotipo$n_serotipo)) == 3){

        if(sum(unique(xy_serotipo$n_serotipo)) == 6){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 1),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "1 Serotipo")  +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 2),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "2 Serotipos") +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 3),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "3 Serotipos") + map
        } else if(sum(unique(xy_serotipo$n_serotipo)) == 7){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 1),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "1 Serotipo")  +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 2),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "2 Serotipos") +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 4),
                                 zcol = "serotype_combination",
                                 col.regions = "#E01E5A",
                                 layer.name = "4 Serotipos") + map
        } else if(sum(unique(xy_serotipo$n_serotipo)) == 8){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 1),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "1 Serotipo")  +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 3),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "3 Serotipos") +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 4),
                                 zcol = "serotype_combination",
                                 col.regions = "#E01E5A",
                                 layer.name = "4 Serotipos") + map
        } else if(sum(unique(xy_serotipo$n_serotipo)) == 9){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 2),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "2 Serotipos")  +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 3),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "3 Serotipos") +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 4),
                                 zcol = "serotype_combination",
                                 col.regions = "#E01E5A",
                                 layer.name = "4 Serotipos") + map
        } else{}

    } else if(length(unique(xy_serotipo$n_serotipo)) == 2){
        if(sum(unique(xy_serotipo$n_serotipo)) == 3){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 1),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "1 Serotipo")  +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 2),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "2 Serotipos") + map
        } else if(sum(unique(xy_serotipo$n_serotipo)) == 4){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 1),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "1 Serotipo")  +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 3),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "3 Serotipos") + map
        } else if(sum(unique(xy_serotipo$n_serotipo)) == 5){
            if(max(unique(xy_serotipo$n_serotipo)) == 4){
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 1),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "1 Serotipo")  +
                    mapview::mapview(x = xy_serotipo |>
                                         dplyr::filter(n_serotipo == 4),
                                     zcol = "serotype_combination",
                                     col.regions = "#E01E5A",
                                     layer.name = "4 Serotipos") + map
            } else{
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 2),
                                 zcol = "serotype_combination",
                                 col.regions = palette,
                                 layer.name = "2 Serotipos")  +
                    mapview::mapview(x = xy_serotipo |>
                                         dplyr::filter(n_serotipo == 3),
                                     zcol = "serotype_combination",
                                     col.regions = palette,
                                     layer.name = "3 Serotipos") + map

            }
        } else if(sum(unique(xy_serotipo$n_serotipo)) == 6){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 2),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "2 Serotipos")  +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 4),
                                 zcol = "serotype_combination",
                                 col.regions = "#E01E5A",
                                 layer.name = "4 Serotipos") + map
        } else if(sum(unique(xy_serotipo$n_serotipo)) == 7){
            mapview::mapview(x = xy_serotipo |>
                                 dplyr::filter(n_serotipo == 3),
                             zcol = "serotype_combination",
                             col.regions = palette,
                             layer.name = "3 Serotipos") +
                mapview::mapview(x = xy_serotipo |>
                                     dplyr::filter(n_serotipo == 4),
                                 zcol = "serotype_combination",
                                 col.regions = "#E01E5A",
                                 layer.name = "4 Serotipos") + map

        } else {

        }


    } else if(length(unique(xy_serotipo$n_serotipo)) == 1){
        mapview::mapview(x = xy_serotipo,
                         zcol = "serotype_combination",
                         col.regions = palette,
                         layer.name = "1 Serotipo") + map
    } else {}

}
