#' Cases of dengue per week.
#'
#' @param path is the path of the sinave files
#' @param mpo is the name of municipality
#' @param year is the year the dataset
#' @name name is the name
#' @return a ggplot2 object
#' @export
#'
#' @examples 1+1

cases_week <- function(path = 'DENGUE2_.txt', edo = NULL, mpo = NULL, year = year(Sys.Date())){
    # Hospitalizados ####
    if (edo == NULL & mpo == NULL){
    a <- readr::read_delim(path, 
                           delim = '|',
                           locale = readr::locale(encoding = 'latin1'),
                           comment = '',
                           quote = '') %>%
        dplyr::filter(ANO == year, ESTATUS_CASO  %in% c(1, 2))
        } else if (edo != NULL & mpo == NULL){
        a <- readr::read_delim(path, 
                           delim = '|',
                           locale = readr::locale(encoding = 'latin1'),
                           comment = '',
                           quote = '') %>%
        dplyr::filter(ANO == year, ESTATUS_CASO  %in% c(1, 2)) %>%
        dplyr::filter(stringr::str_detect(DES_EDO_REP, pattern = edo))
        } else if (edo == NULL & mpo != NULL){
        a <- readr::read_delim(path, 
                           delim = '|',
                           locale = readr::locale(encoding = 'latin1'),
                           comment = '',
                           quote = '') %>%
        dplyr::filter(ANO == year, ESTATUS_CASO  %in% c(1, 2)) %>%
        dplyr::filter(stringr::str_detect(DES_MPO_REP, pattern = mpo))
        } else if (edo != NULL & mpo != NULL){
        a <- readr::read_delim(path, 
                           delim = '|',
                           locale = readr::locale(encoding = 'latin1'),
                           comment = '',
                           quote = '') %>%
        dplyr::filter(ANO == year, ESTATUS_CASO  %in% c(1, 2)) %>%
        dplyr::filter(stringr::str_detect(DES_EDO_REP, pattern = edo),
                      stringr::str_detect(DES_MPO_REP, pattern = mpo))
    }
    a <- a %>%
        dplyr::count(SEM, MANEJO) %>%
        dplyr::mutate(MANEJO = ifelse(MANEJO == 1,
                                      "Hospitalizado",
                                      "Ambulatorio")) %>%
        ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = SEM,
                                        y = n,
                                        col = MANEJO),
                           linewidth = 2) +
        ggplot2::geom_point(ggplot2::aes(x = SEM,
                                         y = n,
                                         fill = MANEJO),
                            color = "white",
                            shape = 21,
                            size = 5) +
        #ggplot2::facet_wrap("DES_MPO_RES",scales = "free_y") +
        ggplot2::theme(legend.position = c(.1, .85),
                       legend.background = ggplot2::element_blank()) +
        ggplot2::scale_color_manual(" ",
                                    values = c("#D7DF71", "#0E5F76")) +
        ggplot2::scale_fill_manual(" ",
                                   values = c("#D7DF71", "#0E5F76")) +
        ggplot2::scale_x_continuous(breaks = seq(from = 0,
                                                 to = 52,
                                                 by = 2)) +
        ggplot2::ylab("Número de Casos") +
        ggplot2::xlab("Semanas Epidemiológicas") +
        ggplot2::theme(strip.text = ggplot2::element_text(size = 12),
                       legend.background = ggplot2::element_blank(),
                       legend.key = ggplot2::element_blank(),
                       legend.box = ggplot2::element_blank())

    # estimados ####
        if (edo == NULL & mpo == NULL){
    b <- readr::read_delim(path, 
                           delim = '|',
                           locale = readr::locale(encoding = 'latin1'),
                           comment = '',
                           quote = '') %>%
        dplyr::filter(ANO == year, ESTATUS_CASO  %in% c(1, 2))
        } else if (edo != NULL & mpo == NULL){
        b <- readr::read_delim(path, 
                           delim = '|',
                           locale = readr::locale(encoding = 'latin1'),
                           comment = '',
                           quote = '') %>%
        dplyr::filter(ANO == year, ESTATUS_CASO  %in% c(1, 2)) %>%
        dplyr::filter(stringr::str_detect(DES_EDO_REP, pattern = edo))
        } else if (edo == NULL & mpo != NULL){
        b <- readr::read_delim(path, 
                           delim = '|',
                           locale = readr::locale(encoding = 'latin1'),
                           comment = '',
                           quote = '') %>%
        dplyr::filter(ANO == year, ESTATUS_CASO  %in% c(1, 2)) %>%
        dplyr::filter(stringr::str_detect(DES_MPO_REP, pattern = mpo))
        } else if (edo != NULL & mpo != NULL){
        b <- readr::read_delim(path, 
                           delim = '|',
                           locale = readr::locale(encoding = 'latin1'),
                           comment = '',
                           quote = '') %>%
        dplyr::filter(ANO == year, ESTATUS_CASO  %in% c(1, 2)) %>%
        dplyr::filter(stringr::str_detect(DES_EDO_REP, pattern = edo),
                      stringr::str_detect(DES_MPO_REP, pattern = mpo))
    }
    b <- b %>%
        count(DES_MPO_REP, SEM, ESTATUS_CASO) %>%
        tidyr::pivot_wider(names_from = ESTATUS_CASO,
                           values_fill = 0,
                           values_from = n) %>%
        dplyr::full_join(
            dplyr::tibble(DES_MPO_REP = NA, SEM = NA, `1` = NA, `2` = NA, `3` = NA)
            ) %>%
        dplyr::filter(!is.na(DES_MPO_REP)) %>%
        replace(is.na(.), 0 ) %>%
        dplyr::rename(probable = `1`,
                      confirmado =`2`,
                      descartado = `3`) %>%
        dplyr::mutate(Porcentaje_Positividad = round((confirmado/(confirmado + descartado))*100,
                                                     digits = 0),
                      estimados = ((Porcentaje_Positividad/100) * probable) + confirmado) %>%
        dplyr::select(DES_MPO_REP, SEM,
                      probable, confirmado,
                      estimados) %>%
        tidyr::pivot_longer(cols = c(probable, confirmado,
                                     estimados),
                            names_to = "status",
                            values_to = "n") |>
        ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = SEM,
                                        y = n,
                                        col = status),
                           linewidth = 2) +
        ggplot2::geom_point(ggplot2::aes(x = SEM,
                                         y = n,
                                         fill = status),
                            color = "white",
                            shape = 21,
                            size = 5) +
        ggplot2::theme(legend.position = c(.1, .85),
                       legend.background = ggplot2::element_blank()) +
        ggplot2::scale_color_manual("",
                                    values = c("#E01E5A", "#2EB67D", "#ECB32D")) +
        ggplot2::scale_fill_manual("",
                                   values = c("#E01E5A", "#2EB67D","#ECB32D")) +
        ggplot2::scale_x_continuous(breaks = seq(from = 0,
                                                 to = 52,
                                                 by = 2)) +
        ggplot2::theme(strip.text = ggplot2::element_text(size = 20),
                       legend.background = ggplot2::element_blank(),
                       legend.key = ggplot2::element_blank(),
                       legend.box = ggplot2::element_blank()) +
        ggplot2::scale_x_continuous(breaks = seq(from = 0,
                                                 to = 52,
                                                 by = 2)) +
        ggplot2::ylab("Número de Casos") +
        ggplot2::xlab("Semanas Epidemiológicas")

    a / b

}
