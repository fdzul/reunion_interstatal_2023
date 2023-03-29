risk_transmision_activa <- function(betas, 
                                    geocoded_data,
                                    locality,
                                    cve_edo){
    
    # 1.1. aplicar la prueba de knox
    z <- geocoded_data |>
        as.data.frame() |>
        sf::st_as_sf(coords = c("long", "lat"),
                     crs = 4326) |>
        dplyr::mutate(onset = lubridate::ymd(FEC_INI_SIGNOS_SINT))
    
    # Step 1.2 Extract the locality
    loc <- rgeomex::extract_locality(cve_edo = cve_edo,
                                     locality = locality)
    library(magrittr)
    
    # Step 1.3. knox test
    knox_res <- denhotspots::knox(x = z[loc,], 
                                  crs = "+proj=eqc", 
                                  ds = 400, 
                                  dt = 20, 
                                  sym = 1000, 
                                  sp_link = FALSE,
                                  planar_coord = FALSE)
    
    ############
    # Step 1.4. calculate the week factor ####
    z <- knox_res$x |>
        dplyr::mutate(week = lubridate::epiweek(onset)) |>
        dplyr::mutate(week_factor = ifelse(week <= 10, "1-10",
                                           ifelse(week > 10 & week <= 20, "11-20",
                                                  ifelse(week > 20 & week <= 25, "21-25",
                                                         ifelse(week > 25 & week <= 30, "26-30",
                                                                ifelse(week > 30 & week <= 35, "31-35",
                                                                       ifelse(week > 35 & week <= 40, "36-40",
                                                                              ifelse(week > 40 & week <= 45, "41-45",
                                                                                     ifelse(week > 45 & week <= 53, "46-53",
                                                                                            NA))))))))) |>
        sf::st_as_sf(coords = c("x", "y"),
                     crs = "+proj=eqc") |>
        sf::st_transform(crs = 4326)
    
    # Step 1.5. load the Space-Time link ####
    st_link <- knox_res$space_time_link |>
        sf::st_set_crs(value = 4326)
    
    # Step 1.6. extract the dengue cases of space links ####
    w <- z[st_link,] |>
        dplyr::mutate(week = lubridate::epiweek(onset)) |>
        dplyr::mutate(week_factor = ifelse(week <= 10, "1-10",
                                           ifelse(week > 10 & week <= 20, "11-20",
                                                  ifelse(week > 20 & week <= 25, "21-25",
                                                         ifelse(week > 25 & week <= 30, "26-30",
                                                                ifelse(week > 30 & week <= 35, "31-35",
                                                                       ifelse(week > 35 & week <= 40, "36-40",
                                                                              ifelse(week > 40 & week <= 45, "41-45",
                                                                                     ifelse(week > 45 & week <= 53, "46-53",
                                                                                            NA))))))))) %>%
        dplyr::mutate(x = sf::st_coordinates(geometry)[,1],
                      y = sf::st_coordinates(geometry)[,2])
    
    # Step 1.7. add the week  to space-time link  ####
    st_link_week <- sf::st_join(x = st_link,
                                y = w[, c("week_factor")])
    
    
    #######
    # Step 2.1 convert the df to sf object 
    betas <- betas |>
        dplyr::mutate(long = x,
                      lat = y) |>
        sf::st_as_sf(coords = c("long", "lat"),
                     crs = 4326)
    
    # Step 2.2 extract the eggs hotspots of locality
    betas <- betas[loc, ] 
    
    # Step 2.3. calculate the intensity 
    intensity_function <- function(x){
        y <- x |>
            dplyr::mutate(hotspots_binary = ifelse(hotspots == "Hotspots", 1, 0)) |>
            as.data.frame() |>
            dplyr::select(x, y, week, hotspots_binary) |>
            tidyr::pivot_wider(id_cols = c(x, y),
                               names_from = "week",
                               #names_prefix = "hotspots",
                               values_from = "hotspots_binary") |>
            as.data.frame() 
        
        y$intensity <- rowSums(y |> dplyr::select(-1, -2))
        y$per_intensity <- round((y$intensity/ncol(y |> dplyr::select(-1, -2, -intensity)))*100,digits = 1)
        y |> dplyr::select(x, y, intensity, per_intensity)
    }
    
    betas <- betas |>
        sf::st_drop_geometry() |>
        dplyr::group_by(year) |>
        tidyr::nest() |>
        dplyr::mutate(intensity = purrr::map(data, intensity_function)) |>
        dplyr::select(-data) |>
        tidyr::unnest(cols = c(intensity))
    
    #########
    ggplot2::ggplot() +
        ggplot2::geom_sf(data = loc,
                         fill = "gray92",
                         col = "white",
                         lwd = 0.01) +
        ggplot2::geom_tile(data = betas |> 
                               dplyr::filter(intensity >= 1),
                           #col ="white",
                           ggplot2::aes(x = x,
                                        y = y,
                                        fill = intensity)) +
        ggplot2::scale_fill_viridis_c() +
        ggplot2::theme_void() +
        ggplot2::geom_sf(data = z,
                         col = "gray40",
                         shape = 19,
                         size = .5) +
        ggnewscale::new_scale_fill() +
        ggplot2::geom_sf(data = w,
                         ggplot2::aes(fill = week_factor),
                         col = "white",
                         shape = 21,
                         stroke = .1,
                         size = 2.5) +
        fishualize::scale_fill_fish_d("Semana",
                                      option = "Scarus_hoefleri",
                                      direction = -1) +
        ggnewscale::new_scale_color() +
        ggplot2::geom_sf(data = st_link_week,
                         ggplot2::aes(col = week_factor),
                         lwd = .5) +
        fishualize::scale_colour_fish_d("Space-Time Links",
                                        option = "Scarus_hoefleri",
                                        direction = 1)
    
}