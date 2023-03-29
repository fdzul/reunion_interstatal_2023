intensity_eggs <- function(betas,
                           locality,
                           cve_ent,
                           palette,
                           name){
    # Step 1. extract the locality 
    
    locality <- rgeomex::extract_locality(cve_edo = cve_ent, 
                                          locality = locality)
    
    # step 2. extract the eggs hotspots of locality ####
    
    # Step 2.1 convert the df to sf object 
    x <- betas |>
        dplyr::mutate(long = x,
                      lat = y) |>
        sf::st_as_sf(coords = c("long", "lat"),
                     crs = 4326)
    
    # Step 2.2 extract the eggs hotspots of locality
    x <- x[locality, ] 
    
    # Step 3. calculate the intensity 
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
    
    x <- x |>
        sf::st_drop_geometry() |>
        dplyr::group_by(year) |>
        tidyr::nest() |>
        dplyr::mutate(intensity = purrr::map(data,intensity_function)) |>
        dplyr::select(-data) |>
        tidyr::unnest(cols = c(intensity))
    
    
    # step 4 plot the map ####
    ggplot2::ggplot() +
        ggplot2::geom_tile(data = x,
                           ggplot2::aes(x = x,
                                        y = y,
                                        fill = intensity)) +
        ggplot2::scale_color_gradientn("Intensidad",
                                       colors = c("gray100", palette(n = max(x$intensity),
                                                                     name = name)),
                                       breaks = seq(0, max(x$intensity), 2),
                                       aesthetics = c("fill")) +
        ggplot2::geom_sf(data = locality,  
                         alpha = 1, 
                         fill = NA,
                         col = "black", 
                         lwd = 0.5) +
        ggplot2::facet_wrap(facets = "year") +
        ggplot2::theme_void() +
        ggplot2::theme(legend.position = "bottom") +
        ggplot2::theme(legend.key.size = ggplot2::unit(.8, "cm"),
                       legend.key.width = ggplot2::unit(.5,"cm"),
                       legend.margin= ggplot2::margin(0,0,0,0),
                       legend.box.margin= ggplot2::margin(-20,0,0,0)) +
        ggplot2::theme(legend.text = ggplot2::element_text(colour = "black",
                                                           face  = "bold"),
                       legend.title = ggplot2::element_text(colour = "darkred",
                                                            face  = "bold")) +
        ggplot2::theme(strip.text = ggplot2::element_text(size = 11,
                                                          face = "bold"))
    
}