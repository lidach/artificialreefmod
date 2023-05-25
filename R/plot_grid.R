#' @title plot_grid
#'
#' @description function to plot spatial model (includes grid of sites, artificial and natural reef locations, landing sites, and sea depth in meters)
#'
#' @param coord_file data frame with depth, longitude, and latitude coordinates for each cell/site
#' @param AR_sites data frame with longitude and latitude coordinates of artificial reefs
#' @param landing_site_file data frame with longitude and latitude coordinates for landing sites
#' @param sea_col color for the coord_file (color of the ocean)
#' @param xlim_map x axis ranges for overall map (longitude); c(long1, long2)
#' @param ylim_map y axis ranges for overall map (latitude); c(lat1, lat2)
#' @param maxdepth maximum depth in map for plotting
#' @param save_plot save the plot to local directory? TRUE/FALSE
#' @param plot_legend
#' @param scenario_name
#' @param savedir directory name to save plot to, only works if save_plot = TRUE

plot_grid <-
function(coord_file,
	AR_sites,
	landing_site_file,
	sea_col,
	xlim_map,
	ylim_map,
	maxdepth,
	save_plot = FALSE,
	plot_legend = TRUE,
	scenario_name,
	savedir)
{
	coord_file_sp <- coord_file[,c("long","lat","depth")]
	colnames(coord_file_sp) <- c("x","y","depth")
	coordinates(coord_file_sp) <- ~x+y
	proj4string(coord_file_sp) <- CRS("+init=epsg:4326")
	coord_file_sp <- spTransform(coord_file_sp, CRS("+init=epsg:4326"))
	gridded(coord_file_sp) = TRUE
	r <- raster(coord_file_sp)
	projection(r) = CRS("+init=epsg:4326")
	newmap <- getMap(resolution = "high")

	NR_long <- coord_file[which(coord_file$HB>0),"long"]
	NR_lat <- coord_file[which(coord_file$HB>0),"lat"]

	AR_long <- coord_file[AR_sites,"long"]
	AR_lat <- coord_file[AR_sites,"lat"]

	## plot
	if(save_plot == TRUE)	jpeg(filename = paste0(savedir,"/",scenario_name,"_map.jpeg"), pointsize = 16, res = 300, units = "in", width = 8, height = 8)
		if(plot_legend == FALSE) par(oma = c(0,0,0,0))
		sp::plot(r, col = sea_col, legend = FALSE)
		sp::plot(newmap, add = TRUE, xlim = xlim_map, ylim = ylim_map, asp = 1, col = "gray95")
		if(plot_legend == TRUE)	sp::plot(r, add = TRUE, legend.only = TRUE, breaks = ceiling(seq(0,maxdepth,length=length(sea_col))), col = sea_col, legend.width = 1, legend.shrink = 0.5,
											axis.args = list(at = seq(0,maxdepth,50),
											labels = seq(0,maxdepth,50),
											cex.axis = 0.7))
		points(NR_long, NR_lat, cex = 1.5, pch = 0)
		points(AR_long, AR_lat, cex = 0.9, pch = 17)
		points(landing_site_file$long, landing_site_file$lat, cex = 1.3, pch = 13, col = "red")
		if(plot_legend) legend("bottomright", legend = c("natural reef","artificial reef", "landing site"), pch = c(0,17,13), col = c("black","black","red"), bg = 'white')
	if(save_plot == TRUE)  dev.off()
} # end of plot_grid function


#' @title res_site_plot
#'
#' @description function to plot results of models spatially
#'
#' @param coord_file data frame with depth, longitude, and latitude coordinates for each cell/site
#' @param AR_sites data frame with longitude and latitude coordinates of artificial reefs
#' @param landing_site_file data frame with longitude and latitude coordinates for landing sites
#' @param res_col color for the coord_file (color of the ocean)

res_site_plot <-
function(coord_file,
	AR_sites,
	res,
	res_col,
	cuts,
	save_plot = FALSE,
	scenario_name,
	savedir,
	plot_legend)
{
	coord_file_res <- cbind(coord_file[,c("long","lat")], res)
	colnames(coord_file_res) <- c("x","y","output")
	coordinates(coord_file_res) <- ~x+y
	proj4string(coord_file_res) <- CRS("+init=epsg:4326")
	coord_file_res <- spTransform(coord_file_res, CRS("+init=epsg:4326"))
	gridded(coord_file_res) = TRUE
	r <- raster(coord_file_res)
	projection(r) = CRS("+init=epsg:4326")

	NR_long <- coord_file[which(coord_file$HB>0),"long"]
	NR_lat <- coord_file[which(coord_file$HB>0),"lat"]

	AR_long <- coord_file[AR_sites,"long"]
	AR_lat <- coord_file[AR_sites,"lat"]

	## plot
	if(save_plot == TRUE) jpeg(filename = paste0(savedir,"/",scenario_name,"_siteoutput.jpeg"), pointsize = 16, res = 300, units = "in", width = 8, height = 8)
		sp::plot(r, col = res_col, breaks = cuts, asp = 1.5, box = FALSE, legend = FALSE, axes = FALSE)
		points(NR_long, NR_lat, cex = 1.5, pch = 0)
		points(AR_long, AR_lat, cex = 1, pch = 17)
		if(plot_legend) legend("bottomright", legend = c("natural reef","artificial reef"), pch = c(0,17), col = c("black","black"), bg = 'white')
	if(save_plot == TRUE) dev.off()
} # end of res_site_plot function