# Wang-Spindles-in-Axons
	Read Me
	The goal of this project is to create MATLAB code capable of detecting the spindle within the hippocampus and generate plots based on the result. 
	The “FD_Generation.m” is to generate a data structure recording the feed direction for each axon based on the previous study. The structure “Index_fd” generated in the workspace should be saved and will later be used. 
	The “General_Script.m” and “Spindle.m” should be used next. “Spindle.m” is the function used in the “General_Script.m”. These two MATLAB codes filter for spindles and generate the data that will be utilized next. The workspace of this code should be saved.
	The “plot_semilog_histogram.m” generates semi-log histograms separated based on the feed direction and regions. The plots will be normalized across each feed direction. This program also does the ANOVA test between each region separated by the feed direction. The workspace of “General_Script.m” should be loaded along with the structure “Index_fd”. The number of axons and spindles within that histogram and the median and mode of the fit will be displayed in the command window.
	The “plot_crossfd_ANOVA.m” is used to perform the ANOVA test across feed directions. The workspace of “General_Script.m” should be loaded along with the structure “Index_fd” as well.
	The “plot_polarhis.m” and “plot_polarhis_stats.m” generate the polar histogram for the phase angle of the spike within a spindle. The “plot_polarhis.m” generates the plot and the “plot_polarhis_stats.m” generates the effect size for the histogram. The workspace of “General_Script.m” should be loaded along with the structure “Index_fd” as well. The MATLAB toolbox of “Circular Statistics Toolbox” is needed.
	The “pot_fft.m” generates a power VS frequency graph for the 10-16Hz using the fast Fourier transform separated by region and feed direction. The workspace of “General_Script.m” should be loaded along with the structure “Index_fd”.
	The “plot_tunnel_tunnel_Avrg_Corr.m” finds the average correlation between each tunnel within the same regions. The workspace of “General_Script.m” should be loaded along with the structure “Index_fd”.
	The “General_Script_CA3.m” filter for spindles within the CA3 subregion. The result should be saved separately and will be used later. 
	The “well_tunnel_corr.m” generates a graph that illustrates the correlation between each well in the CA3 subregion and the tunnels connecting to it in all possible permutations. The Feedforward channels will be marked in red while the Feedback channels will be in blue. 
	

