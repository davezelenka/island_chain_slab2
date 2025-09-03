Research on boundary-layers at slab interfaces at island chain subduction zones.

Island Arc Slab2.0 Global distribution of Slab2 models colored by depth to the slab surface: https://earthquake.usgs.gov/slab2/
Inspected with https://interactive-earth.com/earthquakes
Earthquake catalogue: https://earthquake.usgs.gov/earthquakes/search/
Slab2.0 Models found here: https://www.sciencebase.gov/catalog/item/5aa1b00ee4b0b1c392e86467/

USGS data conformed to the following structure:
['date_time', 'lon360', 'latitude', 'depth_km', 'magnitude']

================================

slab2_cross_section_extractor.py was used to determine distances from interface and slab center.

python slab2_cross_section_extractor.py --grd sco/sco_slab2_dep_02.23.18.grd --thk sco/sco_slab2_thk_02.23.18.grd --eq sco/central_islands.csv --csv sco/distances_surface.csv
python slab2_cross_section_extractor.py --grd izu/izu_slab2_dep_02.24.18.grd --thk izu/izu_slab2_thk_02.24.18.grd --eq izu/mariannas.csv --csv izu/distances_surface.csv 
python slab2_cross_section_extractor.py --grd car/car_slab2_dep_02.24.18.grd --thk car/car_slab2_thk_02.24.18.grd --eq car/puerto_rico.csv --csv car/distances_surface.csv 
python slab2_cross_section_extractor.py --grd kur/kur_slab2_dep_02.24.18.grd --thk kur/kur_slab2_thk_02.24.18.grd --eq kur/kuril.csv --csv kur/distances_surface.csv 
python slab2_cross_section_extractor.py --grd van/van_slab2_dep_02.23.18.grd --thk van/van_slab2_thk_02.23.18.grd --eq van/navolou.csv --csv van/distances_surface.csv 
python slab2_cross_section_extractor.py --grd sum/sum_slab2_dep_02.23.18.grd --thk sum/sum_slab2_thk_02.23.18.grd --eq sum/java.csv --csv sum/distances_surface.csv 
python slab2_cross_section_extractor.py --grd alu/alu_slab2_dep_02.23.18.grd --thk alu/alu_slab2_thk_02.23.18.grd --eq alu/adak.csv --csv alu/distances_surface.csv 
python slab2_cross_section_extractor.py --grd ker/ker_slab2_dep_02.24.18.grd --thk ker/ker_slab2_thk_02.24.18.grd --eq ker/tonga.csv --csv ker/distances_surface.csv 
python slab2_cross_section_extractor.py --grd phi/phi_slab2_dep_02.26.18.grd --thk phi/phi_slab2_thk_02.26.18.grd --eq phi/philippines.csv --csv phi/distances_surface.csv 
python slab2_cross_section_extractor.py --grd ryu/ryu_slab2_dep_02.26.18.grd --thk ryu/ryu_slab2_thk_02.26.18.grd --eq ryu/okinawa.csv --csv ryu/distances_surface.csv 

================================

eq_plot_v2.py was used to create histograms comparing interface with center for each location at defined depth intervals. Terminal output includes statistical summaries.

python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv sco/distances_surface.csv --title "Central Islands (shallow)" --output charts/histograms/sco_histogram_25_70.png --max -25 --min -70 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv izu/distances_surface.csv --title "Mariana Islands (shallow)" --output charts/histograms/izu_histogram_25_70.png --max -25 --min -70 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv car/distances_surface.csv --title "Puerto Rico (shallow)" --output charts/histograms/car_histogram_25_70.png --max -25 --min -70 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv kur/distances_surface.csv --title "Kuril (shallow)" --output charts/histograms/kur_histogram_25_70.png --max -25 --min -70 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv van/distances_surface.csv --title "Navolou (shallow)" --output charts/histograms/van_histogram_25_70.png --max -25 --min -70 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv sum/distances_surface.csv --title "Java (shallow)" --output charts/histograms/sum_histogram_25_70.png --max -25 --min -70 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv alu/distances_surface.csv --title "Adak (shallow)" --output charts/histograms/alu_histogram_25_70.png --max -25 --min -70 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv ker/distances_surface.csv --title "Tonga (shallow)" --output charts/histograms/ker_histogram_25_70.png --max -25 --min -70 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv phi/distances_surface.csv --title "Philippines (shallow)" --output charts/histograms/phi_histogram_25_70.png --max -25 --min -70 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv ryu/distances_surface.csv --title "Okinawa (shallow)" --output charts/histograms/ryu_histogram_25_70.png --max -25 --min -70 --dpi 600

python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv sco/distances_surface.csv --title "Central Islands (intermediate)" --output charts/histograms/sco_histogram_70_300.png --max -70 --min -300 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv izu/distances_surface.csv --title "Mariana Islands (intermediate)" --output charts/histograms/izu_histogram_70_300.png --max -70 --min -300 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv car/distances_surface.csv --title "Puerto Rico (intermediate)" --output charts/histograms/car_histogram_70_300.png --max -70 --min -300 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv kur/distances_surface.csv --title "Kuril (intermediate)" --output charts/histograms/kur_histogram_70_300.png --max -70 --min -300 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv van/distances_surface.csv --title "Navolou (intermediate)" --output charts/histograms/van_histogram_70_300.png --max -70 --min -300 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv sum/distances_surface.csv --title "Java (intermediate)" --output charts/histograms/sum_histogram_70_300.png --max -70 --min -300 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv alu/distances_surface.csv --title "Adak (intermediate)" --output charts/histograms/alu_histogram_70_300.png --max -70 --min -300 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv ker/distances_surface.csv --title "Tonga (intermediate)" --output charts/histograms/ker_histogram_70_300.png --max -70 --min -300 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv phi/distances_surface.csv --title "Philippines (intermediate)" --output charts/histograms/phi_histogram_70_300.png --max -70 --min -300 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv ryu/distances_surface.csv --title "Okinawa (intermediate)" --output charts/histograms/ryu_histogram_70_300.png --max -70 --min -300 --dpi 600

python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv sco/distances_surface.csv --title "Central Islands (deep)" --output charts/histograms/sco_histogram_300_1000.png --max -300 --min -1000 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv izu/distances_surface.csv --title "Mariana Islands (deep)" --output charts/histograms/izu_histogram_300_1000.png --max -300 --min -1000 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv car/distances_surface.csv --title "Puerto Rico (deep)" --output charts/histograms/car_histogram_300_1000.png --max -300 --min -1000 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv kur/distances_surface.csv --title "Kuril (deep)" --output charts/histograms/kur_histogram_300_1000.png --max -300 --min -1000 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv van/distances_surface.csv --title "Navolou (deep)" --output charts/histograms/van_histogram_300_1000.png --max -300 --min -1000 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv sum/distances_surface.csv --title "Java (deep)" --output charts/histograms/sum_histogram_300_1000.png --max -300 --min -1000 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv alu/distances_surface.csv --title "Adak (deep)" --output charts/histograms/alu_histogram_300_1000.png --max -300 --min -1000 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv ker/distances_surface.csv --title "Tonga (deep)" --output charts/histograms/ker_histogram_300_1000.png --max -300 --min -1000 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv phi/distances_surface.csv --title "Philippines (deep)" --output charts/histograms/phi_histogram_300_1000.png --max -300 --min -1000 --dpi 600
python eq_plot_v2.py --x-axis "Earthquake distances (km)" --csv ryu/distances_surface.csv --title "Okinawa (deep)" --output charts/histograms/ryu_histogram_300_1000.png --max -300 --min -1000 --dpi 600

================================

slab2_cross_section_extractor.py creates chart of cross-section.


python slab2_cross_section_extractor.py --input_dep sco/sco_slab2_dep_02.23.18.grd --input_thk sco/sco_slab2_thk_02.23.18.grd --earthquakes sco/central_islands.csv --orientation w_e  --center_lat -58 --center_lon 333.7 --output_plot charts/cross_sections/central_islands.png --title "Central Islands" --buffer_degrees 1 --dpi 600
python slab2_cross_section_extractor.py --input_dep izu/izu_slab2_dep_02.24.18.grd --input_thk izu/izu_slab2_thk_02.24.18.grd --earthquakes izu/mariannas.csv --orientation w_e  --center_lat 17.2 --center_lon 146 --output_plot charts/cross_sections/mariana_islands.png --title "Mariana Islands" --buffer_degrees 1 --dpi 600 
python slab2_cross_section_extractor.py --input_dep car/car_slab2_dep_02.24.18.grd --input_thk car/car_slab2_thk_02.24.18.grd --earthquakes car/puerto_rico.csv --orientation n_s  --center_lat 18.5 --center_lon 293.6 --output_plot charts/cross_sections/puerto_rico.png --title "Puerto Rico" --buffer_degrees 1 --dpi 600 
python slab2_cross_section_extractor.py --input_dep kur/kur_slab2_dep_02.24.18.grd --input_thk kur/kur_slab2_thk_02.24.18.grd --earthquakes kur/kuril.csv --orientation w_e  --center_lat 46.7 --center_lon 152.7 --output_plot charts/cross_sections/kuril.png --title "Kuril" --buffer_degrees 1 --dpi 600 
python slab2_cross_section_extractor.py --input_dep van/van_slab2_dep_02.23.18.grd --input_thk van/van_slab2_thk_02.23.18.grd --earthquakes van/navolou.csv --orientation w_e  --center_lat -18.8 --center_lon 168.9 --output_plot charts/cross_sections/navolou.png --title "Navolou" --buffer_degrees 1 --dpi 600 
python slab2_cross_section_extractor.py --input_dep sum/sum_slab2_dep_02.23.18.grd --input_thk sum/sum_slab2_thk_02.23.18.grd --earthquakes sum/java.csv --orientation n_s  --center_lat -8.3 --center_lon 110.6 --output_plot charts/cross_sections/java.png --title "Java" --buffer_degrees 1 --dpi 600 
python slab2_cross_section_extractor.py --input_dep alu/alu_slab2_dep_02.23.18.grd --input_thk alu/alu_slab2_thk_02.23.18.grd --earthquakes alu/adak.csv --orientation n_s  --center_lat 51.98 --center_lon 184.1 --output_plot charts/cross_sections/adak.png --title "Adak" --buffer_degrees 1 --dpi 600 
python slab2_cross_section_extractor.py --input_dep ker/ker_slab2_dep_02.24.18.grd --input_thk ker/ker_slab2_thk_02.24.18.grd --earthquakes ker/tonga.csv --orientation w_e  --center_lat -21.9 --center_lon 185.2 --output_plot charts/cross_sections/tonga.png --title "Tonga" --buffer_degrees 1 --dpi 600 
python slab2_cross_section_extractor.py --input_dep phi/phi_slab2_dep_02.26.18.grd --input_thk phi/phi_slab2_thk_02.26.18.grd --earthquakes phi/philippines.csv --orientation w_e  --center_lat 9.2 --center_lon 125.1 --output_plot charts/cross_sections/philippines.png --title "Philippines" --buffer_degrees 1 --dpi 600 
python slab2_cross_section_extractor.py --input_dep ryu/ryu_slab2_dep_02.26.18.grd --input_thk ryu/ryu_slab2_thk_02.26.18.grd --earthquakes ryu/okinawa.csv --orientation w_e  --center_lat 26.3 --center_lon 128 --output_plot charts/cross_sections/okinawa.png --title "Okinawa" --buffer_degrees 1 --dpi 600 


Bonus Cascadia cross-section out to Mount Rainier:
python slab2_cross_section_extractor.py --input_dep cas/cas_slab2_dep_02.24.18.grd --input_thk cas/cas_slab2_thk_02.24.18.grd --earthquakes cas/cascadia.csv --orientation w_e  --center_lat 46.8 --center_lon 237 --output_plot cas/cascadia.png --title "Cascadia" --buffer_degrees 1 --dpi 600 

================================

Boxplot visualization requires earthquake_data_full.csv.
python earthquake_bloxplot_visualization.py

================================

Distance from center scatterplot.

python distance_from_center_scatter_plot.py --title "Adak" --location "Adak" --csv "earthquake_data_full.csv" --output "charts/distance_from_center_scatter/Adak_scatter.png" --dpi 600
python distance_from_center_scatter_plot.py --title "Central Islands" --location "Central Islands" --csv "earthquake_data_full.csv" --output "charts/distance_from_center_scatter/Central_Islands_scatter.png" --dpi 600
python distance_from_center_scatter_plot.py --title "Java" --location "Java" --csv "earthquake_data_full.csv" --output "charts/distance_from_center_scatter/Java_scatter.png" --dpi 600
python distance_from_center_scatter_plot.py --title "Kuril" --location "Kuril" --csv "earthquake_data_full.csv" --output "charts/distance_from_center_scatter/Kuril_scatter.png" --dpi 600
python distance_from_center_scatter_plot.py --title "Mariana Islands" --location "Mariana Islands" --csv "earthquake_data_full.csv" --output "charts/distance_from_center_scatter/Mariana_Islands_scatter.png" --dpi 600
python distance_from_center_scatter_plot.py --title "Navolou" --location "Navolou" --csv "earthquake_data_full.csv" --output "charts/distance_from_center_scatter/Navolou_scatter.png" --dpi 600
python distance_from_center_scatter_plot.py --title "Okinawa" --location "Okinawa" --csv "earthquake_data_full.csv" --output "charts/distance_from_center_scatter/okinawa_scatter.png" --dpi 600
python distance_from_center_scatter_plot.py --title "Philippines" --location "Philippines" --csv "earthquake_data_full.csv" --output "charts/distance_from_center_scatter/Philippines_scatter.png" --dpi 600
python distance_from_center_scatter_plot.py --title "Puerto Rico" --location "Puerto Rico" --csv "earthquake_data_full.csv" --output "charts/distance_from_center_scatter/Puerto_Rico_scatter.png" --dpi 600
python distance_from_center_scatter_plot.py --title "Tonga" --location "Tonga" --csv "earthquake_data_full.csv" --output "charts/distance_from_center_scatter/Tonga_scatter.png" --dpi 600