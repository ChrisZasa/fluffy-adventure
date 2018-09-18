# Color definitions for MTXQC plots

#color cinnamic acid
color_ca = c('black','red','red')
names(color_ca) = c('within','below','above')

#color normalization plot (CA vs. SumOfArea)
color_norm = c('#1a1a1a','tomato3','dodgerblue3')
names(color_norm) = c('within','above','below')

#color linearity check pSIRM
color_linearity = c("#3288bd","#d53e4f","#f46d43","#bababa", "#f7f7f7")
names(color_linearity) = c('linear','below','above','na', 'NaCal')

#color pathways
color_pathway = c('#d53e4f','#f28b5b', '#fee08b','#3288bd','#e4e899','#9bcc93','#fbf6c1')
names(color_pathway) = c('glyc','glut','ppp','tca','other','aa','nucleobase')  