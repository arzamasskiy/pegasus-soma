# line styles for ColorBrewer YlGnBu
# for use with sequential data
# provides 8 yellow-green-blue colors of increasing saturation
# compatible with gnuplot >=4.2
# author: Anna Schneider

# line styles
set style line 1 lc rgb '#FFFFD9' # very light yellow-green-blue
set style line 2 lc rgb '#EDF8B1' # 
set style line 3 lc rgb '#C7E9B4' # 
set style line 4 lc rgb '#7FCDBB' # light yellow-green-blue
set style line 5 lc rgb '#41B6C4' # 
set style line 6 lc rgb '#1D91C0' # medium yellow-green-blue
set style line 7 lc rgb '#225EA8' #
set style line 8 lc rgb '#0C2C84' # dark yellow-green-blue

# palette
set palette defined (  0 '#0C2C84',\
                      1 '#225EA8',\
                      2 '#1D91C0',\
                      3 '#41B6C4',\
                      4 '#7FCDBB',\
                      5 '#C7E9B4',\
                      6 '#EDF8B1',\
                      7 '#FFFFD9')
