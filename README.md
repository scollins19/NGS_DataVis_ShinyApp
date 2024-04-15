# NGS_DataVis_ShinyApp

### READ ME 

### R-Shiny app for plotting heatmaps, average profiles and boxplots of bigwig files

A user-friendly R-shiny app to visualise NGS datasets (bigwig files) using the ```computematrix``` functionality of ```deeptools```.

The app computes coverage matrices of  one or more bigwig files over specific regions from one or multiple bed files. The matrices are then imported to R as a dataframe to be used with ``ggplot2``.

The app has two functionalities, plotting NGS data centered on specific regions, or plotting over gene bodies where the length of each gene is divided into an equal number of bins. In the below sections, I will run through the functions of the app, the format of input files, and most certainly what not to do.



-----------------------------------------------------------------------

#### Requirements 

Deeptools (https://deeptools.readthedocs.io/en/develop/content/installation.html). Recommended installation with ```conda```

```
conda install -c bioconda deeptools
```



R packages:

``` 
require(stringi)
require(shiny)
require("dplyr")
require("ggplot2")
require("reshape2")
require("rtracklayer")
require("BSgenome.Hsapiens.UCSC.hg19")
require("plyranges")
require("tidyverse")
require("RColorBrewer")
require("ggpubr")
require(rstatix)
require(gridExtra)
require(cowplot)
require("shinyjs")
require(shinyWidgets)
require(EnsDb.Hsapiens.v75)
require(stringr)
require(shinythemes)
require(shinyFiles)
```



---------------------------------------------------------------------------------------------------------------

#### Installation tips

To successfully run the app, you will need to change some lines of the configuration file ```deeptools_functions.r``` to include the paths of your working directory, the path to 

```
## PATH OF WORKING DIRECTORY
working_dir <- "/home/scollins/Heatmaps_and_profiles/"
```

```
## PATH TO `COMPUTE MATRIX` FUNCTION OF DEEPTOOLS
deeptoolspath <- "/home/DATA_HOME/anaconda3/envs/deeptools/bin/computeMatrix"
```

```
## VOLUMES TO MOUNT WITH RSHINY APP 
volumes <- c(DATA_HOME='/home/DATA/')
```





-------------------------------------------

## Using the app:

### Plots centred on specific regions (DSBs)

#### Input bigwigs:

Input bigwigs can be directly selected using the file finder.

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415141825075.png" alt="image-20240415141825075" style="zoom:67%;" />

The bigwigs that you have selected will appear underneath.

The **labels** for each bigwig need to be in the same order as the bigwig files (in this case: ```Treated``` followed by ```Untreated```), **separated by a comma with no spaces**.



<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415142024388.png" alt="image-20240415142024388" style="zoom:67%;" />



#### Input bed files

Bed files (the file that has the positional information) is also selected from the filefinder.

The bed files themselves (the file that has the positional information) need to be **tab separated** and need to have **at least** four columns and **no column names**.

The first three columns corresponds to the chromosome, start, and end positions, and the fourth column is the name of that region. (if you are using a custom file and don't have specific names, just name them 1,2,3,4...etc )

For the example of the DSB bed file:

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240227165530586.png" alt="image-20240227165530586" style="zoom:50%;" />



The labels for bed files follow the same rules as the bigwig files (i.e., **separated by a comma**)

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415142149598.png" alt="image-20240415142149598" style="zoom:67%;" />

---------------------------------------------------------------------



#### Other input options

*Window size*: number in basepairs to plot e.g. 10000 is equivalent to -/+5kb from the DSB

*bin size*: the size of each in basepairs. A smaller bin size = longer waiting time. As a common rule of thumb, I usually divide my window size by 100 or 200.

*Filename*: The name prefix of the final PDF plots. Try to add as much information as you can.

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415142300345.png" alt="image-20240415142300345" style="zoom:67%;" />

Then you can proceed to click ```Compute Matrix``` and wait until it is done. If the size of the bed file will affect the computational time. 

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240228141131918.png" alt="image-20240228141131918" style="zoom:50%;" />

-----------------------------------------





#### Output plots - Heatmaps

Once your matrix is finishes calculating, navigate to the ```Heatmap``` tab.

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240228141949885.png" alt="image-20240228141949885" style="zoom: 67%;" />

***Plot Title***: Title to go above the heatmap.

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415142620128.png" alt="image-20240415142620128" style="zoom:67%;" />

***Condition to order heatmap rows***: Which bigwig to order the heamap rows by signal to be applied to all other conditions. Needs to one of the bigwig labels that you originally entered when computing the matrix (In this case: ```Treated``` or ```Untreated```)

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415142720732.png" alt="image-20240415142720732" style="zoom:67%;" />

***Bigwig order:*** the order which you want the regions to be represented.  Need to be the same as the bigwig file labels originally entered and need to be comma separated with no spaces (In this case, it will be: ```Untreated,Treated```)

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415142804148.png" alt="image-20240415142804148" style="zoom:67%;" />

***Bed file order:*** the order which you want the regions to be represented. Need to be the same as the bed file labels originally entered and need to be comma separated with no spaces (In this case, it will be: ```regions_of_interest,control_regions```)

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415142920166.png" alt="image-20240415142920166" style="zoom:67%;" />

***Colour palette:*** Select any brewer colour palette (https://r-graph-gallery.com/38-rcolorbrewers-palettes.html):

Colour palettes you can use:

<img src="https://r-graph-gallery.com/38-rcolorbrewers-palettes_files/figure-html/thecode-1.png" alt="R Color Brewer's palettes – the R Graph Gallery" style="zoom: 33%;" />



***Switch Colour Palette direction***: Option to change the direction of the colour palette. Useful when using a diverging colour scale like "RdBu" which goes from red to white to blue. Changing the direction of this palette will force the colours to go from blue to white to red instead.

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415112009229.png" alt="image-20240415112009229" style="zoom:67%;" />



***Scale max/min value:*** Can be left as default or you can input any number to force the colour scale. Typing 'Default' will return the colour scale to it's original values.

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415112217069.png" alt="image-20240415112217069" style="zoom:67%;" />



Then select ```Make Heatmap```.

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415143123994.png" alt="image-20240415143123994" style="zoom:67%;" />



If you change the colour, order, or scale of the heatmap, you need to select ```Make Heatmap``` each time.

If you are happy with your plot, you can download the heatmap by changing the desired width and height in the Download options and then selecting ```Download Heatmap```.

The downloaded file will already have the pre-filled prefix that you have chosen, along with the window size and what type of plot it is.

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415143155837.png" alt="image-20240415143155837" style="zoom:67%;" />

-------------------------------------------



#### Output plots - Average Profiles and boxplots

The Average Profile and Boxplot tabs have similar input options as the Heatmap. Additional options include:

***Y-axis label:*** Change the label that will be shown on the y-axis of the plot

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415114027970.png" alt="image-20240415114027970" style="zoom:67%;" />

***Select Layout:*** This option allows you to decide whether you want to group the plot by condition, or by bed file. By default, the plot colours different conditions, and then splits the display by bed file, like so:

<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415143251125.png" alt="image-20240415143251125" style="zoom:67%;" />



If we change the ```Select Layout``` option to ```variable```, we get the following plot which is now coloured by bed file, and split by condition.



<img src="/home/scollins/.config/Typora/typora-user-images/image-20240415143313474.png" alt="image-20240415143313474" style="zoom:67%;" />





The options are the same for both average profiles and boxplots.

When comparing 2 conditions and 1 bed file, the outputted p-val is a paired Wilcox test (aka. Wilcoxon signed-rank test)

When comparing 1 condition and 2 bed files, the outputted p-val is an unpaired Wilcox test (aka. Wilcoxon rank sum test, aka Mann-Whitney test)
