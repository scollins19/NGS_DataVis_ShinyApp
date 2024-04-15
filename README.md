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

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/01d861c6-8d12-4a55-ae7e-65c45aea90e6)


The bigwigs that you have selected will appear underneath.

The **labels** for each bigwig need to be in the same order as the bigwig files (in this case: ```Treated``` followed by ```Untreated```), **separated by a comma with no spaces**.



![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/7b06e988-2785-4afa-ba41-7f23eded2dc5)



#### Input bed files

Bed files (the file that has the positional information) is also selected from the filefinder.

The bed files themselves (the file that has the positional information) need to be **tab separated** and need to have **at least** four columns and **no column names**.

The first three columns corresponds to the chromosome, start, and end positions, and the fourth column is the name of that region. (if you are using a custom file and don't have specific names, just name them 1,2,3,4...etc )

Bed file example:

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/3beb1bb8-19b7-4d0d-8359-8965c86a920d)



The labels for bed files follow the same rules as the bigwig files (i.e., **separated by a comma**)

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/dae9d583-1b43-4983-8169-03ad7f8a397b)


---------------------------------------------------------------------



#### Other input options

*Window size*: number in basepairs to plot e.g. 10000 is equivalent to -/+5kb from the DSB

*bin size*: the size of each in basepairs. A smaller bin size = longer waiting time. As a common rule of thumb, I usually divide my window size by 100 or 200.

*Filename*: The name prefix of the final PDF plots. Try to add as much information as you can.

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/da84d310-3958-42a3-b7f3-8c1f03c39822)


Then you can proceed to click ```Compute Matrix``` and wait until it is done. If the size of the bed file will affect the computational time. 

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/39f018e0-8f22-4c19-88d4-9da6cef3eacd)


-----------------------------------------





#### Output plots - Heatmaps

Once your matrix is finishes calculating, navigate to the ```Heatmap``` tab.

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/01fe77f7-5481-49a7-ac86-c324dc4c2fc6)


***Plot Title***: Title to go above the heatmap.

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/832cbded-ec7f-437b-914a-8bf9e8734b1a)


***Condition to order heatmap rows***: Which bigwig to order the heamap rows by signal to be applied to all other conditions. Needs to one of the bigwig labels that you originally entered when computing the matrix (In this case: ```Treated``` or ```Untreated```)

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/7e8e7d00-d14c-411e-b20c-eb06749fa85b)


***Bigwig order:*** the order which you want the regions to be represented.  Need to be the same as the bigwig file labels originally entered and need to be comma separated with no spaces (In this case, it will be: ```Untreated,Treated```)

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/394230c7-85ee-47e8-8edb-4eef0f646761)


***Bed file order:*** the order which you want the regions to be represented. Need to be the same as the bed file labels originally entered and need to be comma separated with no spaces (In this case, it will be: ```regions_of_interest,control_regions```)

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/da66fa70-3ff3-4ba9-9f71-23f617865d80)

***Colour palette:*** Select any brewer colour palette (https://r-graph-gallery.com/38-rcolorbrewers-palettes.html):

Colour palettes you can use:

<img src="https://r-graph-gallery.com/38-rcolorbrewers-palettes_files/figure-html/thecode-1.png" alt="R Color Brewer's palettes – the R Graph Gallery" style="zoom: 33%;" />



***Switch Colour Palette direction***: Option to change the direction of the colour palette. Useful when using a diverging colour scale like "RdBu" which goes from red to white to blue. Changing the direction of this palette will force the colours to go from blue to white to red instead.

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/7f34aa9b-8c7b-4d6e-acfc-13289f5a4caf)




***Scale max/min value:*** Can be left as default or you can input any number to force the colour scale. Typing 'Default' will return the colour scale to it's original values.

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/c4846cad-a7c0-41c4-8c3e-79a861919dd5)




Then select ```Make Heatmap```.

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/a938fe9b-cd36-4886-a331-4d9a1ad33fa7)



If you change the colour, order, or scale of the heatmap, you need to select ```Make Heatmap``` each time.

If you are happy with your plot, you can download the heatmap by changing the desired width and height in the Download options and then selecting ```Download Heatmap```.

The downloaded file will already have the pre-filled prefix that you have chosen, along with the window size and what type of plot it is.

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/f62b7b19-f19f-475c-bf42-cf13f26eefc3)


-------------------------------------------



#### Output plots - Average Profiles and boxplots

The Average Profile and Boxplot tabs have similar input options as the Heatmap. Additional options include:

***Y-axis label:*** Change the label that will be shown on the y-axis of the plot

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/5c08ebeb-6ad6-4318-a663-ac8d63d7ac8a)


***Select Layout:*** This option allows you to decide whether you want to group the plot by condition, or by bed file. By default, the plot colours different conditions, and then splits the display by bed file, like so:

![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/673ed192-8dae-4fab-96f3-b9e06bac3e31)




If we change the ```Select Layout``` option to ```variable```, we get the following plot which is now coloured by bed file, and split by condition.



![image](https://github.com/scollins19/NGS_DataVis_ShinyApp/assets/44778109/58665a56-e36b-4a57-ae55-2a78c4fcec83)




The options are the same for both average profiles and boxplots.

When comparing 2 conditions and 1 bed file, the outputted p-val is a paired Wilcox test (aka. Wilcoxon signed-rank test)

When comparing 1 condition and 2 bed files, the outputted p-val is an unpaired Wilcox test (aka. Wilcoxon rank sum test, aka Mann-Whitney test)
