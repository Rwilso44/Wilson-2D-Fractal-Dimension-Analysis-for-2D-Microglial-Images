setwd("/Users/rebeccawilson/Desktop/Microglia Morphology Short Communication/All Clean Masks for Analysis")

# Load the functions (full path, so this works from anywhere)
source("/Users/rebeccawilson/Desktop/Microglia Morphology Short Communication/All Clean Masks for Analysis/boxcount.R")

# Point R at the folder with your masks (outputs are saved here too)
setwd("/Users/rebeccawilson/Desktop/Microglia Morphology Short Communication/Known Images")

# Check R can see the file
#my_file <- "Clean_Mask of C1-C03L16 P21 M CD Sld16 AM Composite Full-Microglia 2.tif"
#file.exists(my_file)          # should print TRUE

# Run everything for this image
for (f in list.files(pattern = "\\.tif$")) analyze_image(f)
analyze_image(my_file)