# Part 0: Preliminary

# Part 0.4: Sets the location of the data to be used and where the packages should be put
install.packages("rstudioapi"); library("rstudioapi");
datadir <- paste(dirname(rstudioapi::getSourceEditorContext()$path), "/Dataset", sep = "")
setwd(datadir)
package_loc <- paste(datadir, "lib", sep = "/")
