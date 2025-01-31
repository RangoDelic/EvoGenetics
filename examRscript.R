library(adegenet)
library(tidyverse)
library(dartR)

# reading the raw file
examData<- read.dataframe("F.raw", header=TRUE, sep=" ")
examData

Fgen<-read.PLINK(file="F.raw", map.file = "F.map")
Fgen

Fgen@pop

Fgen@ind.names

gl.compliance.check(Fgen)

nPop(Fgen)
gl.report.hwe(Fgen_fix,subset="each")
gl.report.hwe(Fgen,subset = "each")

FgenMX<-as.matrix(Fgen)
FgenMX<-Fgen[, -c(1:6)]
gl.report.hwe(FgenMX)

str(Fgen)
#renaming SNP columns
colnames(Fgen_df)<-gsub("_.*","", colnames(Fgen_df))
Fgen_fix <- new("genlight", as.matrix(Fgen_df))
FgenPCA<-gl.pcoa(Fgen_fix)

gl.pcoa(FgenMX)
plot(FgenPCA)

#converting all to numeric
Fgen_numeric <- Fgen_df
Fgen_numeric[, -1]<-apply(Fgen_numeric[, -1], 2, as.numeric)

#PCA plot
FgenPCA<-gl.pcoa(Fgen_fix)
gl.compliance.check(Fgen_fix)
gl.pcoa.plot(x=Fgen, glPca=FgenPCA)


#coordinates
coords<-read.table("coordinates2.txt", header=T)
as.data.frame(coords)
head(coordinates_df)
coords<-as.data.frame(coords)

na.omit(coords)
gl.ibd(Fgen, coordinates = coords[,2-3], distance = "Fst")


Fgen_fix@other$latlon <- coords[, c("Latitude", "Longitude")]

print(Fgen_fix@ind.names)
print(Fgen_fix@pop)

pop_mapping <- data.frame(
  Sample_Code = c("SSWB28", "SSWB24", "SSWB42", "SSWB31", "SSWB32", "WB36", 
                  "SSMR01", "MR01", "CA01", "SSWB51", "SSWB44", "SSWB54"),
  Population  = c("DOR", "DOR", "DOR", "LOG", "LOG", "LOG",
                  "MAR", "MAR", "POR", "POR", "POR", "POR")
)
# Create latitude and longitude table for populations
pop_coords <- data.frame(
  Population = c("DOR", "LOG", "MAR", "POR"),
  lat  = c(40.4168, 38.7169, 33.5731, 39.3999),
  long = c(-3.7038, -9.1399, -7.5898, -8.2245)
)
# Extract population code from sample names
individuals <- data.frame(id = Fgen_fix@ind.names)
individuals$Sample_Code <- sub(".*_(SSWB[0-9]+|MR01|CA01|WB36).*", "\\1", individuals$id)
# Merge individuals with population assignments
individuals <- left_join(individuals, pop_mapping, by = "Sample_Code")
# Merge with coordinates
coordinates_df <- left_join(individuals, pop_coords, by = "Population")
# Save to file
write.table(coordinates_df[, c("id", "lat", "long")], "coordinates2.txt", 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
print(head(coordinates_df))

coords<- read.table("coordinates2.txt",header=TRUE,sep="\t")

print(head(coords))

gl.ibd(Fgen, Dgeo="geographic", distance="Fst")
