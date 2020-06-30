library("SKAT")
library("argparse")

# create parser object
parser <- ArgumentParser()

parser$add_argument("--bfile", help="Prefix of PLINK binary files.")
parser$add_argument("--cfile", help="Covariate file.")
parser$add_argument("--mfile", help="Mapping file.")
parser$add_argument("--opref", help="Output prefix.")
parser$add_argument("--ncov", help="Number of covariates to use.")

args <- parser$parse_args()

start_time = Sys.time()
cat(paste('... Execution started at ', start_time))

# Set the pointers to the data files.
File.Bed <- paste(args$bfile, "bed", sep=".")
File.Bim <- paste(args$bfile, "bim", sep=".")
File.Fam <- paste(args$bfile, "fam", sep=".")
File.Cov <- args$cfile
File.SetID <- args$mfile

# Files generated during execution.
File.SSD <- paste0(args$opref, ".SSD")
File.Info <- paste0(args$opref,".SSD.info")
File.Out <- paste0(args$opref, "_results.txt")
File.Ent <- paste0(args$opref, "_eff_no_tests.txt")
File.Time <- paste0(args$opref, "_profiling.txt")

# Generate the File.SSD file and SSD.Info file.
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, 
                   File.Info)

# Read the FAM and COV files.
FAM <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary=FALSE, cov_header=TRUE)
y <- FAM$Phenotype

SSD.INFO<-Open_SSD(File.SSD, File.Info)

#Z <-Get_Genotypes_SSD(SSD.INFO, SSD.INFO$SetInfo)

# ------------------------------------------------------------------------------
# create the NULL model
if (args$ncov == 0)
{
  cat('Using no covariates.')
  obj<-SKAT_Null_Model(Phenotype ~ 1, out_type="C", data=FAM, Adjustment=FALSE)
} else if (args$ncov == 1)
{
  cat('Using 1 covariate.')
  obj<-SKAT_Null_Model(Phenotype ~ PC1, out_type="C", data=FAM, 
                       Adjustment=FALSE)
} else if (args$ncov == 2)
{
  cat('Using 2 covariates.')
  obj<-SKAT_Null_Model(Phenotype ~ PC1 + PC2, out_type="C", data=FAM, 
                       Adjustment=FALSE)
} else if (args$ncov == 3)
{
  cat('Using 3 covariates.')
  obj<-SKAT_Null_Model(Phenotype ~ PC1 + PC2 + PC3, out_type="C", data=FAM, 
                       Adjustment=FALSE)
} else if (args$ncov == 10) {

  cat('Using 10 covariates.')
  obj<-SKAT_Null_Model(Phenotype ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, out_type="C", data=FAM,
                       Adjustment=FALSE)
}

# ------------------------------------------------------------------------------
# Test the results.
out.skat<-SKAT.SSD.All(SSD.INFO, obj, method="SKATO")
write.table(as.data.frame(out.skat$results),
            file=File.Out, quote=F, sep="\t", row.names=F)

# ------------------------------------------------------------------------------
# Intermediate profiling (testing).
tmp_stamp = Sys.time()
elapsed_test = difftime(time1=tmp_stamp, time2=start_time, units="secs")
cat(paste('... Testing finished after ', elapsed_test, ' secs.'))

# Total time since start.
total_time = difftime(time1=Sys.time(), time2=start_time, units="secs")

# create the output string.
out_str = c(paste("Total:", total_time, "sec."), 
            paste("\t Testing:", elapsed_test, "sec."))

# write profiling to the file.
fout<-file(File.Time)
writeLines(out_str, fout)
close(fout)


# ------------------------------------------------------------------------------
# Compute the number of effective tests.
eff_no_tests = Get_EffectiveNumberTest(out.skat$results$MAP, alpha=0.05)
write.table(as.data.frame(eff_no_tests),
            file=File.Ent, quote=F, sep="\t", row.names=F)


# ------------------------------------------------------------------------------
# Compute the complete runtime.
elapsed_eff = difftime(time1=Sys.time(), time2=tmp_stamp, units="secs")
cat(paste('... Eff. no. tests finished after ', elapsed_eff, ' secs.'))

total_time = difftime(time1=Sys.time(), time2=start_time, units="secs")

# create the output string.
out_str = c(paste("Total:", total_time, "sec."), 
            paste("\t Testing:", elapsed_test, "sec."), 
            paste("\t Eff. number of tests:", elapsed_eff, "sec." ))

fout<-file(File.Time)
writeLines(out_str, fout)
close(fout)


