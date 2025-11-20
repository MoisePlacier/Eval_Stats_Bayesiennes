# Function to save mcmc draws (use the boa/coda packages)
# buildMCMC : Function of BRugs - build an "mcmc.list" class object
# as.matrix : Function of coda - transform an "mcmc.list" class object in a "matrix" class object

f.save.mcmc <- function(var, write)
{

x <- as.data.frame(as.matrix(buildMCMC(var)))
name_file = var
if(eval(write))
{
write.table( x, file = name_file, 
             append = FALSE, row.names = TRUE, col.names = TRUE, quote = FALSE )
}
else x=x
}
