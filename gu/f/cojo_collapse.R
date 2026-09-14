args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=7L)stop('Expected helper file, input, output, size, CHR, POS, P')
# Load only the requested function; sourcing the shared file would reset the RNG.
e<-new.env(parent=baseenv())
for(expr in parse(args[1])) if(is.call(expr) && identical(expr[[1]],as.name('<-')) && identical(expr[[2]],as.name('loci_collapse'))) eval(expr,e)
if(!exists('loci_collapse',e,inherits=FALSE))stop('loci_collapse missing from shared helper')
environment(e$loci_collapse)<-globalenv()
x<-e$loci_collapse(args[2],args[4],CHR=args[5],POS=args[6],P=args[7],output=args[3])
cat(nrow(x$members),'signals ->',nrow(x$loci),'non-overlapping source loci\n')
