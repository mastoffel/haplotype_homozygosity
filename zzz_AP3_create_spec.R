# make a alpha peel spec file
library(glue)

create_spec <- function(spec, nsnps, genos, ped, out, startsnp, endsnp) {
        
        out <- glue("
         nSnp           , {nsnps}
         InputFilePath  , {genos}
         pedigree       , {ped}
         OutputFilePath , {out}
         nCycles        , 10
         runType        , multi
         startsnp       , {startsnp}
         endsnp         , {endsnp}
         
             ")
        
        writeLines(out, spec)
}
