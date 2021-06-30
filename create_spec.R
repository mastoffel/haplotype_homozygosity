# make a alpha peel spec file
library(glue)

create_spec <- function(spec, nsnps, genos, ped, out) {
        
        out <- glue("
         nSnp           , {nsnps}
         InputFilePath  , {genos}
         pedigree       , {ped}
         OutputFilePath , {out}
         nCycles        , 20
         runType        , multi
         startsnp       , 500
         endsnp         , 520
             ")
        
        write.lines(out, spec)
}
