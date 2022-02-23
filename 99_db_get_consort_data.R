library(RJDBC)
library(dplyr)
library(dbplyr)
library(lubridate)
library(tidyverse)
library(janitor)

dbname <- "../sheep/data/db/StKilda_Data.accdb"
driver <- "net.ucanaccess.jdbc.UcanloadDriver"
driverpath <- "../sheep/data/db/UCanAccess/loader/ucanload.jar"
options <- paste0("jdbc:ucanaccess://", dbname, ";memory=false")

con <- DBI::dbConnect(JDBC(driver, driverpath), options)
# src <- src_dbi(con)

tbls <- dbGetTables(con)
flds <- dbGetFields(con, "tblPregnancies")
consorts <- dbGetQuery(con, "Select * from Consorts")
pregs <- dbGetQuery(con, "Select * from tblPregnancies")
dbDisconnect(con)

write_delim(consorts, "data/consorts.txt", " ")
write_delim(pregs, "data/pregnancies.txt", " ")
