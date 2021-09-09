#########################################################################################
## Joint KFF-Stanford Project to Estimate COVID-19 Vaccination Coverage by Race/Ethnicity
## Contact: mreitsma@stanford.edu
#########################################################################################

rm(list = ls())

## Load packages
library(data.table)
library(ggplot2)
library(DBI)
library(RSQLite)
library(lubridate)
library(zoo)

setwd("~/COVID/vaccination_analysis/KFF Analysis/dashboard/")

## State FIPS codes for merging
fips <- as.data.table(tidycensus::fips_codes)
fips <- unique(fips[,.(state_name, state_code)])
fips <- fips[, state_code:=as.numeric(state_code)]

#### Prepare population data from 2019 Single-Year ACS
df <- fread("./data/usa_00026.csv")
df <- df[,.(SERIAL, STATEFIP, GQ, PERNUM, PERWT, AGE, RACE, RACED, HISPAN, HISPAND, PUMA)]
df <- df[GQ%in%c(1, 2, 5)] ## Remove institutionalized populations
df <- df[AGE>=12] ## Restrict to eligible population

## Assign race/ethnicity
df <- df[, ethnicity:=ifelse(HISPAN%in%c(1,2,3,4), "Hispanic", "Not Hispanic")]
df <- df[RACE==1, race:="White"]
df <- df[RACE==2, race:="Black"]
df <- df[RACE==3, race:="American Indian or Alaska Native"]
df <- df[RACE%in%c(4,5), race:="Asian"]
df <- df[RACE%in%c(8, 9), race:="Other"]
df <- df[RACE%in%c(7), race:="NEC"]
df <- df[RACED%in%c(630, 680, 681, 682, 685, 689, 690, 699), race:="Native Hawaiian or Other Pacific Islander"]
df <- df[is.na(race), race:="Asian"]

## Account for race == "NEC" based on census approach
race_breakdowns <- copy(df)
race_breakdowns <- race_breakdowns[race!="NEC"]
race_breakdowns <- race_breakdowns[, tot_race_eth:=sum(PERWT, na.rm=T), by = c("ethnicity", "race", "STATEFIP")]
race_breakdowns <- race_breakdowns[, tot_race:=sum(PERWT, na.rm=T), by = c("ethnicity", "STATEFIP")]
race_breakdowns <- race_breakdowns[, pct_race:=tot_race_eth/tot_race]
race_breakdowns <- unique(race_breakdowns[,.(STATEFIP, ethnicity, race, pct_race)])

nec_expand <- df[race=="NEC"]
nec_expand <- nec_expand[, race:=NULL]
nec_expand <- merge(nec_expand, race_breakdowns, by = c("STATEFIP", "ethnicity"), allow.cartesian=T)
nec_expand <- nec_expand[, PERWT:=PERWT*pct_race]

df <- df[race!="NEC"]
df <- rbind(df, nec_expand[, pct_race:=NULL])

df <- df[ethnicity=="Hispanic", race_grp:="Hispanic"]
df <- df[is.na(race_grp) & race == "White", race_grp:="White"]
df <- df[is.na(race_grp) & race == "Black", race_grp:="Black"]
df <- df[is.na(race_grp) & race == "American Indian or Alaska Native", race_grp:="American Indian or Alaska Native"]
df <- df[is.na(race_grp) & race == "Asian", race_grp:="Asian"]
df <- df[is.na(race_grp) & race == "Other", race_grp:="Other"]
df <- df[is.na(race_grp) & race == "Native Hawaiian or Other Pacific Islander", race_grp:="Native Hawaiian or Other Pacific Islander"]

## Hispanic population breakdown by race
df <- merge(df, fips, by.y = "state_code", by.x="STATEFIP")
hisp_breakdown <- copy(df)
hisp_breakdown <- hisp_breakdown[, tot_race_eth:=sum(PERWT, na.rm=T), by = c("race", "ethnicity", "state_name")]
hisp_breakdown <- hisp_breakdown[, tot_eth:=sum(PERWT, na.rm=T), by = c("ethnicity", "state_name")]
hisp_breakdown <- hisp_breakdown[, tot_race:=sum(PERWT, na.rm=T), by = c("race", "state_name")]
hisp_breakdown <- hisp_breakdown[, pct_race:=tot_race_eth/tot_eth]
hisp_breakdown <- hisp_breakdown[, pct_eth:=tot_race_eth/tot_race]

race_eth <- copy(hisp_breakdown)
race_eth <- unique(race_eth[,.(race, ethnicity, state_name, tot_eth, tot_race)])
race_eth <- race_eth[, pct_eth:=tot_eth/sum(tot_eth, na.rm=T), by = c("state_name", "race")]
race_eth <- race_eth[, pct_race:=tot_race/sum(tot_race, na.rm=T), by = c("state_name", "ethnicity")]
race_eth <- race_eth[ethnicity == "Hispanic"]
setnames(race_eth, "race", "race_grp")
race_eth <- race_eth[, .(state_name, race_grp, ethnicity, pct_eth, pct_race)]
race_eth <- melt(race_eth, id.vars = c("race_grp", "state_name", "ethnicity"), value.var = c("pct_eth", "pct_race"))
race_eth <- race_eth[variable=="pct_eth", race_grp:="Hispanic"]
race_eth <- race_eth[, c("variable", "ethnicity"):=NULL]
setnames(race_eth, "value", "pct_pop")
race_eth <- unique(race_eth)

race_breakdown_by_ethnicity <- copy(hisp_breakdown)
race_breakdown_by_ethnicity <- unique(race_breakdown_by_ethnicity[, .(tot_race_eth, state_name, race, ethnicity)])
race_breakdown_by_ethnicity <- race_breakdown_by_ethnicity[, pct_race_eth:=tot_race_eth/sum(tot_race_eth, na.rm=T), by = c("state_name")]

hisp_breakdown <- hisp_breakdown[ethnicity=="Hispanic"]
hisp_breakdown <- unique(hisp_breakdown[,.(race, state_name, pct_race)])
setnames(hisp_breakdown, "race", "race_grp")

## Merge state name
df <- df[AGE>=12, pop_12_allrace:=sum(PERWT, na.rm=T), by = c("state_name")]
df <- df[AGE>=12, pop_12:=sum(PERWT, na.rm=T), by = c("race_grp", "state_name")]
df <- df[!is.na(pop_12)]
df <- unique(df[,.(pop_12, pop_12_allrace, state_name, race_grp)])

#### Read in KFF extractions
kff_full <- NULL
for (i in c("03-01-2021", "03-15-2021", "03-29-2021", "04-05-2021", "04-12-2021",
            "04-19-2021", "04-26-2021", "05-03-2021", "05-10-2021",
            "05-17-2021", "05-24-2021", "06-07-2021", "06-14-2021", "06-21-2021",
            "06-28-2021", "07-06-2021", "07-19-2021", "08-02-2021",
            "08-16-2021")) {
  kff <- as.data.table(readxl::read_xlsx("./data/Vaccine Distribution 3_15-6_7.xlsx", sheet = paste0(i)))
  setnames(kff, "...1", "state_name")
  kff <- kff[!is.na(state_name)]
  kff <- kff[, c("Total People Vaccinated", "Total People Vaccinated with known Race data", "Total People Vaccinated with known ethnicity"):=NULL]
  num_cols <- colnames(kff)[colnames(kff)%like%"%"]
  kff <- kff[, c(num_cols):=lapply(.SD, as.numeric), .SDcols = c(num_cols), with = F]
  kff <- kff[, c("% of Vaccinations with Known Race", "% of Vaccinations with Known Ethnicity"):=NULL]
  
  kff <- melt(kff, id.vars = c("state_name", "Race Categories Include Hispanic Individuals"))
  kff <- kff[, race_grp:=tstrsplit(x = variable, split = " % of Vaccinations")]
  kff <- kff[, variable:=NULL]
  kff <- kff[race_grp=="Hispanic" & is.na(value), flip_toggle:=1]
  kff <- kff[, flip_toggle:=mean(flip_toggle, na.rm=T), by = "state_name"]
  kff <- kff[flip_toggle==1, `Race Categories Include Hispanic Individuals`:=NA]
  kff <- kff[, flip_toggle:=NULL]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp=="Hispanic", value_adj:=value]
  
  kff <- kff[race_grp!="% of Vaccinations with Unknown Race" & race_grp!="% of Vaccinations with Unknown Ethnicity"]
  kff <- kff[race_grp=="Other", value:=NA] # Other unrealistically high, consider these data unknown for computing uptake rates
  
  kff <- merge(kff, race_eth, by = c("race_grp", "state_name"), all.x=T)
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp!="Hispanic", value:=value/sum(value, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp=="Hispanic", fixed_hisp:=value_adj]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", fixed_hisp:=mean(fixed_hisp, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp!="Hispanic" & !is.na(value), pct_pop:=pct_pop/sum(pct_pop, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", rr_race:=value/pct_pop]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & race_grp=="Hispanic", rr_eth:=rr_race]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", rr_hisp:=mean(rr_eth, na.rm=T), by = "state_name"]
  
  kff <- merge(kff, race_breakdown_by_ethnicity[ethnicity=="Hispanic"], by.x = c("state_name", "race_grp"), by.y=c("state_name", "race"), all.x=T)
  
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & ethnicity=="Hispanic", temp_joint:=pct_race_eth*rr_hisp*rr_race]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & ethnicity=="Hispanic", temp_joint_hisp_sum:=sum(temp_joint, na.rm=T), by = "state_name"]
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes", final_joint:=temp_joint/(temp_joint_hisp_sum/fixed_hisp)]
  
  kff <- kff[`Race Categories Include Hispanic Individuals`=="Yes" & !is.na(final_joint), value_adj:=value-final_joint]
  
  kff <- kff[`Race Categories Include Hispanic Individuals`=="" | is.na(`Race Categories Include Hispanic Individuals`), value_adj:=value/sum(value, na.rm=T), by = "state_name"]
  kff <- kff[, Date:=mdy(paste0(i))]
  kff_full <- rbind(kff_full, kff, fill = T)
}

kff_full <- merge(kff_full, df, by = c("race_grp", "state_name"), all=T)
kff_full <- kff_full[!is.na(value_adj), pop_pct:=pop_12/sum(pop_12, na.rm=T), by = c("state_name", "Date")]

#### Read in and clean CDC state vaccination time series
vax_stats <- fread("./data/COVID-19_Vaccinations_in_the_United_States_Jurisdiction.csv") #https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-Jurisdi/unsk-b7fc
vax_stats <- vax_stats[, Date:=mdy(Date)]
setnames(vax_stats, "Location", "state")
fips <- as.data.table(tidycensus::fips_codes)
fips <- unique(fips[,.(state_name, state, state_code)])
vax_stats <- merge(vax_stats, fips, by = c("state"))

# Impute 12+ for full time series
vax_stats <- vax_stats[, pop12:=Administered_Dose1_Recip_12Plus/(Administered_Dose1_Recip_12PlusPop_Pct/100)]
vax_stats <- vax_stats[, pop12:=mean(pop12, na.rm=T), by = "state_name"]
vax_stats <- vax_stats[Administered_Dose1_Recip_12Plus==0, Administered_Dose1_Recip_12Plus:=Administered_Dose1_Recip]
vax_stats <- vax_stats[Administered_Dose1_Recip_12PlusPop_Pct==0, Administered_Dose1_Recip_12PlusPop_Pct:=(Administered_Dose1_Recip_12Plus/pop12)*100]

vax_stats <- vax_stats[,.(Date, state_name, Administered_Dose1_Recip_12Plus, Administered_Dose1_Recip_12PlusPop_Pct)]
vax_stats <- vax_stats[, lapply(.SD, as.numeric), by = c("state_name", "Date")]
vax_stats <- merge(vax_stats, unique(kff_full[,.(state_name, pop_12_allrace)]), by = "state_name")
vax_stats <- vax_stats[state_name=="Idaho", Administered_Dose1_Recip_12PlusPop_Pct:=(Administered_Dose1_Recip_12Plus/pop_12_allrace)*100]

## IMPUTE DISCONTINUITIES IN TOTAL VAX ##
vax_stats <- vax_stats[, doses_raw:=Administered_Dose1_Recip_12PlusPop_Pct]

for (i in as.list(seq.Date(as.Date(min(vax_stats$Date)), as.Date(max(vax_stats$Date)) , "days"))) {
  print(ymd(i))
  vax_stats <- vax_stats[, current_temp:=NULL]
  vax_stats <- vax_stats[Date == ymd(i), current_temp:=doses_raw]
  vax_stats <- vax_stats[, current_temp:=mean(current_temp, na.rm=T), by = c("state_name")]
  vax_stats <- vax_stats[Date < ymd(i) & doses_raw > current_temp, doses_raw:=NA]
}

vax_stats <- vax_stats[, row_id:=.I]

out <- NULL
for (s in unique(vax_stats$state_name)) {
  temp <- vax_stats[state_name==s]
  for (i in as.list(seq.Date(as.Date(min(vax_stats$Date)), as.Date(max(vax_stats$Date)) , "days"))) {
    if (is.na(temp$doses_raw[temp$Date==i])) {
      print(paste0("Interpolating: ", s, ", ", i))
      
      start_date <- max(temp$Date[temp$Date<i & !is.na(temp$doses_raw)])
      end_date <- min(temp$Date[temp$Date>i & !is.na(temp$doses_raw)])
      
      if (is.na(nchar(ymd(start_date))) | is.na(nchar(ymd(end_date)))) {
        print("Cannot interpolate due to lack of start/end data")
      } else {
        ## Compute growth factor
        b <- 1-((temp$pop_12_allrace[temp$Date==end_date]-temp$doses_raw[temp$Date==end_date])/(temp$pop_12_allrace[temp$Date==start_date]-temp$doses_raw[temp$Date==start_date]))^(1/(interval(ymd(start_date),ymd(end_date))/days(1)))
        start_val <- temp$pop_12_allrace[temp$Date==start_date]-temp$doses_raw[temp$Date==start_date]
        temp <- temp[Date == i, interpolated:=pop_12_allrace-(start_val*(1-b)^(interval(ymd(start_date),ymd(i))/days(1)))]
        rm(start_val, b)
      }
      rm(start_date, end_date)
    }
    out <- rbind(out, temp[Date==i], fill = T)
  }
}

ggplot(data = out[state_name == "New Hampshire" & Date >= "2021-04-15"], aes(x = Date)) +
  geom_point(aes(y = doses_raw), color = "red") + geom_point(aes(y = interpolated), color = "blue") +
  geom_line(data = out[state_name == "New Hampshire" & !is.na(doses_raw) & Date >= "2021-04-15"], aes(y = doses_raw))

out <- out[is.na(doses_raw), doses_raw:=interpolated]
vax_stats <- copy(out)
vax_stats <- vax_stats[,Administered_Dose1_Recip_12PlusPop_Pct:=doses_raw]
vax_stats <- vax_stats[,.(state_name, Date, Administered_Dose1_Recip_12PlusPop_Pct)]

vax_stats <- merge(vax_stats, kff_full, by = c("state_name", "Date"))
vax_stats <- vax_stats[, pop_share_agg:=pop_12/sum(pop_12, na.rm=T), by = c("state_name", "Date")]
vax_stats <- vax_stats[is.na(value_adj), missing:=1]
vax_stats <- vax_stats[!is.na(value_adj), missing:=0]
vax_stats <- vax_stats[is.na(value_adj), pct_missing:=pop_share_agg]
vax_stats <- vax_stats[, pct_missing:=sum(pct_missing, na.rm=T), by = c("state_name", "Date")]

vax_stats <- vax_stats[, value_adj:=value_adj*(1-pct_missing)]
vax_stats <- vax_stats[is.na(value_adj), value_adj:=pop_share_agg]
vax_stats <- vax_stats[, value_adj:=value_adj/sum(value_adj, na.rm=T), by = c("state_name", "Date")]
vax_stats <- vax_stats[, doses_adj:=(Administered_Dose1_Recip_12PlusPop_Pct/100)*pop_12_allrace]
vax_stats <- vax_stats[, race_doses:=doses_adj*value_adj]

## Deal with 100+ round 1
vax_stats <- vax_stats[race_doses>pop_12, add_doses:=race_doses-pop_12]
vax_stats <- vax_stats[, add_doses:=sum(add_doses, na.rm=T), by = c("state_name", "Date")]
vax_stats <- vax_stats[race_doses>pop_12, exceeded:=1]
vax_stats <- vax_stats[race_doses>pop_12, race_doses:=pop_12]
vax_stats <- vax_stats[is.na(exceeded), second_pct:=value_adj/sum(value_adj), by = c("state_name", "Date")]
vax_stats <- vax_stats[is.na(exceeded), race_doses:=race_doses+(second_pct*add_doses)]

vax_stats <- vax_stats[, vaccinated_pct_12:=(race_doses/pop_12)*100]

pdf("./figures/historical_coverage_8.16.21.pdf", width = 8, height = 6)
for (i in unique(kff_full$state_name)) {
  plot <- ggplot(data = vax_stats[state_name==i & race_grp%in%c("Hispanic", "White", "Black", "Asian") & missing!=1], aes(x = Date, y = vaccinated_pct_12, color = race_grp)) + 
    geom_line(data = vax_stats[state_name==i & race_grp=="White"], aes(y = (doses_adj/pop_12_allrace)*100), size = 1.2, color = "dark gray") +
    geom_point(size = 2) + geom_line( size = 1.3) +
    theme_bw() +
    scale_color_manual(breaks= c("Hispanic", "White", "Black", "Asian"), values = c("All" = "#767676", "Hispanic" = "#c42e31", "White" = "#832543", "Asian" = "#6399AC", "Black" = "#e5a825")) +
    labs(x = "", y = "Coverage", color = "Race/Ethnicity", title = paste0(i)) +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 70, linetype = "dashed")+
    scale_y_continuous(limits = c(0, 100))
  print(plot)
}
dev.off()

## PER UNVACCINATED ##
backup <- copy(vax_stats)
vax_stats <- copy(backup)

## SMOOTH DISCONTINUITIES ##

## STEP 1: Fill
vax_stats <- vax_stats[missing!=1, race_doses_approx:=race_doses]

for (i in as.list(seq.Date(as.Date(min(vax_stats$Date)), as.Date(max(vax_stats$Date)) , "days"))) {
  print(ymd(i))
  vax_stats <- vax_stats[, current_temp:=NULL]
  vax_stats <- vax_stats[Date == ymd(i), current_temp:=race_doses_approx]
  vax_stats <- vax_stats[, current_temp:=mean(current_temp, na.rm=T), by = c("state_name", "race_grp")]
  # vax_stats <- vax_stats[Date < ymd(i) & race_doses_approx > current_temp, race_doses_approx:=current_temp]
  vax_stats <- vax_stats[Date < ymd(i) & race_doses_approx > current_temp, race_doses_approx:=NA]
}

vax_stats <- vax_stats[, row_id:=.I]

out <- NULL
for (s in unique(vax_stats$state_name)) {
  temp <- vax_stats[state_name==s & !is.na(race_doses_approx)] # This makes sure that the race is not always missing for the state
  for (r in unique(temp$race_grp)) {
    temp <- vax_stats[state_name==s & race_grp==r]
    for (i in c("2021-03-15", "2021-03-29", "2021-04-05", "2021-04-12",
                "2021-04-19", "2021-04-26", "2021-05-03", "2021-05-10",
                "2021-05-17", "2021-05-24", "2021-06-07", "2021-06-14", "2021-06-21",
                "2021-06-28", "2021-07-06", "2021-07-19", "2021-08-02",
                "2021-08-16")) {
      if (is.na(temp$race_doses_approx[temp$Date==i])) {
        print(paste0("Interpolating: ", r, ", ", s, ", ", i))
      
        start_date <- max(temp$Date[temp$Date<i & !is.na(temp$race_doses_approx)])
        end_date <- min(temp$Date[temp$Date>i & !is.na(temp$race_doses_approx)])

        if (is.na(nchar(ymd(start_date))) | is.na(nchar(ymd(end_date)))) {
          print("Cannot interpolate due to lack of start/end data")
        } else {
          ## Compute growth factor
            b <- 1-((temp$pop_12[temp$Date==end_date]-temp$race_doses_approx[temp$Date==end_date])/(temp$pop_12[temp$Date==start_date]-temp$race_doses_approx[temp$Date==start_date]))^(1/(interval(ymd(start_date),ymd(end_date))/days(1)))
            start_val <- temp$pop_12[temp$Date==start_date]-temp$race_doses_approx[temp$Date==start_date]
            temp <- temp[Date == i, interpolated:=pop_12-(start_val*(1-b)^(interval(ymd(start_date),ymd(i))/days(1)))]
            rm(start_val, b)
        }
        rm(start_date, end_date)
      }
      out <- rbind(out, temp[Date==i], fill = T)
    }
  }
}

ggplot(data = out[race_grp=="Asian" & state_name == "Alaska"], aes(x = Date)) +
  geom_point(aes(y = race_doses_approx), color = "red") + geom_point(aes(y = interpolated), color = "blue") +
  geom_line(data = out[race_grp=="Asian" & state_name == "Alaska" & !is.na(race_doses_approx)], aes(y = race_doses_approx))

out <- rbind(out, vax_stats[!(row_id %in% unique(out$row_id))], fill = T)
out <- out[, interpolated_indicator:=ifelse(!is.na(interpolated), 1, 0)]
out <- out[is.na(race_doses_approx), race_doses_approx:=interpolated]

out <- out[is.na(race_doses_approx) & !is.na(value), missing:=1]
out <- out[order(state_name, race_grp, Date)]
out <- out[, race_doses_approx_fill:=race_doses_approx]
out <- out[, race_doses_approx_fill:=nafill(race_doses_approx_fill, type = "nocb"), by = c("state_name", "race_grp")]
out <- out[is.na(race_doses_approx) & !is.na(value), race_doses_approx:=race_doses_approx_fill]

out <- out[, vaccinated_pct_12:=(race_doses_approx/pop_12)*100]

pdf("./figures/historical_coverage_8.16.21_adj.pdf", width = 8, height = 6)
for (i in unique(kff_full$state_name)) {
  plot <- ggplot(data = out[state_name==i & race_grp%in%c("Hispanic", "White", "Black", "Asian") & missing!=1], aes(x = Date, y = vaccinated_pct_12, color = race_grp)) + 
    geom_line(data = out[state_name==i & race_grp=="White"], aes(y = (doses_adj/pop_12_allrace)*100), size = 1.2, color = "dark gray") +
    geom_line( size = 1.3) +
    geom_point(size = 3, aes(shape = as.factor(interpolated_indicator))) +
    theme_bw() +
    scale_color_manual(breaks= c("Hispanic", "White", "Black", "Asian"), values = c("All" = "#767676", "Hispanic" = "#c42e31", "White" = "#832543", "Asian" = "#6399AC", "Black" = "#e5a825")) +
    labs(x = "", y = "Coverage", color = "Race/Ethnicity", title = paste0(i)) +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 70, linetype = "dashed")+
    scale_y_continuous(limits = c(0, 100))
  print(plot)
}
dev.off()

vax_stats <- copy(out)

vax_stats <- vax_stats[is.na(race_doses_approx), race_doses_approx:=race_doses]
vax_stats <- vax_stats[order(state_name, race_grp, Date)]
vax_stats <- vax_stats[, l2:=shift(race_doses_approx, n = 2, type = "lag"), by = c("state_name", "race_grp")]
vax_stats <- vax_stats[, diff:=race_doses_approx-l2]
vax_stats <- vax_stats[, diff_rate:=(diff/((pop_12-l2)*28))]

# vax_stats <- vax_stats[, all_race:=doses_adj]
# vax_stats <- vax_stats[, l2_all:=shift(all_race, n = 1, type = "lag"), by = c("state_name", "race_grp")]
# vax_stats <- vax_stats[, diff_all:=all_race-l2_all]
# vax_stats <- vax_stats[, diff_all_rate:=(diff_all/(pop_12_allrace*28))]
# vax_stats_out <- rbind(vax_stats_out, vax_stats)
# vax_stats <- copy(vax_stats_out)

## Edge cases
vax_stats <- vax_stats[diff_rate<0, diff_rate:=0]
# vax_stats <- vax_stats[diff_all_rate<0, diff_all_rate:=0]

backup <- copy(vax_stats)
vax_stats <- copy(backup)

for (i in ymd("2021-08-16"):ymd("2021-12-31")) {
  projection <- vax_stats[Date==i]
  # Per unvaccinated population
  projection <- projection[(pop_12-race_doses_approx)<=0, race_doses_approx:=race_doses_approx]
  projection <- projection[(pop_12-race_doses_approx)>0, race_doses_approx:=race_doses_approx+(diff_rate*(pop_12-race_doses_approx))]
  # projection <- projection[, doses_adj:=doses_adj+(diff_all_rate*(pop_12_allrace))]
  
  projection <- projection[, Date:=Date+1]
  projection <- projection[, project:=1]
  vax_stats <- rbind(vax_stats, projection, fill = T)
}

vax_stats <- vax_stats[, vaccinated_pct_12:=(race_doses_approx/pop_12)*100]
vax_stats <- vax_stats[vaccinated_pct_12>100 | is.na(vaccinated_pct_12), race_doses_approx:=pop_12]
vax_stats <- vax_stats[vaccinated_pct_12>100 | is.na(vaccinated_pct_12), vaccinated_pct_12:=100]

vax_stats <- vax_stats[, tot_doses:=sum(race_doses_approx, na.rm=T), by = c("state_name", "Date")]
vax_stats <- vax_stats[, agg:=(tot_doses/pop_12_allrace)*100]

write.csv(vax_stats, "~/Downloads/vax_stats.csv", na = "", row.names = F)

vax_stats_out <- copy(vax_stats)
vax_stats_out <- vax_stats_out[race_grp%in%c("Asian", "Black", "Hispanic", "White")]

temp <- copy(vax_stats_out)
temp <- temp[, vaccinated_pct_12:=agg]
temp <- temp[, pop_12:=pop_12_allrace]
temp <- temp[, race_grp:="All"]
temp <- unique(temp[,.(state_name, Date, race_grp, vaccinated_pct_12, pop_12, project)])

vax_stats_out <- vax_stats_out[,.(state_name, Date, race_grp, vaccinated_pct_12, pop_12, project, missing)]
vax_stats_out <- rbind(vax_stats_out, temp, fill = T)

vax_stats_out <- vax_stats_out[, vaccinated_pct_12_out:=as.character(round(vaccinated_pct_12, digits = 2))]
vax_stats_out <- vax_stats_out[missing == 1, vaccinated_pct_12_out:="NR"]
vax_stats_out <- vax_stats_out[, missing:=NULL]
vax_stats_out <- vax_stats_out[, vaccinated_pct_12:=NULL]

vax_stats_out <- merge(vax_stats_out, fips, by = "state_name")

vax_stats_nat <- copy(vax_stats)
vax_stats_nat <- vax_stats_nat[, nat_doses_race:=sum(race_doses_approx, na.rm=T), by = c("race_grp", "Date")]
vax_stats_nat <- vax_stats_nat[, nat_pop_race:=sum(pop_12, na.rm=T), by = c("race_grp", "Date")]
vax_stats_nat <- vax_stats_nat[, nat_doses_agg:=sum(race_doses_approx, na.rm=T), by = c("Date")]
vax_stats_nat <- vax_stats_nat[, nat_pop_agg:=sum(pop_12, na.rm=T), by = c("Date")]

vax_stats_nat <- vax_stats_nat[, nat_race:=(nat_doses_race/nat_pop_race)*100]
vax_stats_nat <- vax_stats_nat[, nat_race_agg:=(nat_doses_agg/nat_pop_agg)*100]

vax_stats_nat <- unique(vax_stats_nat[,.(Date, race_grp, nat_race_agg, nat_race, nat_pop_agg, nat_pop_race, project)])

setnames(vax_stats_nat, c("nat_race", "nat_race_agg", "nat_pop_race", "nat_pop_agg"),
         c("vaccinated_pct_12", "agg", "pop_12", "pop_12_allrace"))

temp <- copy(vax_stats_nat)
temp <- temp[, vaccinated_pct_12:=agg]
temp <- temp[, pop_12:=pop_12_allrace]
temp <- temp[, race_grp:="All"]
temp <- unique(temp[,.(Date, race_grp, vaccinated_pct_12, pop_12, project)])

vax_stats_nat <- rbind(vax_stats_nat, temp, fill = T)
vax_stats_nat <- vax_stats_nat[race_grp%in%c("Asian", "Black", "Hispanic", "White", "All")]
vax_stats_nat <- vax_stats_nat[, state_name:="United States"]
vax_stats_nat <- vax_stats_nat[, state_code:="00"]
vax_stats_nat <- vax_stats_nat[, state:="US"]

vax_stats_nat <- vax_stats_nat[,.(state_name, Date, race_grp, pop_12, project, vaccinated_pct_12, state, state_code)]
vax_stats_nat <- vax_stats_nat[, vaccinated_pct_12_out:=as.character(round(vaccinated_pct_12, digits = 2))]
vax_stats_nat <- vax_stats_nat[, vaccinated_pct_12:=NULL]

vax_stats_out <- rbind(vax_stats_out, vax_stats_nat)

write.csv(vax_stats_out, "./output/coverage_time_series.csv", na = "", row.names = F)

suppress_table <- copy(vax_stats)
suppress_table <- suppress_table[Date=="2021-08-16"]
suppress_table <- suppress_table[, nr:=ifelse(missing == 1, 1, 0)]
suppress_table <- unique(suppress_table[,.(state_name, race_grp, nr)])
suppress_table <- suppress_table[race_grp%in%c("Asian", "Black", "Hispanic", "White")]

vax_date_out <- copy(vax_stats)
vax_date_out <- vax_date_out[,.(state_name, Date, race_grp, vaccinated_pct_12)]
vax_date_out <- vax_date_out[race_grp%in%c("Asian", "Black", "Hispanic", "White")]
vax_dates_square <- as.data.table(expand.grid(Date = seq(ymd("2021-03-01"), ymd("2021-12-31"), by = "day"), race_grp = c("Asian", "Black", "Hispanic", "White"), state_name = unique(vax_date_out$state_name)))
vax_date_out <- merge(vax_date_out, vax_dates_square, by = c("Date", "race_grp", "state_name"), all = T)
vax_date_out <- vax_date_out[order(state_name, race_grp, Date)]
vax_date_out <- vax_date_out[, vaccinated_pct_12_approx:=na.approx(vaccinated_pct_12), by = c("state_name", "race_grp")]

vax_date_out <- vax_date_out[vaccinated_pct_12_approx>=80, date_:=min(Date, na.rm=T), by = c("state_name", "race_grp")]
vax_date_out <- vax_date_out[, date_:=mean(date_, na.rm=T), by = c("state_name", "race_grp")]
vax_date_out <- unique(vax_date_out[,.(state_name, race_grp, date_)])
vax_date_out <- vax_date_out[is.na(date_), date_:=ymd("2022-01-01")]
vax_date_out <- vax_date_out[, date_70:=month(date_, label = T)]
vax_date_out <- vax_date_out[date_<"2021-09-01", date_70:="Reached 80%"]
vax_date_out <- vax_date_out[date_70 == "Sep", date_70:="Sep. 2021"]
vax_date_out <- vax_date_out[date_70 == "Oct", date_70:="Oct. 2021"]
vax_date_out <- vax_date_out[date_70 == "Nov", date_70:="Nov. 2021"]
vax_date_out <- vax_date_out[date_70 == "Dec", date_70:="Dec. 2021"]
vax_date_out <- vax_date_out[date_70 == "Jan", date_70:="Jan. 2022 or later"]
vax_date_out <- merge(vax_date_out, suppress_table, by = c("state_name", "race_grp"), all = T)
vax_date_out <- vax_date_out[nr==1, date_70:="NR"]
vax_date_out <- vax_date_out[, nr:=NULL]
vax_date_out <- vax_date_out[, date_:=NULL]

vax_date_out <- merge(vax_date_out, fips, by = "state_name")
setnames(vax_date_out, "date_70", "date_80")

write.csv(vax_date_out, "./output/map_date_80.csv", na = "", row.names = F)

pdf("./figures/projections_state_8.16.21.pdf", width = 8, height = 6)
for (i in unique(kff_full$state_name)) {
  plot <- ggplot(data = vax_stats[state_name==i & race_grp%in%c("Hispanic", "White", "Black", "Asian") & is.na(project) & missing==0],
                 aes(x = Date, y = vaccinated_pct_12, color = race_grp)) + 
    geom_line(data = vax_stats[state_name==i & race_grp=="White" & is.na(project)], aes(y = agg), size = 1, color = "dark gray") +
    geom_point(size = 2) + geom_line( size = 1.1) +
    geom_line(data = vax_stats[state_name==i & race_grp=="White" & (project == 1 | Date == "2021-08-16")], aes(y = agg), size = 1.2, color = "dark gray",
              linetype = "dashed", alpha = 0.9) +
    geom_line(data = vax_stats[state_name==i & race_grp%in%c("Hispanic", "White", "Black", "Asian") & (project == 1 | Date == "2021-08-16") & missing==0],
              size = 1.3,
              linetype = "dashed", alpha = 0.9) +
    theme_bw() +
    scale_color_manual(breaks= c("Hispanic", "White", "Black", "Asian"),
                       values = c("All" = "#767676", "Hispanic" = "#c42e31", "White" = "#832543", "Asian" = "#6399AC", "Black" = "#e5a825")) +
    labs(x = "", y = "Coverage Among Eligible (12+)", color = "Race/Ethnicity", title = paste0(i),
         caption = "Gray line shows coverage among whole eligible (12+) population.\nSolid lines show historical (observed) data, dashed lines show projections.") +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 80, linetype = "dashed") +
    scale_y_continuous(limits = c(0, 100)) +
    scale_shape_manual(values = c(19, 1))
  print(plot)
}
dev.off()
