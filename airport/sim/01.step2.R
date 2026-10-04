# 🚩 01.step2
pacman::p_load(data.table, tidyverse, countrycode, lubridate, rworldmap)

dir = "D:/data/avia"

airport <- read.csv(paste0(dir, "/iata-icao.csv"), header = TRUE) %>%
	mutate(iso3c = countrycode(country_code, "iso2c", "iso3c")) %>% dplyr::select( - c(country_code, region_name, icao))
guests <- read.csv(paste0(dir, "/OAG/2019_01.csv"), header = T) %>%
	filter(Time.series >= as.Date("2019-01-21") & Time.series <= as.Date("2019-01-27")) %>%
	merge(airport, by.x = "Origin", by.y = "iata") %>% merge(airport, by.x = "Destination", by.y = "iata") %>%
	filter(iso3c.x != "CHN" & iso3c.y == "CHN") %>% rename(iso3c = iso3c.x) %>%
	dplyr::select(iso3c, Seats) %>% group_by(iso3c) %>% summarize(Seats = sum(Seats))
pop <- subset(world_bank_pop, indicator == "SP.POP.TOTL", select = c("country", "2017")) %>% rename(iso3c = country, pop = "2017")
unrate <- read.table(paste0(dir, "/cfr.txt"), header = TRUE) %>% mutate(estimate = as.numeric(estimate), date = as.Date(date))


dat.raw <- data.table::fread(paste0(dir, "/covid19_confirmed_global.csv")) %>%
	dplyr::select( - c("Province/State", "Lat", "Long")) %>% rename("country" = "Country/Region") %>%
	melt(., id.vars = "country", variable.name = "date", value.name = "new_cases") %>% mutate(date := mdy(date)) %>%
	group_by(date, country) %>% summarize(new_cases = sum(as.numeric(new_cases))) %>%
	#  此处不能用: group_by(country) %>% mutate(new_cases = new_cases - lag(new_cases))
	mutate(iso3c := countrycode(country, 'country.name', 'iso3c', custom_match = c("Micronesia" = "FSM", "Kosovo" = "XXK"))) %>%
	merge(unrate, by = c("iso3c", "date"), all.x = TRUE) %>% merge(pop, by = "iso3c", all.x = TRUE) %>%
	mutate(truerate = (new_cases / estimate) / pop, week = ceiling(as.numeric(date - as.Date("2020-01-19")) / 7)) %>%
	dplyr::select(iso3c, truerate, week) %>%
	group_by(iso3c, week) %>% summarize(truerate = mean(truerate, na.rm = TRUE)) %>%
	merge(guests, by = "iso3c", all.x = TRUE) %>% mutate(cases = round(Seats * truerate)) # merge 都要用 all.x=TRUE


dat1 <- subset(dat.raw, week == 1); sum(dat1$Seats, na.rm = TRUE)
dat1 <- dat.raw %>% filter(week == 10) %>% drop_na()
	sPDF <- joinCountryData2Map(dat1, joinCode = "ISO3", nameJoinColumn = "iso3c")
	mapPlot <- mapCountryData(sPDF, nameColumnToPlot = "truerate", catMethod = "quantiles", colourPalette = "heat", addLegend = TRUE)
dat1 <- dat.raw %>% filter(week >= 10, week <= 150)
	sum(dat1$Seats, na.rm = TRUE) * 0.8256; sum(dat1$cases, na.rm = TRUE) * 0.8256
	# 2020年1月19日星期日作为起点，1月22日开始有JHU感染数据，3月29日星期日【第10周】开始实施“五个一”【2023年1月8日结束】。 
	# week<=37【2020年9月底】，week<=67【2021年5月初】，week<=102【2021年12月底】，week<=150【2022年11月底】


