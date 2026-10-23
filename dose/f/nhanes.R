#!/usr/bin/env Rscript
# nhanes: score, patterns, stroke, plot, flow
# Launch with dose.sh; source this file to register definitions for reuse.
if (!exists("dose_step_code", mode="function")) {
    source(file.path(Sys.getenv("DOSE_SCRIPT_DIR", if (.Platform$OS.type == "windows")
        "D:/scripts/dose/f" else "/mnt/d/scripts/dose/f"), "0.common.R"), local=TRUE)
}

# ===========================================================================
# nhanes / score
# ===========================================================================
dose_register_step("nhanes", "score", expression({
# NHANES external replication of the CURRENT UKB DOSE definition.
# Intermediate cache: .cache/nhanes. Existing scoring is reused unless dose.nhanes.rebuild=TRUE.
# Selected primary: 2003-2018, both reliable days; day 1 is a sensitivity analysis.
# Current DOSE9half; sample sizes are derived after eligibility, not fixed targets.
# Two-day DOSE = score(mean intake), NOT mean(daily scores).
# Mapping uses the current UKB item portions and standalone NHANES food codes.
# Recipe ingredients cannot be recovered from whole-dish grams: see mapping audit.

suppressPackageStartupMessages({library(dplyr); library(haven)})
NH_ROOT <- dose_output_dir()
NH_DATA <- Sys.getenv("NHANES_DATA", dose_files_dir())
NH_OUT <- file.path(NH_ROOT, ".cache/nhanes")
NH_FPED_DIR <- Sys.getenv("NHANES_FPED_DIR", file.path(dose_files_dir(), "fped"))
nh_usda_file <- function(filename, required=TRUE) {
    direct <- file.path(NH_FPED_DIR, filename)
    if (file.exists(direct)) return(direct)
    paths <- list.files(NH_FPED_DIR, recursive=TRUE, full.names=TRUE)
    paths <- paths[tolower(basename(paths)) == tolower(filename)]
    if (length(paths) == 1L) return(paths[[1L]])
    if (!required && !length(paths)) return(character())
    stop("Expected one USDA file ", filename, " under ", NH_FPED_DIR,
         "; found ",length(paths),".")
}
nh_required_usda <- function() c(
    as.vector(outer(c("fped_dr1iff_", "fped_dr2iff_", "fped_dr1tot_", "fped_dr2tot_"),
                    c("0506", "0708", "0910", "1112", "1314", "1516", "1718"),
                    function(prefix, year) paste0(prefix, year, ".sas7bdat"))),
    "equiv9400.txt", "equiv0102.txt", "equiv0304.sas7bdat")
# Only WJFRT comes from the package; CDC, FPED and MPED still use local source data.
nh_wjfrt <- local({
    reference <- NULL
    function() {
        if (is.null(reference)) {
            if (!requireNamespace("dietaryindex", quietly=TRUE))
                stop("Install a compatible dietaryindex package providing NHANES_20032004$WJFRT.")
            ref_env <- new.env(parent=emptyenv())
            utils::data("NHANES_20032004", package="dietaryindex", envir=ref_env)
            candidate <- ref_env$NHANES_20032004$WJFRT
            if (!is.data.frame(candidate) || !nrow(candidate) ||
                !all(c("FoodCode", "ModCode", "WHOLEFRT", "FRTJUICE") %in% names(candidate)))
                stop("dietaryindex must provide a compatible NHANES_20032004$WJFRT data frame.")
            reference <<- candidate
        }
        reference
    }
})
nh_check_usda <- function() {
    paths <- list.files(NH_FPED_DIR, recursive=TRUE, full.names=TRUE)
    required <- nh_required_usda()
    missing <- required[!tolower(required) %in% tolower(basename(paths))]
    if (length(missing)) stop("Missing ", length(missing), " USDA input files under ",
        NH_FPED_DIR, ":\n", paste(missing, collapse="\n"))
    invisible(vapply(required, nh_usda_file, character(1)))
}
NH_CYCLES <- data.frame(cycle=paste(seq(1999,2017,2),seq(2000,2018,2),sep="-"),
                        suffix=c("", paste0("_",LETTERS[2:10])))
nh_scoring <- new.env(parent=globalenv())
local({
    old <- options(dose.score.run=FALSE); on.exit(options(old))
    dose_source_step("score", "build", envir=nh_scoring)
})
nh_csv <- function(x, name) if(isTRUE(getOption("dose.nhanes.audit",FALSE))) write.csv(x,file.path(NH_OUT,name),row.names=FALSE,na="")
nh_xpt_path <- function(cycle,stem) {
    filename <- paste0(stem,NH_CYCLES$suffix[match(cycle,NH_CYCLES$cycle)],".xpt")
    path <- file.path(NH_DATA,"cdc",cycle,filename)
    if (file.exists(path)) return(path)
    candidates <- list.files(file.path(NH_DATA,"raw",cycle), recursive=TRUE, full.names=TRUE)
    candidates <- candidates[tolower(basename(candidates)) == tolower(filename)]
    if (length(candidates) != 1L)
        stop("Expected one NHANES input for ",cycle,"/",filename,"; found ",length(candidates))
    candidates[[1L]]
}
nh_xpt <- function(cycle,stem) {
    path <- nh_xpt_path(cycle,stem)
    x <- as.data.frame(haven::read_xpt(path)); names(x) <- toupper(names(x))
    # Older XPT descriptions contain Windows-1252 bytes marked as UTF-8.
    for (nm in names(x)) if(is.character(x[[nm]])) {
        bad <- !is.na(x[[nm]]) & !validUTF8(x[[nm]])
        x[[nm]][bad] <- iconv(x[[nm]][bad],from="windows-1252",to="UTF-8")
    }
    x
}
nh_col <- function(x, candidates, required=TRUE) {
    hit <- intersect(candidates,names(x))
    if (!length(hit)) {
        if (required) stop("Missing column: ",paste(candidates,collapse=" / "))
        return(rep(NA_real_,nrow(x)))
    }
    as.numeric(x[[hit[1]]])
}
nh_binary <- function(x) ifelse(x==1,1L,ifelse(x==2,0L,NA_integer_))

# Covariate definitions aligned to ukb/00data/analysis_config.rds.
# Retain raw/proxy labels: NHANES cannot exactly reproduce every UKB instrument.
nh_aligned_covariates <- function(cycle) {
    demo<-nh_xpt(cycle,"DEMO"); alq<-nh_xpt(cycle,"ALQ"); paq<-nh_xpt(cycle,"PAQ")
    bpq<-nh_xpt(cycle,"BPQ"); diq<-nh_xpt(cycle,"DIQ")
    get<-function(x,nm) nh_col(x,nm,FALSE)[match(demo$SEQN,x$SEQN)]
    yn<-function(x) factor(nh_binary(x),0:1,c("No","Yes"))
    n<-nrow(demo); alcohol<-rep(NA_character_,n)
    if(cycle>="2017") {
        ever<-get(alq,"ALQ111"); freq<-get(alq,"ALQ121")
        alcohol[ever %in% 2]<-"never"
        alcohol[ever %in% 1 & freq %in% 0]<-"previous"
        alcohol[freq %in% 1:10]<-"current"
    } else {
        # ALD100 is the 2001-2002 name. Lifetime <12 includes rare drinkers.
        annual<-get(alq,c("ALQ101","ALQ100","ALD100")); lifetime<-get(alq,"ALQ110")
        freq<-get(alq,"ALQ120Q")
        alcohol[annual %in% 2 & lifetime %in% 2]<-"never"
        alcohol[(annual %in% 1 | lifetime %in% 1) & freq %in% 0]<-"previous"
        alcohol[is.finite(freq)&freq>0&freq<777]<-"current"
    }
    # Moderate leisure activity is the common observable domain across 1999-2018.
    # Before 2007: past 30 days. From 2007: typical week, with observed days 1-7.
    pa<-get(paq,if(cycle<"2007") "PAD320" else "PAQ665")
    pa_days<-rep(NA_real_,n)
    if(cycle>="2007") {
        pd<-get(paq,"PAQ670");pa_days[pa %in% 2]<-0
        pa_days[pa %in% 1 & pd %in% 1:7]<-pd[pa %in% 1 & pd %in% 1:7]
    }
    pa_cat<-ifelse(pa %in% 1,"Yes",ifelse(pa %in% 2:3,"No","Unknown"))
    # Only documented negative skip paths imply no medication; refusal stays NA.
    med<-function(now,negative) {
        ans<-nh_binary(now)
        ans[is.na(now)&negative]<-0L
        factor(ans,0:1,c("No","Yes"))
    }
    bp<-get(bpq,"BPQ020");bp_pres<-get(bpq,"BPQ040A")
    chol_pres<-get(bpq,"BPQ090D");chol<-get(bpq,"BPQ080");checked<-get(bpq,"BPQ060")
    dm<-get(diq,"DIQ010");insulin<-get(diq,"DIQ050");pills<-get(diq,c("DIQ070","DID070"))
    mi<-nh_binary(insulin);mp<-nh_binary(pills)
    mi[is.na(insulin)&dm %in% 2]<-0L;mp[is.na(pills)&dm %in% 2]<-0L
    drugdm<-ifelse(mi %in% 1 | mp %in% 1,1L,ifelse(mi %in% 0 & mp %in% 0,0L,NA_integer_))
    edu<-get(demo,"DMDEDUC2");eth<-get(demo,"RIDRETH1")
    data.frame(SEQN=demo$SEQN,cycle,
      ethnicity_ukb=factor(ifelse(eth %in% 3,"White",ifelse(eth %in% c(1,2,4,5),"Others",NA)),c("White","Others")),
      education_ukb=factor(ifelse(edu %in% 1:2,"low",ifelse(edu %in% 3:4,"middle",ifelse(edu %in% 5,"high",NA))),c("low","middle","high")),
      alcohol_status=factor(alcohol,c("never","previous","current")),
      pa_moderate=factor(pa_cat,c("No","Yes","Unknown")),pa_moderate_days=pa_days,
      drug_htn=med(get(bpq,"BPQ050A"),bp %in% 2 | bp_pres %in% 2),
      drug_lipid=med(get(bpq,"BPQ100D"),chol_pres %in% 2 | (cycle<"2011" & chol %in% 2) | checked %in% 2),
      # UKB drug.dm uses code 3 in Fields 6153/6177: INSULIN, not all diabetes drugs.
      drug_dm=factor(mi,0:1,c("No","Yes")),
      drug_diabetes_any=factor(drugdm,0:1,c("No","Yes")))
}

nh_food_rules <- function() {
    # First matching rule wins; each code maps to at most one current UKB item.
    # Portion is ALWAYS looked up from current food_group_map.csv.
    r <- function(item,prefix,include,exclude="") data.frame(item,prefix,include,exclude)
    rules <- bind_rows(
      r("p104450","6","^apples?([, ]|$)","juice|sauce|butter|dried|dehydrated|candied|baby"),
      r("p104560","6","^pears?([, ]|$)","juice|sauce|dried|dehydrated|candied|baby"),
      r("p104420","6","^prunes?([, ]|$)","juice|baby"),
      r("p104430","6","^(raisins?|dates?|figs?|apricots?|apples?|pears?|peaches?|cranberries|mango|papaya|fruit)","juice|sauce|baby|candied"),
      r("p104350","[78]","^tomatoes?([, ]|$).*canned","sauce|paste|juice|soup|with.*(okra|corn)|dried|puree"),
      r("p104340","[78]","^tomatoes?([, ]|$)","sauce|paste|juice|soup|with.*(okra|corn)|dried|puree|canned"),
      r("p104140","7","^broccoli([, ]|$)","soup|casserole|salad|with.*(rice|noodles|pasta)|and cauliflower"),
      r("p104160","7","^(cabbage|kale)([, ]|$)","soup|coleslaw|salad|roll|stuffed|pickled"),
      r("p104240","[78]","^lettuce([, ]|$)","salad|with.*(tomato|meat|chicken|egg|cheese)"),
      r("p104300","[78]","^spinach([, ]|$)","soup|salad|dip|casserole|souffle"),
      r("p104310","7","^brussels sprouts?([, ]|$)"),
      r("p104370","7","^watercress([, ]|$)","soup|salad"),
      r("p103270","4","^tofu([, ]|$)","soup|with.*(rice|noodles)|sandwich"),
      r("p104000","4","^(beans?.*baked|baked beans?)", "soup"),
      r("p104120","7","^(beans?,.*(green|string|snap|wax)|green beans?|string beans?)", "soup|salad|casserole"),
      r("p104280","7","^(peas?,.*green|green peas?|peas?, cooked|peas?, raw|peas?, canned|peas?, frozen)","soup|salad|carrot"),
      r("p104110","[47]","^(broad beans?|fava beans?|beans?,.*(broad|fava))", "soup|salad"),
      r("p104010","4","^(beans?|refried beans?|lentils?|chickpeas?|chick peas?|garbanzo|cowpeas?|black.?eyed peas?|split peas?|soybeans?)([, ]|$)","soup|salad|chips|snack|flour|raw|dry,|uncooked|sprouts|with meat"),
      r("p26131","1","^milk,.*(reduced fat|2%|semi.?skim)","dry|evaporated|condensed|chocolate|flavor|malt|cocoa"),
      r("p26133","1","^milk,.*(low.?fat|fat.?free|nonfat|skim|1%)", "dry|evaporated|condensed|chocolate|flavor|malt|cocoa"),
      r("p26102","1","^yogh?urt,.*(low.?fat|nonfat|fat.?free|skim|reduced fat)","frozen|dip|drink"),
      r("p102870","1","^(cheese, cottage|cottage cheese)"),
      r("p102850","1","^cheese.*(spread|cream).*", "frosting|dip"),
      r("p102810","1","^cheese,", "sauce|sandwich|spread|cream|cottage"),
      r("p102950","3","^egg.*omelet", "sandwich"),
      r("p102980","3","^egg.*(scrambled|substitute|white|deviled|pickled)","sandwich|salad"),
      r("p102940","3","^egg[s]?,", "omelet|scrambled|substitute|white|deviled|pickled|salad|sandwich|dried"),
      r("p102960","[35]","^egg.*sandwich", "bacon|sausage|ham"),
      r("p102970","[23]","^scotch egg"),
      r("p103070","[28]","^(bacon|pork bacon|turkey bacon)([, ]|$)","sandwich|bits|with.*(egg|lettuce)"),
      r("p103080","[28]","^(ham|turkey ham)([, ]|$)","fresh|baby|sandwich|salad"),
      r("p103010","2","(sausage|frankfurter|hot dog)","sandwich|bun|biscuit|pizza|pasta|rice|with.*egg|casserole|corn dog"),
      r("p26111","8","^(margarine|vegetable oil.*spread).*", "butter blend"),
      r("p26112","8","^(margarine|vegetable oil.*spread).*", "butter blend"),
      r("p100220","[169]","^(fruit smoothie|smoothie, fruit)","dry|powder"),
      r("p100230","[19]","^(smoothie|dairy smoothie|yogurt.*smoothie)","dry|powder"),
      r("p100550","[19]","^(hot chocolate|cocoa,.*(beverage|made|prepared|ready))", "dry, unprepared|powder, unprepared"),
      r("p100530","1","^(milk,.*(chocolate|flavor|strawberry)|chocolate milk|milk shake|milkshake|milk, malted)","dry|powder|condensed"),
      r("p100170","9","^(soft drink|carbonated beverage|soda,|fruit flavored drink, carbonated)","concentrate|dry|powder"),
      r("p100180","[69]","^(fruit (drink|punch|juice drink|flavored drink)|juice drink|lemonade|limeade|cranberry juice (cocktail|drink)|orange juice drink|apple juice beverage)","dry|powder|concentrate"),
      r("p26110","8","^(olive oil|oil, olive)([, ]|$)"),
      r("p26151","9","^(wine.*(dessert|fortified|port|sherry)|port wine|sherry)","nonalcoholic|cooler|spritzer"),
      r("p26152","9","^wine.*red","nonalcoholic|cooler|spritzer"),
      r("p26153","9","^wine","nonalcoholic|cooler|spritzer|rice|dessert|fortified|port|sherry|red"),
      r("p100770","562","^(oatmeal|porridge)","dry|uncooked"),
      r("p100800","57","muesli"),
      r("p100810","57","granola|oat crunch"),
      r("p100840","57","bran"),
      r("p100850","57","shredded wheat|whole wheat|whole-wheat"),
      r("p100820","57","sweetened|frosted|sugar|chocolate|cocoa|honey|marshmallow|fruit loops|froot|cinnamon"),
      r("p100830","57","corn flakes|cornflakes|puffed rice|rice crisp|rice krisp|cheerios|puffed wheat"),
      r("p100860","57",".","baby"),
      r("p102720","56","^(pasta|spaghetti|macaroni|noodles).*whole","uncooked|dry"),
      r("p102740","562","^rice, brown","uncooked|dry|with.*(meat|chicken|vegetable|bean|cheese)"),
      r("p102780","56","^(barley|buckwheat|bulgur|quinoa|millet|wild rice).*cooked","with.*(meat|chicken|vegetable|bean|cheese)|uncooked"),
      r("p104320","7","^corn,","soup|salad|pudding|scalloped|baby|bread"),
      # Current WebQ .b combines mixed, wholemeal and seeded (codes 2-4).
      r("baguette.b","51","^bread,.*(french|vienna).*(wheat|whole grain|mixed grain|multigrain|seed)"),
      r("bap.b","51","^roll,.*(wheat|whole grain|mixed grain|multigrain|seed).*(hamburger|hot dog)"),
      r("roll.b","51","^roll,.*(wheat|whole grain|mixed grain|multigrain|brown|seed)"),
      r("bread.b","51","^bread,.*(wheat|whole grain|mixed grain|multigrain|multi-grain|brown|seed)")
    )
    map <- read.csv(file.path(NH_ROOT,"ukb/02foodwas/food_group_map.csv"),check.names=FALSE)
    ii <- match(rules$item,map$name)
    if (anyNA(ii)) stop("Current UKB mapping missing rule items.")
    rules$Food <- map$Food.Group[ii]; rules$portion_g <- as.numeric(map$Portion[ii])
    rules$UKB_item <- map$Food.Item[ii]
    rules<-rules[rules$Food %in% nh_scoring$dose_food_rules()$Food,]
    stopifnot(all(is.finite(rules$portion_g) & rules$portion_g>0),
      setequal(rules$Food,nh_scoring$dose_food_rules()$Food))
    rules
}
nh_map_codes <- function(dict,rules) {
    desc <- tolower(dict$food_desc); code <- sprintf("%08.0f",dict$food_code)
    # These codes record the ingredient's own grams, not whole-sandwich grams.
    desc <- sub(", for use on a sandwich$","",desc)
    rule_id <- rep(NA_integer_,nrow(dict))
    for (j in seq_len(nrow(rules))) {
        z <- grepl(paste0("^",rules$prefix[j]),code) & grepl(rules$include[j],desc,perl=TRUE)
        exclusion_text <- desc
        if(rules$item[j] %in% c("p100220","p100230","p100550","p100530","p100170","p100180")) {
            # Reconstituted/prepared drinks are liquid grams, even when made from powder.
            prepared <- grepl("reconstituted|made from dry mix|prepared with|made with",desc) & !grepl("unprepared",desc)
            exclusion_text[prepared] <- gsub("dry|powder|concentrate","",desc[prepared])
        }
        if(nzchar(rules$exclude[j])) z <- z & !grepl(rules$exclude[j],exclusion_text,perl=TRUE)
        # Additional conditions keep overlapping dairy/spread rules explicit.
        if(rules$item[j]=="p104430") z <- z & grepl("dried|dehydrated|^raisins?|^dates?",desc)
        if(rules$item[j] %in% c("p102850","p102810")) z <- z & grepl("low.?fat|reduced fat|fat.?free|nonfat|skim",desc)
        if(rules$item[j]=="p26111") z <- z & grepl("reduced|low.?fat|light|diet|fat.?free",desc)
        rule_id[is.na(rule_id) & z] <- j
    }
    dict$rule_id <- rule_id
    dict$item <- rules$item[rule_id]; dict$Food <- rules$Food[rule_id]
    dict$portion_g <- rules$portion_g[rule_id]
    dict
}
nh_dictionary <- function(cycle) {
    d <- nh_xpt(cycle,if(cycle<"2003") "DRXFMT" else "DRXFCD")
    if("START" %in% names(d)) {
        d <- d[grepl("IFDCD",d$FMTNAME),]
        ans <- data.frame(food_code=as.numeric(d$START),food_desc=d$LABEL)
    } else {
        long<-intersect(c("DRXFCLD","DRXFDLD"),names(d))
        if(!length(long))stop("No full food description for ",cycle)
        ans<-data.frame(food_code=d$DRXFDCD,food_desc=d[[long[1]]])
    }
    ans <- unique(ans); if(anyDuplicated(ans$food_code)) stop("Duplicate dictionary code: ",cycle)
    ans$cycle <- cycle; ans
}
nh_day <- function(cycle,day,codes,foods) {
    early <- cycle<"2003"; pre <- if(early) "DRX" else paste0("DR",day)
    tot <- nh_xpt(cycle,paste0(pre,"TOT")); iff <- nh_xpt(cycle,paste0(pre,"IFF"))
    status <- nh_col(tot,c(paste0(pre,"DRSTZ"),"DRDDRSTZ","DRDDRSTS"))
    weight_var <- if(cycle<"2003") "WTDR4YR" else "WTDRD1"
    d <- data.frame(SEQN=tot$SEQN,cycle,day,recall_status=status,
        energy_kcal=nh_col(tot,paste0(pre,"TKCAL")),alcohol_g=nh_col(tot,paste0(pre,"TALCO")),
        wt_day1=nh_col(tot,weight_var),wt_two=nh_col(tot,"WTDR2D",FALSE),weight_var)
    if(anyDuplicated(d$SEQN)) stop("Duplicate dietary totals.")
    x <- data.frame(SEQN=iff$SEQN,food_code=nh_col(iff,c(paste0(pre,"IFDCD"),"DRDIFDCD")),
                    grams=nh_col(iff,paste0(pre,"IGRMS")))
    k <- match(x$food_code,codes$food_code)
    x$Food <- codes$Food[k]; x$portion_g <- codes$portion_g[k]
    # Unknown food codes or missing/negative amounts do not become zero intake.
    bad_ids <- unique(x$SEQN[is.na(k)|!is.finite(x$grams)|x$grams<0])
    x$servings <- x$grams/x$portion_g
    for(food in foods) {
        a <- x[x$Food %in% food,]
        sm <- tapply(a$servings,a$SEQN,sum)
        d[[food]] <- unname(sm[as.character(d$SEQN)])
        d[[food]][is.na(d[[food]])] <- 0
    }
    d$mapping_complete <- !d$SEQN %in% bad_ids & d$SEQN %in% x$SEQN
    d$valid <- !is.na(status) & status==1 & is.finite(d$energy_kcal) & d$energy_kcal>0 & d$mapping_complete
    d[!d$valid,foods] <- NA_real_
    counts <- x %>% group_by(food_code) %>% summarise(records=n(),grams=sum(grams),.groups="drop")
    counts$cycle <- cycle; counts$day <- day
    list(data=d,counts=counts)
}
nh_covariates <- function(cycle) {
    d <- nh_xpt(cycle,"DEMO"); mcq <- nh_xpt(cycle,"MCQ"); smq <- nh_xpt(cycle,"SMQ")
    bmx <- nh_xpt(cycle,"BMX"); diq <- nh_xpt(cycle,"DIQ"); bpq <- nh_xpt(cycle,"BPQ")
    get <- function(x,cols,required=TRUE) nh_col(x,cols,required)[match(d$SEQN,x$SEQN)]
    ever <- get(smq,"SMQ020"); now <- get(smq,"SMQ040")
    edu <- nh_col(d,"DMDEDUC2"); eth <- nh_col(d,"RIDRETH1") # same 5 categories across cycles
    cv <- sapply(c("MCQ160B","MCQ160C","MCQ160D","MCQ160E"),function(v)nh_binary(get(mcq,v)))
    cvd <- ifelse(rowSums(cv==1,na.rm=TRUE)>0,1L,ifelse(rowSums(is.na(cv))==0,0L,NA_integer_))
    ans<-data.frame(SEQN=d$SEQN,cycle,age=nh_col(d,"RIDAGEYR"),
      sex=factor(nh_col(d,"RIAGENDR"),levels=1:2,labels=c("Male","Female")),
      ethnicity=factor(eth,levels=1:5),education=factor(ifelse(edu %in% 1:5,edu,NA)),
      income=nh_col(d,"INDFMPIR"),bmi=get(bmx,"BMXBMI"),
      smoking=factor(ifelse(ever==2,"Never",ifelse(ever==1 & now==3,"Former",ifelse(ever==1 & now %in% 1:2,"Current",NA))),levels=c("Never","Former","Current")),
      hypertension=factor(nh_binary(get(bpq,"BPQ020"))),
      diabetes_raw=get(diq,"DIQ010"),
      diabetes=factor(ifelse(get(diq,"DIQ010")==2,"No",ifelse(get(diq,"DIQ010")==3,"Borderline",ifelse(get(diq,"DIQ010")==1,"Yes",NA))),levels=c("No","Borderline","Yes")),cvd=factor(cvd),
      stroke_prevalent=nh_binary(get(mcq,"MCQ160F")),
      psu=nh_col(d,"SDMVPSU"),strata=nh_col(d,"SDMVSTRA"))
    left_join(ans,nh_aligned_covariates(cycle),by=c("SEQN","cycle"))
}
nh_attach_score <- function(d,foods) {
    input <- d[c("SEQN",foods)]; input$eid <- paste(d$cycle,d$SEQN,sep="_")
    attr(input,"intake_unit") <- "servings/day"
    score <- nh_scoring$calculate_dose(input)
    d <- cbind(d,score[setdiff(names(score),"eid")]); d$eid <- input$eid
    d$diet.dose.sum <- d$dose
    d$energy_kJ <- d$energy_kcal*4.184
    stopifnot(identical(d$dose,d$dose9half),all(is.na(d$dose)|d$dose %in% seq(0,9,.5)))
    attr(d,"primary_scoring_method") <- "dose9half"
    attr(d,"score_intake_unit") <- "servings/day"; d
}
run_nhanes_score <- function() {
    dir.create(NH_OUT,recursive=TRUE,showWarnings=FALSE)
    provenance<-nh_provenance()
    rules <- nh_food_rules(); foods <- unique(nh_scoring$dose_food_rules()$Food)
    nh_csv(rules,"food_mapping_rules.csv")
    dd <- cc <- cv <- dicts <- list()
    for(cycle in NH_CYCLES$cycle) {
        message("NHANES food mapping: ",cycle)
        codes <- nh_map_codes(nh_dictionary(cycle),rules); dicts[[cycle]] <- codes
        cv[[cycle]] <- nh_covariates(cycle)
        for(day in if(cycle<"2003") 1 else 1:2) {
            ans <- nh_day(cycle,day,codes,foods); key <- paste(cycle,day)
            dd[[key]] <- ans$data; cc[[key]] <- ans$counts
        }
    }
    days <- bind_rows(dd); cov <- bind_rows(cv)
    audit <- left_join(bind_rows(cc),bind_rows(dicts),by=c("cycle","food_code"))
    nh_csv(audit,"food_code_audit.csv")
    nh_csv(audit[!is.na(audit$Food),],"food_code_matches.csv")
    nh_csv(days %>% group_by(cycle,day) %>% summarise(N=n(),reliable=sum(recall_status==1,na.rm=TRUE),
        mapping_incomplete=sum(!mapping_complete),valid=sum(valid),.groups="drop"),"recall_QC.csv")
    # No selecting only positive-food days: absent foods contribute 0 on EVERY valid day.
    one <- days[days$day==1,]; one$weight <- one$wt_day1*ifelse(one$cycle<"2003",4/20,2/20)
    one <- left_join(one,cov,by=c("SEQN","cycle")); one$recall_mode <- "day1_1999_2018"
    one <- nh_attach_score(one,foods)
    both <- days %>% filter(cycle>="2003") %>% group_by(SEQN,cycle) %>%
      filter(n()==2,all(valid),all(is.finite(wt_two) & wt_two>0)) %>%
      summarise(across(all_of(c(foods,"energy_kcal","alcohol_g")),mean),wt_two=first(wt_two),.groups="drop")
    both$weight <- both$wt_two/8; both$valid <- TRUE
    both <- left_join(both,cov,by=c("SEQN","cycle")); both$recall_mode <- "two_day_2003_2018"
    both <- nh_attach_score(both,foods)
    if(!identical(provenance,nh_provenance()))stop("Sources changed during scoring; rerun before saving cache.")
    attr(one,"provenance")<-attr(both,"provenance")<-provenance
    dose_save_rds(one,file.path(NH_OUT,"nhanes_day1_score.rds"))
    dose_save_rds(both,file.path(NH_OUT,"nhanes_two_day_score.rds"))
    dose_save_rds(list(rules=rules,codes=audit,recall_qc=days %>% group_by(cycle,day) %>%
      summarise(N=n(),Valid=sum(valid),Incomplete=sum(!mapping_complete),.groups="drop"),
      covariates=cov,provenance=provenance),file.path(NH_OUT,"scoring_audit.rds"))
    invisible(list(day1=one,two_day=both))
}
# Source fingerprints and raw-file metadata are captured WHEN the cache is built.
nh_provenance<-function() {
    source<-c(file.path(dose_script_dir(),c("nhanes.R","score.R")),
      file.path(NH_ROOT,"ukb/02foodwas/food_group_map.csv"))
    raw<-sort(unique(c(list.files(NH_DATA,pattern="\\.(xpt|sas7bdat|txt)$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE),
      list.files(NH_FPED_DIR,pattern="(fped_.*|equiv.*)\\.(sas7bdat|txt)$",recursive=TRUE,full.names=TRUE,ignore.case=TRUE))))
    fi<-file.info(raw)
    list(source_md5=tools::md5sum(source),raw=data.frame(path=raw,size=fi$size,mtime=as.numeric(fi$mtime)),
      dietaryindex=as.character(packageVersion("dietaryindex")),
      wjfrt_source="dietaryindex::NHANES_20032004$WJFRT")
}
nh_read_mped<-function(cycle) {
    if(cycle=="2003-2004")return(as.data.frame(haven::read_sas(nh_usda_file("equiv0304.sas7bdat"))))
    tag<-if(cycle=="1999-2000")"9400" else "0102"
    p<-nh_usda_file(paste0("equiv",tag,".txt"))
    x<-read.fwf(p,widths=c(8,1,6,rep(8,32)),header=FALSE)
    names(x)<-c("DRDIFDCD","EQUIVFLAG","MODCODE","G_TOTAL","G_WHL","G_NWHL","V_TOTAL","V_DRKGR",
      "V_DPYEL","V_POTATO","V_STARCY","V_TOMATO","V_OTHER","F_TOTAL","F_CITMLB","F_OTHER",
      "D_TOTAL","D_MILK","D_YOGURT","D_CHEESE","M_MPF","M_MEAT","M_ORGAN","M_FRANK","M_POULT",
      "M_FISH_HI","M_FISH_LO","M_EGG","M_SOY","M_NUTSD","LEGUMES","DISCFAT_OIL","DISCFAT_SOL","ADD_SUG","A_BEV")
    # Nonzero modifiers in equiv9400 belong to CSFII 1994-1998, not NHANES.
    # USDA MPED documentation, Table 2: 1999-2000 uses the 7,360 base codes.
    x<-x[x$MODCODE==0,]
    stopifnot(!anyDuplicated(x$DRDIFDCD));x
}
# DASH uses total fruit only. Remove the unused whole-fruit/juice join, without
# altering the eight DASH inputs. Disable the 2003 pizza correction before 2003.
# Changes are confined to a local copy of the installed package function.
nh_mped_function<-function(early=FALSE) {
    f<-dietaryindex::DASH_NHANES_MPED
    walk<-function(x) {
        if(!is.call(x))return(x)
        tx<-paste(deparse(x),collapse=" ")
        if(identical(x[[1]],as.name("if")) && grepl('"ModCode" %in% colnames(WJFRT)',tx,fixed=TRUE))
          return(quote(MPED_PER_100_GRAM_4 <- MPED_PER_100_GRAM_3))
        if(as.character(x[[1]])[1] %in% c("<-","=")) {
          lhs<-as.character(x[[2]])[1]
          if(early && lhs=="MPED_PER_100_GRAM_3")return(quote(MPED_PER_100_GRAM_3 <- MPED_PER_100_GRAM_2))
          if(lhs=="selected_columns" && grepl('"WHOLEFRT"',tx,fixed=TRUE))return(quote(selected_columns <- selected_columns))
        }
        for(i in seq_along(x)[-1])x[[i]]<-walk(x[[i]])
        x
    }
    body(f)<-walk(body(f));f
}
nh_build_dash<-function() {
    provenance<-nh_provenance()
    ans<-list(); details<-list()
    for(cycle in NH_CYCLES$cycle) {
      message("DASH: ",cycle)
      suf<-NH_CYCLES$suffix[match(cycle,NH_CYCLES$cycle)]
      if(cycle<"2005") {
        mp<-nh_read_mped(cycle);f<-nh_mped_function(cycle<"2003")
        args<-list(MPED_PER_100_GRAM_PATH=mp,DEMO_PATH=nh_xpt(cycle,"DEMO"),
          NUTRIENT_PATH=nh_xpt(cycle,if(cycle<"2003")"DRXTOT" else "DR1TOT"),
          NUTRIENT_IND_PATH=nh_xpt(cycle,if(cycle<"2003")"DRXIFF" else "DR1IFF"))
        # Check every reliable recall food key before package inner joins.
        for(day in if(cycle<"2003")1 else 1:2) {
          iff<-if(day==1)args$NUTRIENT_IND_PATH else nh_xpt(cycle,"DR2IFF")
          status<-nh_col(iff,c(paste0("DR",day,"DRSTZ"),"DRDDRSTZ","DRDDRSTS"))
          code<-nh_col(iff,c(paste0("DR",day,"IFDCD"),"DRDIFDCD"))
          key<-if(cycle<"2003")code else paste(code,iff[[paste0("DR",day,"MC")]])
          mk<-if(cycle<"2003")mp$DRDIFDCD else paste(mp$DRDIFDCD,mp$DRDIMC)
          stopifnot(!anyDuplicated(mk),all(key[status %in% 1] %in% mk))
        }
        a<-do.call(f,args);b<-NULL
        if(cycle=="2003-2004") {
          args$NUTRIENT_PATH2<-nh_xpt(cycle,"DR2TOT");args$NUTRIENT_IND_PATH2<-nh_xpt(cycle,"DR2IFF")
          b<-do.call(f,args)
          # Verify removal of unused fruit fields leaves native results identical.
          args$WJFRT<-nh_wjfrt()
          native<-do.call(dietaryindex::DASH_NHANES_MPED,args)
          stopifnot(isTRUE(all.equal(b$DASH_ALL,native$DASH_ALL)),identical(b$SEQN,native$SEQN))
        }
      } else {
        tag<-paste0(substr(cycle,3,4),substr(cycle,8,9))
        fp<-function(day)nh_usda_file(paste0("fped_dr",day,"iff_",tag,".sas7bdat"))
        nut<-function(day)nh_xpt_path(cycle,paste0("DR",day,"IFF"))
        a<-dietaryindex::DASH_NHANES_FPED(FPED_IND_PATH=fp(1),NUTRIENT_IND_PATH=nut(1))
        b<-dietaryindex::DASH_NHANES_FPED(FPED_IND_PATH=fp(1),NUTRIENT_IND_PATH=nut(1),FPED_IND_PATH2=fp(2),NUTRIENT_IND_PATH2=nut(2))
      }
      z<-data.frame(SEQN=a$SEQN,dash_day1=a$DASH_ALL)
      if(!is.null(b))z<-full_join(z,data.frame(SEQN=b$SEQN,dash_two=b$DASH_ALL),by="SEQN") else z$dash_two<-NA_real_
      z$cycle<-cycle;stopifnot(!anyDuplicated(z$SEQN),all(z$dash_day1[is.finite(z$dash_day1)] %in% 8:40))
      ans[[cycle]]<-z;details[[cycle]]<-list(day1=a,two=b)
    }
    if(!identical(provenance,nh_provenance()))stop("Sources changed during DASH scoring; rerun before saving cache.")
    z<-bind_rows(ans);attr(z,"provenance")<-provenance
    dose_save_rds(z,file.path(NH_OUT,"nhanes_dash_score.rds"))
    dose_save_rds(details,file.path(NH_OUT,"dash_components.rds"));invisible(z)
}
if(isTRUE(getOption("dose.nhanes.run",TRUE))) {
    nh_check_usda()
    invisible(nh_wjfrt())
    expected<-nh_provenance()
    current<-function(fn) {p<-file.path(NH_OUT,fn);!dose_replace()&&file.exists(p)&&identical(attr(readRDS(p),"provenance"),expected)}
    if(isTRUE(getOption("dose.nhanes.rebuild",FALSE)) || !all(vapply(c("nhanes_day1_score.rds","nhanes_two_day_score.rds"),current,logical(1))))run_nhanes_score()
    if(isTRUE(getOption("dose.nhanes.dash.rebuild",FALSE)) || !current("nhanes_dash_score.rds"))nh_build_dash()
}
}))

# ===========================================================================
# nhanes / patterns
# ===========================================================================
dose_register_step("nhanes", "patterns", expression({
# NHANES additional patterns. Sourced automatically by nhanes.R.
# MIND uses current UKB rules with explicit NHANES food-item/portion adaptations.
# aMED uses dietaryindex MED (0-9), distinct from UKB's 13-component MedDiet.
suppressPackageStartupMessages({library(dplyr);library(haven)})
nh_pattern_provenance<-function()list(base=nh_provenance(),source=tools::md5sum(c(
 file.path(dose_script_dir(), 'nhanes.R'),file.path(dose_script_dir(), 'foodwas.R'),file.path(dose_files_dir(), "foods.xlsx"))))
nh_mind_spec<-function() {
 expr<-dose_step_code("foodwas", "preprocess")
 pick<-function(n) {
  z<-list()
  walk<-function(x) {
   if(is.expression(x)) {for(y in as.list(x))walk(y);return(invisible(NULL))}
   if(!is.call(x))return(invisible(NULL))
   if(identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(n))) {
    z[[length(z)+1L]]<<-x
   } else if(identical(x[[1]],as.name('{'))) {
    for(y in as.list(x)[-1L])walk(y)
   } else if(identical(x[[1]],as.name('dose_cached_phase'))) {
    # Inspect cached code without evaluating its inputs, outputs or other analyses.
    # Do not descend into unrelated calls or function-local definitions.
    walk(x[['code']])
   }
  }
  walk(expr)
  if(!length(z))stop('Missing UKB MIND definition: ',n)
  if(length(z)!=1L)stop('Ambiguous UKB MIND definition: ',n,' (found ',length(z),')')
  z[[1L]]
 }
 e<-new.env(parent=globalenv());eval(pick('mind_map'),e)
 list(map=e$mind_map,calculation=pick('dat.mind'))
}
nh_mind_rules<-function(spec) {
 r<-function(item,prefix,include,exclude='')data.frame(item,prefix,include,exclude)
 # New rules precede reused DOSE rules when food descriptions overlap.
 more<-bind_rows(
  r('sauce.pnutbutter','4','^(peanut|nut|almond|cashew).*butter','sandwich|cookie'),
  r('sauce.hummus','4','^(hummus|houmous)'),r('sauce.guacamole','7','^guacamole'),
  r('sauce.tomato','[27]','^(sauce,.*(tomato|spaghetti|marinara)|tomato sauce|spaghetti sauce)','meat|cheese'),
  r('canned.vegetable','[127]','^soup,.*(vegetable|tomato|broccoli|mushroom|onion|carrot|squash|celery)','dry|condensed|meat|chicken|beef|noodle|bean|lentil|pea'),
  r('canned.pulse','[1247]','^soup,.*(bean|lentil|split pea)','dry|condensed'),
  r('canned.fish','[127]','^(soup|chowder|bisque).*(fish|seafood|clam|lobster|shrimp)','dry|condensed'),
  r('p104080','7','^coleslaw'),r('p104090','7','^(salad,.*(lettuce|vegetable|garden)|lettuce salad|vegetable salad)','meat|chicken|egg|tuna|potato|pasta|macaroni'),
  r('p104470','6','^(strawberr|blueberr|bluberr|raspberr|blackberr|cranberr|boysenberr|mulberr|gooseberr|currant|berries)','juice|dried|sauce|jam|baby'),
  r('p104130','7','^beets?([, ]|$)','soup|salad'),r('p104170','7','^carrots?([, ]|$)','soup|cake|juice|baby'),
  r('p104190','7','^celery([, ]|$)','soup'),r('p104270','7','^parsnips?([, ]|$)'),
  r('p104360','7','^(turnip|rutabaga)([, ]|$)','greens'),r('p104220','7','^garlic([, ]|$)'),
  r('p104230','7','^leeks?([, ]|$)'),r('p104260','7','^onions?([, ]|$)','soup|rings|breaded|battered'),
  r('p104150','7','^squash,.*(winter|butternut|acorn|hubbard)'),r('p104290','7','^peppers?([, ]|$)','stuffed'),
  r('p104330','7','^sweet potato','pie'),r('p104060','7','^(vegetables?, mixed|mixed vegetables?)','soup|with.*meat'),
  r('p104100','6','^avocado'),r('p104200','7','^(zucchini|squash,.*summer|squash,.*zucchini)'),
  r('p104210','7','^cucumber','pickl'),r('p104180','7','^cauliflower','soup'),r('p104250','7','^mushrooms?([, ]|$)','soup|stuffed'),
  r('p104380','7','^(asparagus|artichoke|okra|eggplant|bamboo|pumpkin|vegetable pieces)([, ]|$)','pie|soup|casserole'),
  r('p103010','2','^(bologna|salami|pepperoni|luncheon meat|deli meat|pastrami)','sandwich'),
  r('p103090','2','^(liver|beef,.*liver|pork,.*liver|chicken,.*liver)','sandwich'),
  r('p103020','2','^beef([, ]|$)','sandwich|soup|stew|with.*(rice|noodle|pasta)|liver|sausage|frankfurter'),
  r('p103030','2','^pork([, ]|$)','bacon|ham|sausage|sandwich|soup|stew|with.*(rice|noodle|pasta)|liver'),
  r('p103040','2','^lamb([, ]|$)','sandwich|soup|stew|with.*(rice|noodle|pasta)|liver'),
  r('p103050','2','^(chicken|turkey|poultry).*(fried|breaded|battered|nugget|patty)','stir.fried|sandwich|soup'),
  r('p103060','2','^(chicken|turkey|duck|goose|poultry)([, ]|$)','sandwich|soup|stew|with.*(rice|noodle|pasta)|liver|ham|bacon|sausage|fried|breaded|battered|nugget'),
  r('p103180','2','^(fish|cod|haddock|pollock|whiting|catfish|flounder|tilapia).*(batter|fried)','stir.fried|sandwich|soup'),
  r('p103170','2','^(fish|cod|haddock|pollock|whiting|catfish|flounder|tilapia).*(breaded|stick|nugget)','sandwich|soup'),
  r('p103150','2','^tuna.*canned','salad|sandwich'),
  r('p103160','2','^(salmon|sardine|mackerel|herring|trout|anchov)([, ]|$)','fried|breaded|batter|sandwich|salad'),
  r('p103190','2','^(cod|haddock|pollock|whiting|catfish|flounder|tilapia)([, ]|$)','fried|breaded|batter|sandwich|salad'),
  r('p103200','2','^(shrimp|prawn)([, ]|$)','fried|breaded|batter|sandwich|salad'),
  r('p103210','2','^(lobster|crab)([, ]|$)','fried|breaded|batter|cake|salad'),
  r('p103220','2','^(clam|mussel|oyster|scallop|squid|octopus)([, ]|$)','fried|breaded|batter|salad'),
  r('p103230','2','^(fish|tuna|sea bass|halibut|swordfish|snapper|perch)([, ]|$)','fried|breaded|batter|sandwich|salad|soup'),
  r('p26063','8','^butter([, ]|$)'),r('p26111','8','^margarine.*(reduced|low.fat|light|diet)'),r('p26112','8','^margarine'),
  r('p26154','1','^cream([, ]|$)','sour|substitute|non.dairy|imitation'),
  r('p102860','1','^cheese,.*(cream|spread)','low.fat|reduced fat|fat.free|nonfat'),
  r('p102840','1','^cheese,.*blue','low.fat|reduced fat|fat.free|nonfat'),
  r('p102880','1','^cheese,.*feta','low.fat|reduced fat|fat.free|nonfat'),
  r('p102890','1','^cheese,.*mozzarella','low.fat|reduced fat|fat.free|nonfat|skim'),
  r('p102900','1','^cheese,.*goat','low.fat|reduced fat|fat.free|nonfat'),
  r('p102830','1','^cheese,.*(brie|camembert|soft)','low.fat|reduced fat|fat.free|nonfat'),
  r('p102820','1','^cheese,','low.fat|reduced fat|fat.free|nonfat|cottage|skim|sauce'),
  r('p102450','4','^(seeds|sunflower|pumpkin seed|sesame|flax|chia)([, ]|$)','butter|bread'),
  r('p102420','4','^peanuts?([, ]|$)','butter|candy|chocolate'),
  r('p102440','4','^(almonds?|walnuts?|pecans?|cashews?|pistachios?|hazelnuts?|macadamia|nuts|mixed nuts|brazil nuts)([, ]|$)','butter|candy|chocolate'),
  r('p104020','7','^potato.*(french.fr|fried|hash brown|tater tot)','chips|with.*(chili|cheese)'),r('p102460','[57]','^(potato chips|chips|corn chips|tortilla chips)'),
  r('p102350','5','^cookies?.*chocolate'),r('p102360','5','^cookies?'),r('p102370','5','^(cereal bar|granola bar)'),
  r('p102220','5','^cheesecake'),r('p102180','5','^cake,.*fruit'),r('p102190','5','^cake([, ]|$)'),
  r('p102200','5','^(doughnut|donut)'),r('p101970','5','^pie([, ]|$)'),r('p102010','5','^(pancake|waffle)'),
  r('p102050','5','^croissant'),r('p102060','5','^(danish|pastry)'),r('p102070','5','^scone'),
  r('p102120','1','^(ice cream|ice milk|frozen yogurt|sherbet)'),r('p102140','1','^pudding'),
  r('p102300','[69]','^(raisins|nuts).*chocolate'),r('p102270','917','white chocolate'),
  r('p102290','917','dark.*chocolate|chocolate.*dark'),r('p102280','917','chocolate|choc\\.'),
  r('p102320','917','sugar.free'),r('p102330','917','candy|candies|fudge|nougat|caramel|toffee|marshmallow|reese|twix|snicker|toblerone|m&m'),
  r('p102380','917','.'),r('p102210','5','^pudding,.*sponge'),r('p102030','5','^yorkshire pudding'),
  r('p102170','1','^(soy|soya).*dessert'),r('p101990','5','^(crumble|cobbler)'),
  r('p103170','999','^fish fillet.*fried.*ingredient'),
  r('p101250','5','^(crispbread|crackers.*(rye|whole wheat))'),r('p101260','5','^oatcake'),r('p102770','56','^couscous.*cooked')
 )
 base<-nh_food_rules()[c('item','prefix','include','exclude')]
 rules<-bind_rows(more,base);target<-unique(c(unlist(spec$map),'p104090'));rules<-rules[rules$item %in% target,]
 # MIND distinguishes low-calorie drinks from the sugary-drink items in UKB.
 ix<-rules$item %in% c('p100170','p100180','p100550','p100530')
 rules$exclude[ix]<-paste0(rules$exclude[ix],'|diet|low.calorie|sugar.free|sugar free|reduced sugar')
 map<-as.data.frame(readxl::read_excel(file.path(dose_files_dir(), "foods.xlsx"),sheet='Sheet2'))
 stopifnot(!anyDuplicated(map$name));k<-match(rules$item,map$name)
 rules$portion_g<-as.numeric(map$Portion[k]);rules$Food<-rules$item
 stopifnot(!anyNA(k),all(is.finite(rules$portion_g)&rules$portion_g>0))
 rules
}
nh_mind_calculate<-function(d,spec) {
 e<-new.env(parent=globalenv());e$dat1<-data.frame(eid=d$eid)
 for(g in names(spec$map))e$dat1[[g]]<-rowSums(d[spec$map[[g]]])
 e$dat1$Green_leafy_vegetables<-e$dat1$Green_leafy_vegetables+.5*d$p104090
 e$dat1$Other_vegetables<-e$dat1$Other_vegetables+.5*d$p104090
 weekly<-setdiff(names(spec$map),c('Olive_oil','Butter_and_margarine','Whole_grain'))
 e$dat1[weekly]<-lapply(e$dat1[weekly],function(x)x*7)
 eval(spec$calculation,e)
 ans<-as.data.frame(e$dat.mind);names(ans)[names(ans)=='diet.mind.sum']<-'mind'
 stopifnot(all(is.na(ans$mind)|ans$mind %in% seq(0,15,.5)))
 list(score=ans,intakes=d)
}
nh_amed_mped_function<-function(early) {
 f<-dietaryindex::MED_NHANES_MPED
 # Same unused fruit-field removal and cycle-specific pizza safeguard as DASH.
 walk<-function(x) {
  if(!is.call(x))return(x);tx<-paste(deparse(x),collapse=' ')
  if(identical(x[[1]],as.name('if'))&&grepl('"ModCode" %in% colnames(WJFRT)',tx,fixed=TRUE))return(quote(MPED_PER_100_GRAM_4 <- MPED_PER_100_GRAM_3))
  if(as.character(x[[1]])[1] %in% c('<-','=')) {
   lhs<-as.character(x[[2]])[1]
   if(early&&lhs=='MPED_PER_100_GRAM_3')return(quote(MPED_PER_100_GRAM_3 <- MPED_PER_100_GRAM_2))
   if(lhs=='selected_columns'&&grepl('"WHOLEFRT"',tx,fixed=TRUE))return(quote(selected_columns <- selected_columns))
  }
  for(i in seq_along(x)[-1])x[[i]]<-walk(x[[i]]);x
 }
 body(f)<-walk(body(f));f
}
# Same validity definition as nh_day(), with a single grouped aggregation.
nh_pattern_day<-function(cycle,day,codes,items) {
 early<-cycle<"2003";pre<-if(early)"DRX" else paste0("DR",day)
 tot<-nh_xpt(cycle,paste0(pre,"TOT"));iff<-nh_xpt(cycle,paste0(pre,"IFF"))
 st<-nh_col(tot,c(paste0(pre,"DRSTZ"),"DRDDRSTZ","DRDDRSTS"))
 d<-data.frame(SEQN=tot$SEQN,cycle,day,wt_two=nh_col(tot,"WTDR2D",FALSE))
 energy<-nh_col(tot,paste0(pre,"TKCAL"));fc<-nh_col(iff,c(paste0(pre,"IFDCD"),"DRDIFDCD"));g<-nh_col(iff,paste0(pre,"IGRMS"))
 k<-match(fc,codes$food_code);bad<-unique(iff$SEQN[is.na(k)|!is.finite(g)|g<0])
 d$valid<-st %in% 1 & is.finite(energy)&energy>0 & !d$SEQN %in% bad & d$SEQN %in% iff$SEQN
 x<-data.frame(SEQN=iff$SEQN,item=codes$item[k],serv=g/codes$portion_g[k]);x<-x[!is.na(x$item),]
 # rowsum is keyed on SEQN+item; absent decoded items are zero within valid recalls.
 sm<-rowsum(x$serv,paste(x$SEQN,x$item,sep="|"),reorder=FALSE)
 for(nm in items) {v<-sm[match(paste(d$SEQN,nm,sep="|"),rownames(sm)),1];v[is.na(v)]<-0;v[!d$valid]<-NA_real_;d[[nm]]<-v}
 d
}
nh_build_patterns<-function() {
 prov<-nh_pattern_provenance();spec<-nh_mind_spec();rules<-nh_mind_rules(spec)
 items<-unique(c(unlist(spec$map),'p104090'));daily<-codeaudit<-amed<-amed_details<-list()
 for(cycle in NH_CYCLES$cycle) {
  message('MIND / aMED: ',cycle)
  codes<-nh_map_codes(nh_dictionary(cycle),rules);codeaudit[[cycle]]<-codes
  for(day in if(cycle<'2003')1 else 1:2) {
   dd<-nh_pattern_day(cycle,day,codes,items)
   # UKB sauces are presence indicators per recall, not grams/portion servings.
   for(nm in grep('^sauce\\.',names(dd),value=TRUE))dd[[nm]]<-as.numeric(dd[[nm]]>0)
   dd$eid<-paste(cycle,dd$SEQN,sep='_');daily[[paste(cycle,day)]]<-dd
  }
  suf<-NH_CYCLES$suffix[match(cycle,NH_CYCLES$cycle)]
  if(cycle<'2005') {
   mp<-nh_read_mped(cycle);f<-nh_amed_mped_function(cycle<'2003')
   args<-list(MPED_PER_100_GRAM_PATH=mp,DEMO_PATH=nh_xpt(cycle,'DEMO'),NUTRIENT_PATH=nh_xpt(cycle,if(cycle<'2003')'DRXTOT' else 'DR1TOT'),NUTRIENT_IND_PATH=nh_xpt(cycle,if(cycle<'2003')'DRXIFF' else 'DR1IFF'))
   for(day in if(cycle<'2003')1 else 1:2) {
    iff<-if(day==1)args$NUTRIENT_IND_PATH else nh_xpt(cycle,'DR2IFF')
    st<-nh_col(iff,c(paste0('DR',day,'DRSTZ'),'DRDDRSTZ','DRDDRSTS'))
    fc<-nh_col(iff,c(paste0('DR',day,'IFDCD'),'DRDIFDCD'))
    key<-if(cycle<'2003')fc else paste(fc,iff[[paste0('DR',day,'MC')]])
    mk<-if(cycle<'2003')mp$DRDIFDCD else paste(mp$DRDIFDCD,mp$DRDIMC)
    stopifnot(!anyDuplicated(mk),all(key[st %in% 1] %in% mk))
   }
   a<-do.call(f,args);b<-NULL
   if(cycle=='2003-2004') {
    args$NUTRIENT_PATH2<-nh_xpt(cycle,'DR2TOT');args$NUTRIENT_IND_PATH2<-nh_xpt(cycle,'DR2IFF');b<-do.call(f,args)
    args$WJFRT<-nh_wjfrt()
    native<-do.call(dietaryindex::MED_NHANES_MPED,args)
    stopifnot(identical(b$SEQN,native$SEQN),isTRUE(all.equal(b$MED_ALL,native$MED_ALL)))
   }
  } else {
   tag<-paste0(substr(cycle,3,4),substr(cycle,8,9))
   fp<-function(day)nh_usda_file(paste0('fped_dr',day,'tot_',tag,'.sas7bdat'))
   nt<-function(day)nh_xpt_path(cycle,paste0('DR',day,'TOT'))
   demo<-nh_xpt_path(cycle,'DEMO')
   a<-dietaryindex::MED_NHANES_FPED(FPED_PATH=fp(1),NUTRIENT_PATH=nt(1),DEMO_PATH=demo)
   b<-dietaryindex::MED_NHANES_FPED(FPED_PATH=fp(1),NUTRIENT_PATH=nt(1),DEMO_PATH=demo,FPED_PATH2=fp(2),NUTRIENT_PATH2=nt(2))
  }
  am<-data.frame(SEQN=a$SEQN,amed_day1=a$MED_ALL)
  if(!is.null(b))am<-full_join(am,data.frame(SEQN=b$SEQN,amed_two=b$MED_ALL),by='SEQN')else am$amed_two<-NA_real_
  am$cycle<-cycle;amed[[cycle]]<-am;amed_details[[cycle]]<-list(day1=a,two=b)
 }
 days<-bind_rows(daily);one<-days[days$day==1,]
 two<-one[one$cycle>='2003',];second<-days[days$day==2,];k<-match(two$eid,second$eid)
 good<-two$valid & second$valid[k] & is.finite(two$wt_two) & two$wt_two>0 & is.finite(second$wt_two[k]) & second$wt_two[k]>0
 for(nm in items)two[[nm]]<-(two[[nm]]+second[[nm]][k])/2
 two<-two[good %in% TRUE,]
 d1<-nh_mind_calculate(one,spec);d2<-nh_mind_calculate(two,spec)
 scores<-full_join(data.frame(eid=d1$score$eid,mind_day1=d1$score$mind),data.frame(eid=d2$score$eid,mind_two=d2$score$mind),by='eid')
 am<-bind_rows(amed);am$eid<-paste(am$cycle,am$SEQN,sep='_');scores<-full_join(scores,am,by='eid')
 stopifnot(!anyDuplicated(scores$eid),all(scores$amed_day1[is.finite(scores$amed_day1)] %in% 0:9))
 if(!identical(prov,nh_pattern_provenance()))stop('Pattern sources changed while running; rerun.')
 bundle<-list(scores=scores,mind=list(day1=d1,two=d2,spec=spec,rules=rules,codes=bind_rows(codeaudit)),amed=amed_details,provenance=prov,
   labels=c(mind='MIND-adapted',amed='aMED'),range=list(mind=c(0,15),amed=c(0,9)))
 dose_save_rds(bundle,file.path(NH_OUT,'nhanes_patterns.rds'));invisible(bundle)
}
nh_patterns<-function() {
 p<-file.path(NH_OUT,'nhanes_patterns.rds');prov<-nh_pattern_provenance()
 if(!dose_replace() && file.exists(p)) {a<-readRDS(p);if(identical(a$provenance,prov))return(a)}
 nh_build_patterns()
}
}))

# ===========================================================================
# nhanes / stroke
# ===========================================================================
dose_register_step("nhanes", "stroke", expression({
# Current DOSE9half / DASH, fourteen pooled-cycle comparisons; no omnibus test.
suppressPackageStartupMessages({library(survey);library(dplyr)})
ROOT<-dose_output_dir()
INPUT<-file.path(ROOT,".cache/nhanes")
OUT<-file.path(ROOT,"nhanes")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
UKB_CFG<-readRDS(file.path(ROOT,"ukb/00data/analysis_config.rds"))
COV_MAP<-c(age="age",sex.b="sex",ethnic.c3="ethnicity_ukb",edu.c="education_ukb",
  tdi="income",energy.kJ="energy_kJ",smoke.c="smoking",days_pa_mod.c="pa_moderate",
  bmi="bmi",alcohol.c="alcohol_status",fhist_stroke=NA_character_,
  primary_ahtn="hypertension",t2dm="diabetes",drug.lipid="drug_lipid",drug.htn="drug_htn",drug.dm="drug_dm")
map_cov<-function(x) {
    if(length(setdiff(x,names(COV_MAP))))stop("Unmapped UKB covariates: ",paste(setdiff(x,names(COV_MAP)),collapse=", "))
    unname(na.omit(COV_MAP[x]))
}
M2_BASE<-unique(c(map_cov(setdiff(UKB_CFG$covs.m1,UKB_CFG$covs.alcohol)),
  map_cov(UKB_CFG$covs.m2.add),"factor(cycle)"))
SCORES<-c("dose","dash","mind","amed")
SCORE_LABELS<-c(dose="DOSE9half",dash="DASH",mind="MIND-adapted",amed="aMED")
M2<-list(dose=M2_BASE,dash=unique(c(M2_BASE,map_cov(UKB_CFG$covs.alcohol))),mind=M2_BASE,amed=M2_BASE)
# Inverse weighted empirical CDF; ties stay intact.
weighted_quantile <- function(x,w,p) {
    ok <- is.finite(x)&is.finite(w)&w>0; x<-x[ok]; w<-w[ok]
    o<-order(x); x<-x[o]; cw<-cumsum(w[o])/sum(w)
    vapply(p,function(z)x[which(cw>=z)[1]],numeric(1))
}
extract_or <- function(fit,term,increment=1) {
    s<-summary(fit); b<-coef(fit)[term]; se<-sqrt(vcov(fit)[term,term])
    df<-fit$df.residual
    if(!isTRUE(fit$converged)||!is.finite(b)||!is.finite(se)||df<=0) stop("Invalid fit: ",term)
    crit<-qt(.975,df)
    data.frame(OR=exp(b*increment),LCI=exp((b-crit*se)*increment),
       UCI=exp((b+crit*se)*increment),p=s$coefficients[term,4],df=df)
}


write_nhanes_readme<-function(bundle) {
 r<-bundle$results
 fmt<-function(x) sprintf("%.3f (%.3f–%.3f)",x$OR,x$LCI,x$UCI)
 lines<-c("# NHANES：DOSE、DASH、MIND-adapted、aMED 周期比较", "",
 "## 定义与改动", "",
 "- DOSE 直接调用当前 f/score.R 的 DOSE9half（0–9 分，步长 0.5），份量来自当前 UKB food_group_map.csv。",
 "- 新增橄榄油、全谷物、葡萄酒；移除乳制品、人造黄油；苹果/梨只要求摄入 >0。水果及蔬菜分别为两个成分的均分。",
 "- 本次从原始回顾重算；修正三明治单列菠菜编码漏计；DIQ010=3 保留为 Borderline（No/Borderline/Yes），7/9/缺失仍排除。旧人数不作为目标。",
 "- 1999–2002 无第二日回顾，1999 起始窗口仅用首日。2003/2005/2007 起始窗口均比较首日与两日，分别截至 2016/2018，共 14 个窗口；均完整报告。",
 "- 同一窗口四种评分使用完全相同人群、权重、插补参数；均报告连续每分/每 SD、三组比较和 P for trend，不报告总体检验。",
 "- 结局 MCQ160F：自报既往任何卒中，非新发缺血性卒中；复杂抽样 Logistic 回归报告 OR。", "",
 "## M2 汇总", "", "|摄入|周期|N / 卒中|评分|每 SD OR (95% CI)|P 连续|T3 vs T1 OR (95% CI)|P trend|", "|---|---|---:|---|---|---:|---|---:|")
 for(set in unique(r$Analysis))for(sc in SCORES) {
   a<-subset(r,Analysis==set & Score==sc & Contrast=="Per SD")
   b<-subset(r,Analysis==set & Score==sc & Contrast=="T3")
   lines<-c(lines,sprintf("|%s|%s–%s|%s / %s|%s|%s|%.4g|%s|%.4g|",a$Recall,a$Start,a$End,a$N,a$Cases,SCORE_LABELS[sc],fmt(a),a$p,fmt(b),b$P_trend))
 }
 lines<-c(lines,"","## 调整与抽样","",
 paste0("DOSE M2：",paste(M2$dose,collapse=" + "),"。"),
 paste0("DASH M2：",paste(M2$dash,collapse=" + "),"。"),
 "新版 DOSE 含葡萄酒，按照当前 UKB model_covs() 逻辑不额外调整饮酒；DASH 仍调整饮酒。因此两者人群相同，协变量存在这一个有意的差别。",
 "从全部有正权重的原始回顾记录建立分层/PSU设计后再作成年人分析域筛选，避免先删除不合格回顾改变抽样设计。周期嵌套分层/PSU，加入周期固定效应；孤立 PSU 使用 adjust。",
 "首日：1999–2002 用 WTDR4YR × 4/合并年数，其余用 WTDRD1 × 2/合并年数。两日：WTDR2D/合并周期数。",
 "无年龄上限。民族/教育/吸烟/饮酒众数、BMI/收入均值插补，参照该窗口共同合格人群；体力活动 Unknown 单独保留。高血压或糖尿病状态未知排除。",
 "三组按加权经验分位点划分，同分者不拆分；P trend 使用组内加权中位数；每 SD 为该窗口加权标准差。CI 与 P 使用模型残余自由度的 t 分布。",
 "收入贫困比近似 TDI，中等强度休闲活动近似 UKB 活动天数；民族/教育分类也为近似；没有跨周期一致的卒中特异家族史，因此未纳入。", "",
 "## 新增 MIND 与 aMED 的定义及限制", "",
 "- MIND-adapted：15 个组分、0–15 分。直接提取 f/foodwas.R 的 mind_map 与 dat.mind 评分表达式；份量读取 D:/files/diet/foods.xlsx 的 Sheet2。两日先平均原始摄入，再应用 UKB 固定阈值。",
 "- MIND 食物逐项映射，独立于 DOSE 的组名：例如西兰花在 MIND 属其他蔬菜；MIND 全谷物不直接复用 DOSE Whole_grains。酱料按每天是否摄入计数（与 UKB 的二值输入一致）。",
 "- MIND-adapted 属操作性近似：菜肴中无法识别的配料不拆算；汤按明确命名的鱼/豆/蔬菜汤近似，非罐装汤可使用相同220g份量；油仅计明确报告橄榄油，不能证明其为家庭主要用油。不能称为原始 FFQ MIND 的完全复刻。",
 "- 部分美国食物采用最近似 UKB 项目及份量，如普通硬奶酪40g、普通蛋糕60g、华夫饼归入煎饼110g；食物映射和全部15组摄入/分数保存在 .cache/nhanes/nhanes_patterns.rds。",
 "- aMED：dietaryindex::MED_NHANES_MPED / MED_NHANES_FPED，0–9 分。包含水果、蔬菜、全谷物、豆类、坚果、鱼、红/加工肉、单不饱和/饱和脂肪比、酒精。中位数在每周期、每日、年龄≥2岁的有效回顾参考人群中确定；沿用包的零摄入和中位数并列值处理。两日分数取均值。",
 "- aMED 沿用包实现的酒精10–25 g/day计1分；非性别特异阈值，不应与其他aMED实现混称。它也不同于本项目UKB的13组分MedDiet（medi24）。[方法与实现说明](https://jamesjiadazhan.github.io/dietaryindex_manual/reference/MED_NHANES_FPED.html)。",
 "- MIND 和 aMED 均含酒精成分，M2与DOSE一样不额外调整饮酒；DASH仍调整饮酒。分析比较使用相同人群不等于各评分的组分和测量误差相同。", "", "## 原有 DOSE / DASH", "",
 "DOSE 将明确命名的食物克数除以对应 UKB 份量；无法拆解的混合菜肴成分不推算。橄榄油仅明确报告的油，不能视作所有烹调用油。玫瑰红葡萄酒按普通葡萄酒份量，米酒/酒精饮料混合饮品不计葡萄酒。",
 "Whole_grains 严格沿用当前 UKB 食物组的操作定义：其中包括早餐谷物、甜玉米，以及 mixed/wholemeal/seeded 面包（[UKB coding 81](https://biobank.ndph.ox.ac.uk/ukb/coding.cgi?id=81)），并非 USDA 的纯全谷物当量。DOSE 两日方案先平均摄入再打分。",
 "DASH 使用 dietaryindex 逐日原生实现，两个有效日的得分取平均。1999–2002 用 MPED1，2003–2004 用 MPED2，2005–2018 用 FPED；沿用包默认食品代码、按周期/日分位评分参考人群。",
 "1999–2000 只使用 equiv9400 中 MODCODE=0 的基础食品：非零修订码仅适用于 CSFII 1994–1998（[USDA MPED 文档表2](https://www.ars.usda.gov/ARSUserFiles/80400530/pdf/mped/mped2_doc.pdf)）。早期 MPED 包装函数去除不参与 DASH 的 whole-fruit/juice 拆分；2003 原生函数对照结果一致。2003 披萨修正只用于 2003–2004，早期不套用。可靠回顾的食物编码在计算前核对完整匹配。",
 "MPED/FPED 存在历史分类差异，且包默认的后期饮品/低脂乳品代码列表用于早期回顾可能产生测量差异；包含早期数据的 DASH 属敏感性分析。",
 "旧结果变动同时涉及评分、菠菜、Borderline 与饮酒调整，不能归因于单一因素。多周期探索未校正多重比较；不按最小 P 值选择主分析。既往卒中横断面关联存在反向因果。", "",
 "## 文件与运行", "",
 "图形仅导出 PNG：Fig_NHANES_M2.png（每 SD 周期比较）、Fig_NHANES_M2_tertiles.png（2003–2018 两日均值的三组及连续评分比较）、Sfig1_nhanes.png（两种最长窗口筛选图）、Sfig1_nhanes_primary.png（主分析人群简版筛选图）。数据与说明：NHANES_M2_results.csv、NHANES_analysis.rds、README.md。",
 "评分与映射中间数据位于 .cache/nhanes。每次读取检查评分源代码、UKB规则/份量表、原始文件大小/修改时间及包版本；失效时第一步自动重建。原始文件若被等大小且保留时间戳替换，需显式 rebuild。新增模式缓存还检查辅助脚本、UKB MIND源码和foods.xlsx；第二步自动补算失效的新增模式缓存。",
 "```bash", "./dose.sh nhanes", "```",
 "新增模式前备份：_archive/nhanes_before_patterns_20260919。辅助评分脚本：f/nhanes.R，由第二步自动调用。")
 writeLines(enc2utf8(lines),file.path(OUT,"README.md"),useBytes=TRUE)
}
run_nhanes_stroke<-function() {
 old<-options(survey.lonely.psu="adjust",survey.adjust.domain.lonely=TRUE);on.exit(options(old))
 env<-new.env(parent=globalenv());op<-options(dose.nhanes.run=FALSE)
 dose_source_step("nhanes", "score", envir=env);options(op)
 provenance<-env$nh_provenance()
 read_checked<-function(fn) {
   d<-readRDS(file.path(INPUT,fn))
   if(!identical(attr(d,"provenance"),provenance))stop("Stale cache; rerun nhanes.R: ",fn)
   d
 }
 dose_source_step("nhanes", "patterns", envir=env)
 patterns<-env$nh_patterns();ps<-patterns$scores
 one<-read_checked("nhanes_day1_score.rds");two<-read_checked("nhanes_two_day_score.rds");dash<-read_checked("nhanes_dash_score.rds")
 cov<-readRDS(file.path(INPUT,"scoring_audit.rds"))$covariates
 stopifnot(!anyDuplicated(one$eid),!anyDuplicated(two$eid),!anyDuplicated(dash[c("SEQN","cycle")]))
 windows<-rbind(expand.grid(Start=c(1999,2003,2005,2007),End=c(2016,2018),Recall="Day 1"),
   expand.grid(Start=c(2003,2005,2007),End=c(2016,2018),Recall="Two-day mean"))
 rows<-cohorts<-audit<-list()
 for(j in seq_len(nrow(windows))) {
   q<-windows[j,];set<-paste(if(q$Recall=="Day 1")"day1" else "two_day",q$Start,q$End,sep="_")
   message("M2: ",set)
   d<-as.data.frame(one); d$year<-as.integer(substr(d$cycle,1,4));d<-d[d$year>=q$Start & d$year<q$End,]
   if(q$Recall=="Two-day mean") {
     k<-match(d$eid,two$eid)
     score_fields<-setdiff(names(env$nh_scoring$calculate_dose(
       structure(two[1,c("eid",env$nh_scoring$dose_food_rules()$Food)],intake_unit="servings/day"))),"eid")
     replace<-c(env$nh_scoring$dose_food_rules()$Food,score_fields,"diet.dose.sum","energy_kcal","energy_kJ","alcohol_g","valid","recall_mode")
     for(nm in replace)d[[nm]]<-two[[nm]][k]
     d$weight<-d$wt_two/((q$End-q$Start+1)/2)
   } else d$weight<-d$wt_day1*ifelse(d$year<2003,4,2)/(q$End-q$Start+1)
   dk<-match(paste(d$cycle,d$SEQN),paste(dash$cycle,dash$SEQN))
   d$dash<-dash[[if(q$Recall=="Day 1")"dash_day1" else "dash_two"]][dk]
   pk<-match(d$eid,ps$eid)
   for(sc in c("mind","amed"))d[[sc]]<-ps[[paste0(sc,if(q$Recall=="Day 1")"_day1" else "_two")]][pk]
   design_ok<-with(d,is.finite(weight)&weight>0&!is.na(psu)&!is.na(strata))
   adult<-!is.na(d$age)&d$age>=20
   dietary<-adult & d$valid %in% TRUE & is.finite(d$dose) & design_ok
   known<-dietary & !is.na(d$stroke_prevalent)
   disease<-known & complete.cases(d[c("hypertension","diabetes")])
   reference<-disease & complete.cases(d[c("age","sex","energy_kJ",SCORES)])
   imp<-list()
   for(nm in c("ethnicity_ukb","education_ukb","smoking","alcohol_status","bmi","income")) {
     v<-d[[nm]][reference];value<-if(is.numeric(v))mean(v,na.rm=TRUE) else names(which.max(table(v)))
     stopifnot(length(value)==1,!is.na(value));imp[[nm]]<-list(value=value,n=sum(reference&is.na(d[[nm]])))
     d[[nm]][is.na(d[[nm]])]<-value
   }
   required<-c(setdiff(unique(unlist(M2)),"factor(cycle)"),"stroke_prevalent",SCORES)
   d$complete<-reference & complete.cases(d[required])
   d$stratum_id<-interaction(d$cycle,d$strata,drop=TRUE);d$psu_id<-interaction(d$cycle,d$strata,d$psu,drop=TRUE)
   full<-svydesign(ids=~psu_id,strata=~stratum_id,weights=~weight,nest=TRUE,data=d[design_ok,])
   des<-subset(full,complete);z<-des$variables
   period<-cov[as.integer(substr(cov$cycle,1,4))>=q$Start & as.integer(substr(cov$cycle,1,4))<q$End,]
   flow<-c(Period=nrow(period),Adults=sum(period$age>=20,na.rm=TRUE),Reliable_diet=sum(dietary),Known_stroke=sum(known),Known_diseases=sum(disease),Common_M2=nrow(z))
   stopifnot(all(diff(flow)<=0),nrow(z)>0,!anyNA(z[required]),all(z$dose==z$dose9half))
   audit[[set]]<-list(window=q,imputation=imp,flow=flow,excluded=-diff(flow),cases=sum(z$stroke_prevalent),
     borderline=sum(z$diabetes=="Borderline"),missing_scores=sapply(d[SCORES],function(x)sum(disease & !is.finite(x))),full_design_N=sum(design_ok),design_df=degf(des))
   for(score in SCORES) {
     ds<-des;x<-z[[score]];w<-z$weight;br<-c(min(x),weighted_quantile(x,w,c(1/3,2/3)),max(x))
     stopifnot(!anyDuplicated(br));ds$variables$x<-x
     ds$variables$group<-cut(x,br,include.lowest=TRUE,right=TRUE,labels=c("T1","T2","T3"))
     med<-vapply(levels(ds$variables$group),function(g){i<-ds$variables$group==g;weighted_quantile(x[i],w[i],.5)},numeric(1))
     ds$variables$trend<-unname(med[as.character(ds$variables$group)])
     fc<-svyglm(reformulate(c("x",M2[[score]]),"stroke_prevalent"),ds,family=quasibinomial())
     fg<-svyglm(reformulate(c("group",M2[[score]]),"stroke_prevalent"),ds,family=quasibinomial())
     ft<-svyglm(reformulate(c("trend",M2[[score]]),"stroke_prevalent"),ds,family=quasibinomial())
     pt<-extract_or(ft,"trend")$p
     meta<-data.frame(Analysis=set,Recall=q$Recall,Start=q$Start,End=q$End,Score=score,Score_label=SCORE_LABELS[score],Model="M2")
     for(g in levels(ds$variables$group)) {
       i<-ds$variables$group==g
       e<-if(g=="T1")data.frame(OR=1,LCI=1,UCI=1,p=NA_real_,df=fg$df.residual) else extract_or(fg,paste0("group",g))
       rows[[paste(set,score,g)]]<-cbind(meta,Contrast=g,Minimum=min(x[i]),Maximum=max(x[i]),N=sum(i),Cases=sum(z$stroke_prevalent[i]),
         Weighted_percent=100*sum(w[i])/sum(w),Increment=NA_real_,e,P_trend=pt)
     }
     sdw<-sqrt(as.numeric(coef(svyvar(~x,ds))))
     for(ct in c("Per point","Per SD")) {
       k<-if(ct=="Per SD")sdw else 1
       rows[[paste(set,score,ct)]]<-cbind(meta,Contrast=ct,Minimum=NA_real_,Maximum=NA_real_,N=nrow(z),Cases=sum(z$stroke_prevalent),
         Weighted_percent=NA_real_,Increment=k,extract_or(fc,"x",k),P_trend=NA_real_)
     }
     z[[paste0(score,"_tertile")]]<-ds$variables$group
     audit[[set]][[score]]<-list(cutpoints=br,group_medians=med,sd=sdw,covariates=M2[[score]],converged=all(fc$converged,fg$converged,ft$converged))
   }
   cohorts[[set]]<-z
 }
 results<-bind_rows(rows);stopifnot(nrow(results)==14*length(SCORES)*5)
 bundle<-list(primary=cohorts$two_day_2003_2018,cohorts=cohorts,M2=M2,score_labels=SCORE_LABELS,results=results,audit=audit,patterns_provenance=patterns$provenance,
   endpoint="MCQ160F: self-reported prevalent any stroke",provenance=provenance,session=sessionInfo(),
   code_md5=tools::md5sum(file.path(dose_script_dir(),c("nhanes.R","nhanes.R","nhanes.R","nhanes.R","score.R","score.R"))),
   ukb_config_md5=tools::md5sum(file.path(ROOT,"ukb/00data/analysis_config.rds")))
 write.csv(results,file.path(OUT,"NHANES_M2_results.csv"),row.names=FALSE,na="")
 dose_save_rds(bundle,file.path(OUT,"NHANES_analysis.rds"));write_nhanes_readme(bundle)
 invisible(results)
}
if(isTRUE(getOption("dose.nhanes.models.run",TRUE)))nhanes_results<-run_nhanes_stroke()
}))

# ===========================================================================
# nhanes / plot
# ===========================================================================
dose_register_step("nhanes", "plot", expression({
# Current DOSE9half: pooled-wave comparison and dynamically calculated flow.
suppressPackageStartupMessages({library(ggplot2);library(dplyr)})
ROOT<-dose_output_dir()
OUT<-file.path(ROOT,'nhanes')
run_nhanes_plots<-function() {
 a<-read.csv(file.path(OUT,'NHANES_M2_results.csv'))
 d<-subset(a,Contrast=='Per SD')
 d$Window<-paste(ifelse(d$Recall=='Day 1','Day 1','2-day'),paste0(d$Start,'-',d$End))
 d$Window<-factor(d$Window,levels=rev(unique(d$Window)))
 d$Score<-d$Score_label
 d$Text<-sprintf('%.3f (%.3f-%.3f)  P=%.3g',d$OR,d$LCI,d$UCI,d$p)
 p<-ggplot(d,aes(OR,Window,color=Recall))+
   geom_vline(xintercept=1,linetype=2,color='grey65')+
   geom_errorbar(aes(xmin=LCI,xmax=UCI),orientation='y',width=.16)+geom_point(size=2.5)+
   geom_text(aes(x=1.18,label=Text),color='grey15',hjust=0,size=3.1)+
   facet_wrap(~Score,ncol=2)+scale_color_manual(values=c('Day 1'='#26658A','Two-day mean'='#B96537'))+
   scale_x_log10(limits=c(.60,2.12),breaks=c(.6,.8,1))+labs(
   title='NHANES: four dietary patterns with prevalent stroke',
   subtitle='M2: adjusted odds ratio per weighted SD | All 14 pooled-wave comparisons',
   x='Odds ratio (95% CI)',y=NULL,
   caption='Same participants and survey weights within each window. DOSE / MIND / aMED exclude alcohol adjustment; DASH includes it.\nDay 1: 1999-2018 available; two-day mean: 2003-2018 available. Exploratory comparisons; no multiple-testing correction.')+
   theme_bw(base_size=11)+theme(legend.position='top',panel.grid.minor=element_blank(),panel.grid.major.y=element_blank(),plot.title=element_text(face='bold'))
 ggsave(file.path(OUT,'Fig_NHANES_M2.png'),p,width=14,height=15,dpi=220,bg='white')
 z<-subset(a,Analysis=='two_day_2003_2018' & Contrast %in% c('T1','T2','T3','Per SD'))
 z$Row<-factor(z$Contrast,levels=c('Per SD','T3','T2','T1'))
 z$Text<-ifelse(z$Contrast=='T1','Reference',sprintf('%.3f (%.3f-%.3f); P=%.3g',z$OR,z$LCI,z$UCI,z$p))
 z$Counts<-paste(z$Cases,'/',z$N)
 z$Score<-z$Score_label
 labs<-sapply(split(z,z$Score),function(x)sprintf('%s | Two-day mean 2003-2018\nP trend = %.4g',x$Score[1],x$P_trend[1]))
 p2<-ggplot(z,aes(OR,Row))+geom_vline(xintercept=1,linetype=2,color='grey65')+
   geom_errorbar(aes(xmin=LCI,xmax=UCI),orientation='y',width=.1,color='#26658A')+geom_point(color='#26658A',size=3)+
   geom_text(aes(x=1.2,label=Text),hjust=0,size=2.8)+geom_text(aes(x=3.3,label=Counts),hjust=0,size=2.8)+
   facet_wrap(~Score,ncol=2,labeller=as_labeller(labs))+scale_x_log10(limits=c(.45,4.6),breaks=c(.5,.75,1))+
   labs(title='NHANES 2003-2018: weighted tertiles and continuous scores',subtitle='M2 | Right column: stroke cases / participants',
   x='Adjusted odds ratio (95% CI)',y=NULL,caption='Tertiles preserve score ties. P trend uses weighted within-group medians. Full results and group boundaries are in the CSV.')+
   theme_bw(base_size=12)+theme(panel.grid.minor=element_blank(),panel.grid.major.y=element_blank())
 ggsave(file.path(OUT,'Fig_NHANES_M2_tertiles.png'),p2,width=14,height=15,dpi=220,bg='white')
 unlink(file.path(OUT,'Fig_NHANES_M2.pdf'))
 invisible(p)
}
run_nhanes_flow_plot<-function() {
 a<-readRDS(file.path(OUT,'NHANES_analysis.rds'))
 keys<-c('day1_1999_2018','two_day_2003_2018')
 draw<-function() {
  grid::grid.newpage()
  grid::grid.text('NHANES participant selection',x=.035,y=.98,just=c('left','top'),gp=grid::gpar(fontsize=20,fontface='bold'))
  for(j in 1:2) {
   f<-a$audit[[keys[j]]];n<-format(f$flow,big.mark=',',trim=TRUE);ex<-format(f$excluded,big.mark=',',trim=TRUE)
   cx<-if(j==1).255 else .745;ys<-seq(.83,.18,length.out=6);bh<-.075
   title<-if(j==1)'Day 1 | 1999-2018' else 'Two-day mean | 2003-2018'
   grid::grid.text(title,x=cx,y=.92,gp=grid::gpar(fontsize=15,fontface='bold',col='#26658A'))
   text<-c(paste0(n[1],' participants in selected cycles'),paste0(n[2],' adults aged >=20 years'),
    paste0(n[3],' with valid dietary recall(s)\nand positive dietary weights'),
    paste0(n[4],' with known stroke status'),paste0(n[5],' with known hypertension\nand diabetes status (including Borderline)'),
    paste0(n[6],' in all FOUR dietary-pattern M2 analyses\n',format(f$cases,big.mark=','),' prevalent stroke cases'))
   reasons<-c('aged <20 / age missing','invalid recall(s) / weights','unknown stroke','unknown hypertension / diabetes','missing score / required covariates')
   for(i in 1:6) {
    grid::grid.roundrect(x=cx,y=ys[i],width=.435,height=bh,r=grid::unit(.01,'snpc'),gp=grid::gpar(fill=if(i==6)'#EAF2F7' else 'white',col='#3C5260'))
    grid::grid.text(text[i],x=cx,y=ys[i],gp=grid::gpar(fontsize=11,fontface=if(i==6)'bold' else 'plain'))
    if(i<6) {
     grid::grid.lines(x=c(cx,cx),y=c(ys[i]-bh/2,ys[i+1]+bh/2),arrow=grid::arrow(length=grid::unit(.07,'inches'),type='closed'))
     grid::grid.text(paste0('Excluded: ',ex[i],' (',reasons[i],')'),x=cx+.015,y=(ys[i]+ys[i+1])/2,
       just='left',gp=grid::gpar(fontsize=7.3,col='#555555'))
    }
   }
  }
  grid::grid.text('Counts are sequential and unweighted. No upper age exclusion. Valid recalls require status=1, positive energy and complete food-code/amount processing.\nOther baseline covariates use the documented imputation. Outcome: self-reported prevalent any stroke (MCQ160F).\nDiabetes Borderline is retained as a separate category; All four patterns share each final cohort.',
   x=.035,y=.075,just=c('left','top'),gp=grid::gpar(fontsize=9,lineheight=1.4))
 }
 png(file.path(OUT,'Sfig1_nhanes.png'),width=15,height=11,units='in',res=200);draw();dev.off()
 unlink(file.path(OUT,c('Sfig1_nhanes.tif','Sfig1_nhanes.tiff')))
}
if(isTRUE(getOption('dose.nhanes.plots.run',TRUE))) {run_nhanes_plots();run_nhanes_flow_plot()}
}))

# ===========================================================================
# nhanes / flow
# ===========================================================================
dose_register_step("nhanes", "flow", expression({
# Compact primary-cohort flow; counts and final IDs checked against fitted data.
ROOT <- dose_output_dir()
a <- readRDS(file.path(ROOT, 'nhanes/NHANES_analysis.rds'))
d <- readRDS(file.path(ROOT, '.cache/nhanes/nhanes_two_day_score.rds'))
d <- d[d$cycle >= '2003' & d$cycle < '2019', ]
stopifnot(!anyNA(d$age), all(d$valid), all(is.finite(d$weight) & d$weight > 0),
          !anyNA(d[c('psu','strata')]), !anyDuplicated(d$eid))
young <- d$age < 20
stroke_missing <- !young & is.na(d$stroke_prevalent)
disease_missing <- !young & !stroke_missing & !complete.cases(d[c('hypertension','diabetes')])
keep <- !young & !stroke_missing & !disease_missing
stopifnot(setequal(d$eid[keep], a$primary$eid), sum(keep) == nrow(a$primary))
n <- function(x) format(x, big.mark=',', trim=TRUE, scientific=FALSE)
cases <- sum(d$stroke_prevalent[keep])
output <- file.path(ROOT, 'nhanes/Sfig1_nhanes_primary.png')
grDevices::png(output, width=10, height=6.3, units='in', res=300, family='Arial', bg='white')
grid::grid.newpage()
text <- function(label,x,y,size=14,bold=FALSE,just='centre')
 grid::grid.text(label,x=x,y=y,just=just,gp=grid::gpar(fontsize=size,fontface=if(bold)'bold' else 'plain',col='#202020',lineheight=1.35))
box <- function(x,y,w,h,fill)
 grid::grid.roundrect(x=x,y=y,width=w,height=h,r=grid::unit(.012,'snpc'),gp=grid::gpar(fill=fill,col='#62676D',lwd=1))
arrow <- function(x,y)
 grid::grid.lines(x=x,y=y,arrow=grid::arrow(length=grid::unit(.10,'inches'),type='closed'),gp=grid::gpar(col='#777C83',fill='#777C83',lwd=1.3))
text('NHANES 2003-2018',.045,.94,18,TRUE,'left')
box(.31,.79,.53,.20,'#E5EDF7')
text(paste0(n(nrow(d)), ' participants'),.31,.83,17,TRUE)
text('with two valid 24-h dietary recalls',.31,.775,14)
text('and valid dietary survey weights',.31,.731,13)
arrow(c(.31,.31),c(.69,.30))
arrow(c(.31,.58),c(.50,.50))
box(.785,.50,.37,.30,'#F6F2DF')
text(paste0('Excluded (n = ',n(sum(!keep)),')'),.619,.60,14,TRUE,'left')
text(paste0(n(sum(young)), ' aged <20 years'),.619,.548,13,FALSE,'left')
text(paste0(n(sum(stroke_missing)), ' with missing stroke status'),.619,.496,13,FALSE,'left')
text(paste0(n(sum(disease_missing)), ' with missing hypertension\n     or diabetes status'),.619,.423,13,FALSE,'left')
box(.31,.20,.53,.20,'#E5EDF7')
text(paste0(n(sum(keep)), ' participants included'),.31,.247,17,TRUE)
text('in the DOSE and DASH analyses',.31,.19,14)
text(paste0(n(cases), ' self-reported stroke cases'),.31,.14,14)
text('Counts are unweighted; exclusions are sequential. No upper age limit was applied.',.045,.045,10,FALSE,'left')
grDevices::dev.off()
unlink(file.path(ROOT,'nhanes/Sfig1_nhanes.pdf'))
cat('Saved:',output,'\n')
cat('Flow:',nrow(d),'-',sum(young),'-',sum(stroke_missing),'-',sum(disease_missing),'=',sum(keep),'; cases:',cases,'\n')
}))

if (sys.nframe() == 0L) dose_run_module("nhanes")
