# Run independently after regression; reads mri_lm_results.csv and saves PDF heatmaps.
# Anatomical plots use s1_5_mri_lm.plot.py separately.
# Colour: MRI SD per 1-point higher FLT; stars: BH within FLT and MRI family.
local({
  project <- Sys.getenv("FLT_PROJECT", unset="D:/projects/flt")
  out <- file.path(project, "analysis/01lm/mri")
  fig <- file.path(out, "figures")
  dir.create(fig, recursive=TRUE, showWarnings=FALSE)
  r <- read.csv(file.path(out,"mri_lm_results.csv"), check.names=FALSE)
  required <- c("exposure","food_liking","field_id","family","label","beta","FDR_family","status")
  stopifnot(all(required %in% names(r)), !anyDuplicated(r[c("exposure","field_id")]))
  exposures <- c("p20628","F-Salad_vegetables","p20638")
  foods <- c("Butter on bread","Total salad vegetable","Chicken")
  families <- c("ThicknessDK","FA27","MD27","Volume139","WMH")
  stopifnot(setequal(unique(r$exposure), exposures), all(families %in% r$family))
  r$valid <- r$status=="ok" & is.finite(r$beta) & is.finite(r$FDR_family)
  if (!any(r$valid)) stop("No valid MRI estimates to plot")
  # One symmetric colour scale across every family/page; no coefficient clipping.
  limit <- max(abs(r$beta[r$valid]))
  if (limit==0) limit <- 0.01
  palette <- colorRampPalette(c("#2166AC","#FFFFFF","#B2182B"))(201)
  short_label <- function(x) {
    x <- sub("^Volume of grey matter in ","",x)
    x <- sub("^Mean thickness of ","",x)
    x <- sub("^Mean FA in ","",x)
    x <- sub("^Mean MD in ","",x)
    x
  }
  draw <- function(d, family, page) {
    ids <- unique(d$field_id)
    labs <- short_label(d$label[match(ids,d$field_id)])
    if (family %in% c("ThicknessDK","FA27","MD27") && "atlas_key" %in% names(d)) {
      keys <- d$atlas_key[match(ids,d$field_id)]
      use <- !is.na(keys) & nzchar(keys)
      labs[use] <- keys[use]
    }
    labs <- vapply(labs,function(x) paste(strwrap(x,width=65),collapse="\n"),character(1))
    nr <- length(ids)
    layout(matrix(1:2,1), widths=c(7,1.5))
    par(mar=c(5,27,4,1), family="sans")
    plot(NA, xlim=c(.5,3.5), ylim=c(.5,nr+.5), axes=FALSE,
         xlab="",ylab="",xaxs="i",yaxs="i")
    for (j in seq_along(ids)) for (i in seq_along(exposures)) {
      v <- d[d$field_id==ids[j] & d$exposure==exposures[i],,drop=FALSE]
      y <- nr+1-j
      valid <- nrow(v)==1L && isTRUE(v$valid)
      color <- if(valid) palette[101+round(100*v$beta/limit)] else "#DDDDDD"
      rect(i-.5,y-.5,i+.5,y+.5,col=color,border="white")
      mark <- if(!valid) "NA" else if(v$FDR_family<.001) "***" else if(v$FDR_family<.01) "**" else if(v$FDR_family<.05) "*" else ""
      text(i,y,mark,cex=.9,col=if(valid && abs(v$beta)>limit*.65) "white" else "black")
    }
    axis(1,at=1:3,labels=foods,tick=FALSE,cex.axis=.8)
    axis(2,at=nr:1,labels=labs,las=1,tick=FALSE,cex.axis=.72)
    title(main=sprintf("FLT and MRI | %s | page %d",family,page),cex.main=1)
    mtext("* q<0.05; ** q<0.01; *** q<0.001; grey/NA: unavailable",side=1,line=2.5,cex=.7)
    mtext("BH within each FLT and full MRI family; covs.all + ICV for volumes/WMH",side=1,line=3.5,cex=.65)
    par(mar=c(5,1,4,5))
    plot(NA,xlim=c(0,1),ylim=c(-limit,limit),axes=FALSE,xlab="",ylab="",yaxs="i")
    edges <- seq(-limit,limit,length.out=202)
    rect(0,edges[-202],.35,edges[-1],col=palette,border=NA)
    axis(4,at=seq(-limit,limit,length.out=5),labels=formatC(seq(-limit,limit,length.out=5),digits=3,format="f"),pos=.35,las=1,cex.axis=.75)
    mtext("MRI SD / FLT point",side=3,line=1,cex=.7)
  }
  manifest <- list()
  for (family in families) {
    d <- r[r$family==family,]
    ids <- sort(unique(d$field_id))
    for (page in seq_len(ceiling(length(ids)/30))) {
      selected <- ids[seq.int((page-1)*30+1,min(page*30,length(ids)))]
      chunk <- d[order(match(d$field_id,selected)),]
      chunk <- chunk[chunk$field_id %in% selected,]
      stem <- sprintf("mri_%s_heatmap_%02d",family,page)
      height <- max(4,2.2+.27*length(selected))
      pdf(file.path(fig,paste0(stem,".pdf")),width=14,height=height)
      tryCatch(draw(chunk,family,page),finally=dev.off())
      chunk$plot_page <- stem
      manifest[[length(manifest)+1L]] <- chunk
    }
  }
  write.csv(do.call(rbind,manifest),file.path(fig,"mri_heatmap_data.csv"),row.names=FALSE,na="")
  message("Saved MRI heatmaps: ",fig)
})
