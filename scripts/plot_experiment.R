#!/usr/bin/env Rscript
library(ggplot2)
library(grid)
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."))
dest <- file.path(root,"analysis","experimental")
systems <- c("reference","conf-01","conf-02","dimer-01-01","dimer-01-02","dimer-02-02")
labels <- c("Reference","conf-01","conf-02","Dimer 01/01","Dimer 01/02","Dimer 02/02")
colors <- setNames(c("#14191D","#3974A2","#CE8133","#697B3E","#9B5F85","#424C57"),systems)
types <- setNames(c("solid","solid","22","longdash","dotdash","13"),systems)
theme_set(theme_classic(base_size=12,base_family="sans") +
 theme(plot.title=element_text(face="bold",size=15,color="#26333E"),
 plot.subtitle=element_text(size=10,color="#59656E",margin=margin(b=9)),
 panel.grid.major.y=element_line(color="#E9ECEF",linewidth=.3),
 axis.line=element_line(color="#54616A",linewidth=.4),
 strip.background=element_blank(),strip.text=element_text(face="bold",hjust=0),
 plot.margin=margin(10,15,10,12),legend.position="top",legend.title=element_blank(),
 plot.caption=element_text(hjust=0,size=9,color="#59656E")))
saveplot <- function(name,p,w=7.2,h=5.7) {
 ggsave(file.path(dest,paste0(name,".pdf")),p,width=w,height=h,device=pdf,family="Helvetica",useDingbats=FALSE)
 png(file.path(dest,paste0(name,".png")),width=w,height=h,units="in",res=210,type=if (capabilities("aqua")) "quartz" else "cairo")
 print(p); dev.off()
}
uv <- read.csv(file.path(dest,"uv-normalized-calculations.csv"))
ref <- read.csv(file.path(dest,"uv-reference.csv"))
ref$system <- "reference"
d <- rbind(uv,ref[,names(uv)])
d$system <- factor(d$system,levels=systems)
p <- ggplot(d,aes(wavelength_nm,normalized_intensity,color=system,linetype=system)) +
 geom_line(linewidth=.75) + geom_point(data=d[d$system=="reference",],size=.9,show.legend=FALSE) +
 scale_color_manual(values=colors,labels=labels) + scale_linetype_manual(values=types,labels=labels) +
 scale_x_continuous(limits=c(200,293),breaks=seq(200,290,10),expand=expansion(mult=0)) +
 scale_y_continuous(limits=c(0,1.04),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=0)) +
 guides(color=guide_legend(nrow=2,byrow=TRUE),linetype=guide_legend(nrow=2,byrow=TRUE)) +
 labs(title="Sorbic acid | UV reference and calculations",
 subtitle="Reference: NIST / Timmons (1966)  |  Theory: TD-PBE0/6-31G(d), SMD(ether)",
 x="Wavelength / nm",y="Individually normalized linear intensity",
 caption="Reference: linearized absorption; solvent unspecified. Theory: Gaussian FWHM = 0.30 eV.")
saveplot("uv-experiment-comparison",p,h=5.5)
nmr <- read.csv(file.path(dest,"nmr-H-relative-shifts.csv"))
nmr$system <- factor(nmr$system,levels=systems)
nmr$site <- factor(nmr$site,levels=c("CH3","H15","H13","H12","H14","OH"))
# Equal-height sticks show positions only; they do not invent experimental intensities.
p <- ggplot(nmr,aes(relative_shift_ppm,color=system)) +
 geom_vline(xintercept=c(0,3.95,4.42,5.34,10.36),linewidth=.3,color="#D5DADD",linetype="dotted") +
 geom_segment(aes(xend=relative_shift_ppm,y=0,yend=.85),linewidth=.65) +
 facet_wrap(~system,ncol=1,labeller=as_labeller(setNames(labels,systems))) +
 scale_color_manual(values=colors,guide="none") +
 scale_x_reverse(limits=c(11,-.4),breaks=c(0,2,4,6,8,10),expand=expansion(mult=0)) +
 scale_y_continuous(limits=c(0,1.08),breaks=NULL,expand=expansion(mult=0)) +
 labs(title="Sorbic acid | proton peak-position comparison",
 subtitle="Reference: DMSO-d6, 600 MHz, 25 C  |  Calculations: GIAO-PBE0/6-31G(d), gas",
 x="Shift relative to methyl / ppm",y=NULL,
 caption="Sticks: peak positions. Dotted lines: MHLW (2011), Fig. 1, reference positions.\nTheory: mean(methyl shielding) - site shielding; one anchor per species.") +
 theme(axis.line.y=element_blank(),panel.grid.major.y=element_blank())
ref_labels <- data.frame(system=factor("reference",levels=systems),
 relative_shift_ppm=c(0,3.95,4.65,5.34,10.36),label=c("CH3","H15","H12/13","H14","OH"))
p <- p + geom_text(data=ref_labels,aes(label=label,y=.99),size=3,color="#26333E",inherit.aes=TRUE)
p <- p + geom_segment(data=ref_labels[ref_labels$label=="H12/13",],
 aes(x=4.42,xend=4.62,y=.85,yend=.93),inherit.aes=FALSE,color="#26333E",linewidth=.3)
saveplot("nmr-H-experiment-comparison",p,h=7.4)
p <- ggplot(nmr[nmr$system!="reference" & nmr$site!="CH3",],
 aes(site,error_ppm,color=system,shape=system)) +
 geom_hline(yintercept=0,color="#6E767D",linewidth=.4) +
 geom_point(position=position_dodge(width=.55),size=2.5) +
 scale_color_manual(values=colors[-1],labels=labels[-1]) +
 scale_shape_manual(values=c(16,17,15,1,2),labels=labels[-1]) +
 labs(title="Sorbic acid | relative proton-shift residuals",
 subtitle="Calculated relative shift minus the DMSO-d6 reference; methyl anchor excluded",
 x=NULL,y="Relative-shift residual / ppm",
 caption="Individual dimer fragments retained. OH is excluded from the alkene-site error metric.\nGas-phase models and solution experiment have different hydrogen-bonding environments.") +
 guides(color=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))
saveplot("nmr-H-reference-residuals",p,h=4.8)

carbon <- read.csv(file.path(dest,"nmr-C-relative-shifts.csv"))
carbon$system <- factor(carbon$system,levels=systems)
carbon$site <- factor(carbon$site,levels=c("CH3",paste0("alkene-",1:4),"COOH"))
carbon_ref <- carbon[carbon$system=="reference",]
p <- ggplot(carbon,aes(relative_shift_ppm,color=system)) +
 geom_vline(xintercept=carbon_ref$relative_shift_ppm,linewidth=.3,color="#D5DADD",linetype="dotted") +
 geom_segment(aes(xend=relative_shift_ppm,y=0,yend=.85),linewidth=.65) +
 facet_wrap(~system,ncol=1,labeller=as_labeller(setNames(labels,systems))) +
 scale_color_manual(values=colors,guide="none") +
 scale_x_reverse(limits=c(163,-6),breaks=seq(0,160,20),expand=expansion(mult=0)) +
 scale_y_continuous(limits=c(0,1.08),breaks=NULL,expand=expansion(mult=0)) +
 labs(title="Sorbic acid | carbon peak-position comparison",
 subtitle="Reference: CDCl3, TMS  |  Calculations: GIAO-PBE0/6-31G(d), gas",
 x="Shift relative to methyl carbon / ppm",y=NULL,
 caption="Sticks: peak positions. Dotted lines: Nishiyama et al. (1991), Table II, SOA.\nTheory: mean(methyl-C shielding) - site shielding; one anchor per species.") +
 theme(axis.line.y=element_blank(),panel.grid.major.y=element_blank())
carbon_labels <- data.frame(system=factor("reference",levels=systems),
 relative_shift_ppm=c(0,114,154.3),label=c("CH3","Alkene carbons","COOH"))
p <- p + geom_text(data=carbon_labels,aes(label=label,y=.99),size=3,color="#26333E")
saveplot("nmr-C-experiment-comparison",p,h=7.4)

p <- ggplot(carbon[carbon$system!="reference" & carbon$site!="CH3",],
 aes(site,error_ppm,color=system,shape=system)) +
 geom_hline(yintercept=0,color="#6E767D",linewidth=.4) +
 geom_point(position=position_dodge(width=.55),size=2.5) +
 scale_color_manual(values=colors[-1],labels=labels[-1]) +
 scale_shape_manual(values=c(16,17,15,1,2),labels=labels[-1]) +
 scale_x_discrete(labels=c(paste("Alkene",1:4),"COOH")) +
 labs(title="Sorbic acid | relative carbon-shift residuals",
 subtitle="Calculated relative shift minus CDCl3 reference; methyl anchor excluded",
 x=NULL,y="Relative-shift residual / ppm",
 caption="Four alkene peaks paired in ascending order within each fragment, not by atom assignment.\nIndividual dimer fragments retained. Gas-phase models versus solution reference.") +
 guides(color=guide_legend(nrow=2,byrow=TRUE),shape=guide_legend(nrow=2,byrow=TRUE))
saveplot("nmr-C-reference-residuals",p,h=4.8)
