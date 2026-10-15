#!/usr/bin/env Rscript
library(ggplot2)
library(grid)
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."))
dest <- file.path(root,"analysis","comparison")
systems <- c("conf-01","conf-02","dimer-01-01","dimer-01-02","dimer-02-02")
labels <- c("conf-01","conf-02","Dimer 01/01","Dimer 01/02","Dimer 02/02")
colors <- setNames(c("#3974A2","#CE8133","#697B3E","#9B5F85","#424C57"),systems)
types <- setNames(c("solid","22","longdash","dotdash","13"),systems)
ink <- "#26333E"
theme_set(theme_classic(base_size=12,base_family="sans") +
 theme(plot.title=element_text(face="bold",size=15,color=ink),
 plot.subtitle=element_text(size=10,color="#59656E",margin=margin(b=9)),
 panel.grid.major.y=element_line(color="#E9ECEF",linewidth=.3),
 axis.line=element_line(color="#54616A",linewidth=.4),
 strip.background=element_blank(),strip.text=element_text(face="bold",hjust=0),
 plot.margin=margin(10,15,10,12),legend.position="top",legend.title=element_blank()))
saveplot <- function(name,p,w=7.2,h=5.7) {
 ggsave(file.path(dest,paste0(name,".pdf")),p,width=w,height=h,device=pdf,family="Helvetica",useDingbats=FALSE)
 png(file.path(dest,paste0(name,".png")),width=w,height=h,units="in",res=210,type=if (capabilities("aqua")) "quartz" else "cairo")
 print(p); dev.off()
}
uv <- read.csv(file.path(dest,"uv-curves.csv"))
uv$system <- factor(uv$system,levels=systems)
full <- uv; full$region <- "160-350 nm"
zoom <- uv[uv$wavelength_nm >= 235 & uv$wavelength_nm <= 290,];zoom$region <- "235-290 nm"
d <- rbind(full,zoom); d$region <- factor(d$region,levels=c("160-350 nm","235-290 nm"))
p <- ggplot(d,aes(wavelength_nm,amplitude_per_acid_unit,color=system,linetype=system)) +
 geom_line(linewidth=.75) + facet_wrap(~region,ncol=1,scales="free_x") +
 scale_color_manual(values=colors,labels=labels,drop=FALSE) + scale_linetype_manual(values=types,labels=labels,drop=FALSE) +
 scale_y_continuous(limits=c(0,1.06),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=0)) +
 scale_x_continuous(expand=expansion(mult=c(.01,.01))) +
 guides(color=guide_legend(nrow=2,byrow=TRUE),linetype=guide_legend(nrow=2,byrow=TRUE)) +
 labs(title="Sorbic acid | monomer and dimer UV spectra",
 subtitle="TD-PBE0/6-31G(d), SMD(diethyl ether)  |  Gaussian FWHM = 0.30 eV",
 x="Wavelength / nm",y="Amplitude per acid unit / a.u.",
 caption="Amplitudes per acid unit: dimer intensities divided by two.") +
 theme(plot.caption=element_text(hjust=0,size=9,color="#59656E"))
saveplot("uv-comparison",p,h=6.6)
nmr <- read.csv(file.path(dest,"nmr-shieldings.csv"))
for (element in c("H","C")) {
 width <- if(element=="H") .08 else .80
 limits <- if(element=="H") c(18.5,31.5) else c(25,180)
 curves <- do.call(rbind,lapply(systems,function(name){
   nuclei <- nmr[nmr$system==name & nmr$element==element,]
   x <- seq(limits[1],limits[2],length.out=12000)
   y <- vapply(x,function(s)sum((width/2)/pi/((s-nuclei$shielding_ppm)^2+(width/2)^2)),0)/nuclei$units[1]
   data.frame(system=name,shielding_ppm=x,amplitude=y)
 }))
 curves$intensity <- curves$amplitude/max(curves$amplitude)
 curves$system <- factor(curves$system,levels=systems)
 sticks <- nmr[nmr$element==element,]; sticks$system <- factor(sticks$system,levels=systems)
 p <- ggplot(curves,aes(shielding_ppm,intensity,color=system)) + geom_line(linewidth=.65) +
   geom_segment(data=sticks,aes(x=shielding_ppm,xend=shielding_ppm,y=0,yend=-.06),linewidth=.4) +
   facet_wrap(~system,ncol=1,labeller=as_labeller(setNames(labels,systems))) +
   scale_color_manual(values=colors,guide="none") +
   scale_x_continuous(limits=limits,breaks=if(element=="H")seq(20,30,2) else seq(30,180,30),expand=expansion(mult=.004)) +
   scale_y_continuous(limits=c(-.075,1.06),breaks=c(0,.5,1),expand=expansion(mult=0)) +
   labs(title=paste("Sorbic acid |",if(element=="H")"proton shielding comparison" else "carbon-13 shielding comparison"),
    subtitle=paste0("GIAO-PBE0/6-31G(d), gas phase  |  Lorentzian FWHM = ",width," ppm"),
    x=expression(paste("Isotropic shielding, ",sigma," / ppm")),y="Shared-scale amplitude",
    caption="Equal area per nucleus / acid unit; one common intensity scale.") +
   theme(plot.caption=element_text(hjust=0,size=9,color="#59656E"))
 saveplot(paste0("nmr-",element,"-comparison"),p,h=7.4)
 write.csv(curves,file.path(dest,paste0("nmr-",element,"-curves.csv")),row.names=FALSE)
}
atoms <- read.csv(file.path(dest,"dimer-atoms.csv")); bonds <- read.csv(file.path(dest,"dimer-bonds.csv"))
atoms$system <- factor(atoms$system,levels=systems[3:5]); bonds$system <- factor(bonds$system,levels=systems[3:5])
p <- ggplot() + geom_segment(data=bonds[bonds$kind=="covalent",],aes(x,y,xend=xend,yend=yend),color="#66717B",linewidth=.6) +
 geom_segment(data=bonds[bonds$kind=="hydrogen",],aes(x,y,xend=xend,yend=yend),color="#3974A2",linetype="dashed",linewidth=.8) +
 geom_point(data=atoms,aes(x,y,fill=element,size=element),shape=21,color="#35414A",stroke=.3) +
 scale_fill_manual(values=c(C="#5B6570",H="#FFFFFF",O="#CE8133")) +
 scale_size_manual(values=c(C=2.8,H=1.8,O=3.8),guide="none") +
 facet_wrap(~system,ncol=1,labeller=as_labeller(setNames(labels,systems))) + coord_fixed() +
 labs(title="Sorbic acid | optimized cyclic dimers",subtitle="PBE0/6-31G(d), gas phase  |  Projection onto the molecular plane",fill=NULL,
 caption="Dashed lines: O-H...O hydrogen bonds. Atom positions from optimized geometries.") +
 theme_void(base_size=12,base_family="Helvetica") + theme(legend.position="top",strip.text=element_text(face="bold",size=11),
 plot.title=element_text(face="bold",size=15),plot.subtitle=element_text(size=10,margin=margin(b=6)),
 plot.caption=element_text(size=9,hjust=0),plot.title.position="plot",plot.caption.position="plot",
 plot.margin=margin(12,15,12,15))
saveplot("dimer-structures",p,h=5.8)
message("Saved four comparison figures as vector PDF and PNG.")
