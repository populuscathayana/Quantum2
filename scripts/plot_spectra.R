#!/usr/bin/env Rscript
# Deterministic static scientific figures. No experimental spectrum is inferred.
library(ggplot2)
library(grid)
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- normalizePath(file.path(dirname(script), ".."))
dest <- file.path(root, "analysis", "stage2")
uv <- read.csv(file.path(dest, "uv-transitions.csv"))
nmr <- read.csv(file.path(dest, "nmr-shieldings.csv"))
blue <- "#3974A2"; orange <- "#CE8133"; ink <- "#26333E"
theme_set(theme_classic(base_size=12, base_family="sans") +
  theme(plot.title=element_text(face="bold", size=15, color=ink),
        plot.subtitle=element_text(size=10.5, color="#59656E", margin=margin(b=12)),
        axis.title=element_text(color=ink), axis.text=element_text(color=ink),
        panel.grid.major.y=element_line(color="#E9ECEF", linewidth=.3),
        axis.line=element_line(color="#54616A", linewidth=.4),
        axis.ticks=element_line(color="#54616A", linewidth=.4),
        plot.margin=margin(12,20,12,12)))

save_panels <- function(name, panels, heights, width=7.2, height=5.45) {
  draw_all <- function() {
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(length(panels), 1, heights=unit(heights,"null"))))
    for (i in seq_along(panels)) print(panels[[i]], vp=viewport(layout.pos.row=i))
  }
  pdf(file.path(dest,paste0(name,".pdf")), width=width,height=height,family="Helvetica",useDingbats=FALSE)
  # Use a core PDF font for portable vector exports; PNG uses the platform sans-serif face.
  old <- theme_get(); theme_update(text=element_text(family="Helvetica"))
  # Plot objects captured their themes, so override text on each object as well.
  originals <- panels
  panels <- lapply(panels,function(p) p + theme(text=element_text(family="Helvetica")))
  draw_all(); dev.off(); theme_set(old); panels <- originals
  png(file.path(dest,paste0(name,".png")), width=width,height=height,units="in",res=240,type=if (capabilities("aqua")) "quartz" else "cairo")
  draw_all(); dev.off()
}

lambda <- seq(120,400,by=.025)
signal <- vapply(1239.841984/lambda, function(E) sum(uv$oscillator_strength *
  exp(-4*log(2)*((E-uv$energy_ev)/.30)^2)), 0)
curve <- data.frame(wavelength_nm=lambda, intensity=signal/max(signal))
write.csv(curve, file.path(dest,"uv-broadened.csv"),row.names=FALSE)
p1 <- ggplot(curve,aes(wavelength_nm,intensity)) +
  geom_area(fill=blue,alpha=.12) + geom_line(color=blue,linewidth=.8) +
  annotate("text",x=279,y=.92,label="S2   257.44 nm",hjust=0,size=4,color=ink) +
  annotate("segment",x=259,xend=277,y=.92,yend=.92,color="#64717B",linewidth=.3) +
  scale_x_continuous(limits=c(120,400), breaks=seq(120,400,40),expand=expansion(mult=0)) +
  scale_y_continuous(limits=c(0,1.08),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=0)) +
  labs(title="Sorbic acid | UV absorption", subtitle="conf-02  /  TD-PBE0/6-31G(d), SMD(diethyl ether)  /  20 singlet states",
       x=NULL,y="Normalized amplitude")
p2 <- ggplot(uv,aes(wavelength_nm,oscillator_strength)) +
  geom_segment(aes(xend=wavelength_nm,yend=0),color=orange,linewidth=.65) +
  annotate("text",x=279,y=.77,label="S2   f = 0.9375",hjust=0,size=3.6,color=ink) +
  annotate("text",x=295,y=.22,label="S1   262.36 nm\nf = 0.0001",hjust=0,size=3.5,color=ink) +
  annotate("segment",x=262.36,xend=292,y=.01,yend=.17,color="#64717B",linewidth=.3) +
  scale_x_continuous(limits=c(120,400), breaks=seq(120,400,40),expand=expansion(mult=0)) +
  scale_y_continuous(limits=c(0,1.02),breaks=c(0,.25,.5,.75,1),expand=expansion(mult=0)) +
  labs(subtitle="Energy-domain Gaussian FWHM = 0.30 eV (top); discrete transitions (bottom)",
       x="Wavelength / nm",y="Oscillator strength, f")
save_panels("uv-spectrum",list(p1,p2),c(1.15,1))

nmrplot <- function(elem, limits, width, color, title, ticks, labels) {
  d <- nmr[nmr$element==elem,]
  x <- seq(limits[1],limits[2],length.out=12000)
  y <- vapply(x,function(s) sum((width/2)/pi/((s-d$shielding_ppm)^2+(width/2)^2)),0)
  dcurve <- data.frame(shielding_ppm=x,intensity=y/max(y))
  write.csv(dcurve,file.path(dest,paste0("nmr-",elem,"-broadened.csv")),row.names=FALSE)
  ggplot(dcurve,aes(shielding_ppm,intensity)) +
    geom_line(color=color,linewidth=.8) +
    geom_segment(data=d,aes(x=shielding_ppm,xend=shielding_ppm,y=0,yend=-.055),color=color,linewidth=.5) +
    geom_text(data=labels,aes(x=x,y=y,label=label),inherit.aes=FALSE,color=ink,size=3.6) +
    scale_x_continuous(limits=limits,breaks=ticks,expand=expansion(mult=.005)) +
    scale_y_continuous(limits=c(-.075,1.21),breaks=c(0,.5,1),expand=expansion(mult=0)) +
    labs(title=title,subtitle=paste0("conf-02  /  GIAO-PBE0/6-31G(d), gas phase  /  Lorentzian FWHM = ",width," ppm"),
         x=expression(paste("Isotropic shielding, ",sigma," / ppm")), y="Normalized amplitude")
}
hlabels <- data.frame(x=c(24.62,25.80,26.49,27.08,30.22), y=c(.60,1.10,.60,.85,1.10),
                      label=c("H14","H12 / H13","H15","H16 (OH)","H9-H11 (CH3)"))
clabels <- data.frame(x=c(37.48,51.712,59.3,68.62,81.03,173.83), y=c(1.12,1.12,.84,1.12,1.12,1.12),
                      label=c("C6","C4","C2","C3","C5","C1"))
ph <- nmrplot("H",c(23.8,31.2),.08,blue,"Sorbic acid | proton shieldings",seq(24,31,1),hlabels)
pc <- nmrplot("C",c(25,185),.80,orange,"Sorbic acid | carbon-13 shieldings",seq(40,180,20),clabels) +
  annotate("segment",x=56.5682,xend=58.2,y=.92,yend=.87,color="#64717B",linewidth=.3) +
  labs(caption="Equal-area Lorentzian broadening of static isotropic shieldings.") +
  theme(plot.caption=element_text(hjust=0,size=10,color="#59656E"))
save_panels("nmr-shielding-spectra",list(ph,pc),c(1,1.06),height=5.9)
message("Saved UV and NMR figures (PNG + vector PDF).")
