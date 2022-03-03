setwd('./UO Hallett/Projects/usda-climvar')
source('./Competition/model-fit/AM_partitioning.R')

library(plotrix)


### Prepare data for plotting ----

## AVFA
avfa_eps_plot <- as.data.frame(rbind(invader_no_mech$avfa, ir_eps_0$avfa, ir_eps_lamb$avfa, 
                                     ir_eps_alpha$avfa, ir_eps_int$avfa))
avfa_eps_plot$partition <- c("total","eps_0","eps_lamb","eps_alpha","eps_int")
avfa_eps_plot$invader <- "avfa"

## BRHO
brho_eps_plot <- as.data.frame(rbind(invader_no_mech$brho, ir_eps_0$brho, ir_eps_lamb$brho, 
                                     ir_eps_alpha$brho, ir_eps_int$brho))
brho_eps_plot$partition <- c("total","eps_0","eps_lamb","eps_alpha","eps_int")
brho_eps_plot$invader <- "brho"

## VUMY
vumy_eps_plot <- as.data.frame(rbind(invader_no_mech$vumy, ir_eps_0$vumy, ir_eps_lamb$vumy, 
                                     ir_eps_alpha$vumy, ir_eps_int$vumy))
vumy_eps_plot$partition <- c("total","eps_0","eps_lamb","eps_alpha","eps_int")
vumy_eps_plot$invader <- "vumy"

## LACA
laca_eps_plot <- as.data.frame(rbind(invader_no_mech$laca, ir_eps_0$laca, ir_eps_lamb$laca, 
                                     ir_eps_alpha$laca, ir_eps_int$laca))
laca_eps_plot$partition <- c("total","eps_0","eps_lamb","eps_alpha","eps_int")
laca_eps_plot$invader <- "laca"

## ESCA
esca_eps_plot <- as.data.frame(rbind(invader_no_mech$esca, ir_eps_0$esca, ir_eps_lamb$esca, 
                                     ir_eps_alpha$esca, ir_eps_int$esca))
esca_eps_plot$partition <- c("total","eps_0","eps_lamb","eps_alpha","eps_int")
esca_eps_plot$invader <- "esca"

## TRHI
trhi_eps_plot <- as.data.frame(rbind(invader_no_mech$trhi, ir_eps_0$trhi, ir_eps_lamb$trhi, 
                                     ir_eps_alpha$trhi, ir_eps_int$trhi))
trhi_eps_plot$partition <- c("total","eps_0","eps_lamb","eps_alpha","eps_int")
trhi_eps_plot$invader <- "trhi"


### AVFA pairwise partitioning ----

pdf("./Competition/Figures/AVFA_pairwise_partitioning.pdf", height = 16, width = 6)
par(mfrow=c(4,2))

#vs BRHO
barplot(avfa_eps_plot$brho, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "AVFA into BRHO",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(brho_eps_plot$avfa, ylim = c(-4, 2.5), xlab = "Mechanism", main = "BRHO into AVFA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs VUMY
barplot(avfa_eps_plot$vumy, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "AVFA into VUMY",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(vumy_eps_plot$avfa, ylim = c(-4, 2.5), xlab = "Mechanism", main = "VUMY into AVFA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs LACA
barplot(avfa_eps_plot$laca, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "AVFA into LACA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(laca_eps_plot$avfa, ylim = c(-4, 2.5), xlab = "Mechanism", main = "LACA into AVFA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs ESCA
barplot(avfa_eps_plot$esca, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "AVFA into ESCA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(esca_eps_plot$avfa, ylim = c(-6, 2.5), xlab = "Mechanism", main = "ESCA into AVFA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

dev.off()
### BRHO pairwise partitioning ----

pdf("./Competition/Figures/BRHO_pairwise_partitioning.pdf", height = 16, width = 6)
par(mfrow=c(4,2))

#vs AVFA
barplot(brho_eps_plot$avfa, ylim = c(-3.5, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "BRHO into AVFA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(avfa_eps_plot$brho, ylim = c(-3.5, 2.5), xlab = "Mechanism", main = "AVFA into BRHO",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs VUMY
barplot(brho_eps_plot$vumy, ylim = c(-3.5, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "BRHO into VUMY",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(vumy_eps_plot$brho, ylim = c(-3.5, 2.5), xlab = "Mechanism", main = "VUMY into BRHO",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs LACA
barplot(brho_eps_plot$laca, ylim = c(-3.5, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "BRHO into LACA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(laca_eps_plot$brho, ylim = c(-3.5, 2.5), xlab = "Mechanism", main = "LACA into BRHO",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs ESCA
barplot(brho_eps_plot$esca, ylim = c(-3.5, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "BRHO into ESCA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(esca_eps_plot$brho, ylim = c(-3.5, 2.5), xlab = "Mechanism", main = "ESCA into BRHO",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

dev.off()

### VUMY pairwise partitioning ----

pdf("./Competition/Figures/VUMY_pairwise_partitioning.pdf", height = 16, width = 6)
par(mfrow=c(4,2))

#vs AVFA
barplot(vumy_eps_plot$avfa, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "VUMY into AVFA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(avfa_eps_plot$vumy, ylim = c(-4, 2.5), xlab = "Mechanism", main = "AVFA into VUMY",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs BRHO
barplot(vumy_eps_plot$brho, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "VUMY into BRHO",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(brho_eps_plot$vumy, ylim = c(-4, 2.5), xlab = "Mechanism", main = "BRHO into VUMY",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs LACA
barplot(vumy_eps_plot$laca, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "VUMY into LACA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(laca_eps_plot$vumy, ylim = c(-4, 2.5), xlab = "Mechanism", main = "LACA into VUMY",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs ESCA
barplot(vumy_eps_plot$esca, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "VUMY into ESCA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(esca_eps_plot$vumy, ylim = c(-6, 2.5), xlab = "Mechanism", main = "ESCA into VUMY",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

dev.off()
### LACA pairwise partitioning ----

pdf("./Competition/Figures/LACA_pairwise_partitioning.pdf", height = 16, width = 6)
par(mfrow=c(4,2))

#vs AVFA
barplot(laca_eps_plot$avfa, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "LACA into AVFA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(avfa_eps_plot$laca, ylim = c(-4, 2.5), xlab = "Mechanism", main = "AVFA into LACA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs BRHO
barplot(laca_eps_plot$brho, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "LACA into BRHO",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(brho_eps_plot$laca, ylim = c(-4, 2.5), xlab = "Mechanism", main = "BRHO into LACA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs VUMY
barplot(laca_eps_plot$vumy, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "LACA into VUMY",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(vumy_eps_plot$laca, ylim = c(-4, 2.5), xlab = "Mechanism", main = "VUMY into LACA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs ESCA
barplot(laca_eps_plot$esca, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "LACA into ESCA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(esca_eps_plot$laca, ylim = c(-6, 2.5), xlab = "Mechanism", main = "ESCA into LACA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

dev.off()
### ESCA pairwise partitioning ----

pdf("./Competition/Figures/ESCA_pairwise_partitioning.pdf", height = 16, width = 6)
par(mfrow=c(4,2))

#vs AVFA
barplot(esca_eps_plot$avfa, ylim = c(-6, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "ESCA into AVFA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(avfa_eps_plot$esca, ylim = c(-6, 2.5), xlab = "Mechanism", main = "AVFA into ESCA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs BRHO
barplot(esca_eps_plot$brho, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "ESCA into BRHO",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(brho_eps_plot$esca, ylim = c(-4, 2.5), xlab = "Mechanism", main = "BRHO into ESCA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs VUMY
barplot(esca_eps_plot$vumy, ylim = c(-5, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "ESCA into VUMY",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(vumy_eps_plot$esca, ylim = c(-5, 2.5), xlab = "Mechanism", main = "VUMY into ESCA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

#vs LACA
barplot(esca_eps_plot$laca, ylim = c(-4, 2.5), ylab = "Partition of GRWR", xlab = "Mechanism", main = "ESCA into LACA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

barplot(laca_eps_plot$esca, ylim = c(-4, 2.5), xlab = "Mechanism", main = "LACA into ESCA",
        names.arg = c(expression(paste(italic("r")["i"],"-",italic("r")["r"])),expression(Delta^0),expression(Delta^lambda),expression(Delta^alpha),expression(Delta^"int")),
        col = c("#f07e41","gray45","gray45","gray45","gray45"))

dev.off()
### Variation-independent mechanism comparisons ----
eps_0_plot <- as.data.frame(t(as.data.frame(ir_eps_0)))
eps_0_plot$invader <- rownames(eps_0_plot)
                            
eps_0_plot <- eps_0_plot %>%
        pivot_longer(cols = brho:trhi, names_to = "resident", values_to = "GRWR")
eps_0_plot$resident[eps_0_plot$invader == "brho"] <- colnames(brho_eps_plot)[1:5]
eps_0_plot$resident[eps_0_plot$invader == "vumy"] <- colnames(vumy_eps_plot)[1:5]
eps_0_plot$resident[eps_0_plot$invader == "laca"] <- colnames(laca_eps_plot)[1:5]
eps_0_plot$resident[eps_0_plot$invader == "esca"] <- colnames(esca_eps_plot)[1:5]
eps_0_plot$resident[eps_0_plot$invader == "trhi"] <- colnames(trhi_eps_plot)[1:5]

eps_0_plot <- eps_0_plot %>%
        mutate(invader = recode(invader, avfa = "A.fatua", brho = "B.hordeaceus", vumy = "V.myuros",
                                laca = "L.californica", trhi = "T.hirtum", esca = "E.californica"),
               resident = recode(resident, avfa = "A.fatua", brho = "B.hordeaceus", vumy = "V.myuros",
                                laca = "L.californica", trhi = "T.hirtum", esca = "E.californica")) %>%
        mutate(invader = factor(invader, levels = c("A.fatua", "V.myuros", "B.hordeaceus", 
                                                    "T.hirtum", "L.californica", "E.californica")))

pdf('./Competition/Figures/figures/VarInd_LDGR.pdf', height = 6.5, width = 5.5)
ggplot(eps_0_plot, aes(x = GRWR, y = resident, fill = resident)) + 
        geom_bar(stat = "identity") +
        xlab("Invader LDGR - Variation Independent Mechanisms") +  
        geom_vline(xintercept = 0, linetype = "dashed") + 
        facet_grid(rows = vars(invader)) + 
        ylab("") +
        scale_fill_brewer(palette = "Accent") +
        theme_bw() + 
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
dev.off()

### All average mechanism comparison ----

## Just resident
mech <- data.frame(mechanism = factor(c("Variation-Independent","RNL in Lambda","RNL in Alpha","E-C Covariance"),
                                          levels = c("Variation-Independent","RNL in Alpha","RNL in Lambda","E-C Covariance")),
                       average = c(mean(unlist(coexist_eps_0)),
                                   mean(unlist(coexist_eps_lamb)),
                                   mean(unlist(coexist_eps_alpha)),
                                   mean(unlist(coexist_eps_int))),
                       se = c(std.error(unlist(coexist_eps_0)),
                              std.error(unlist(coexist_eps_lamb)),
                              std.error(unlist(coexist_eps_alpha)),
                              std.error(unlist(coexist_eps_int))))

pdf('./Competition/Figures/mechanism_contribution.pdf', width = 4, height = 4)
ggplot(mech, aes(x = 1:4, y = average)) + 
        geom_point(size = 2) + 
        geom_errorbar(aes(ymin = average - 2*se, ymax = average + 2*se), width = 0.2) + 
        geom_hline(yintercept = 0, linetype = "dashed") + 
        ylab("Mean Contribution to Log LDGR") + xlab("Mechanism") +
        scale_x_continuous(breaks = 1:4, labels = c(expression(Delta^"0"), expression(Delta^lambda), expression(Delta^alpha),expression(Delta^"cov(E,C)"))) + 
        coord_cartesian(xlim = c(0.5,4.5)) +
        theme_classic() +
        theme(axis.text.x = element_text(size = 14))
dev.off()

## Invader - Resident comparison

mech <- data.frame(mechanism = factor(c("Variation-Independent","RNL in Lambda","RNL in Alpha","E-C Covariance"),
                                      levels = c("Variation-Independent","RNL in Alpha","RNL in Lambda","E-C Covariance")),
                   average = c(mean(unlist(ir_eps_0)),
                               mean(unlist(ir_eps_lamb)),
                               mean(unlist(ir_eps_alpha)),
                               mean(unlist(ir_eps_int))),
                   se = c(std.error(unlist(ir_eps_0)),
                          std.error(unlist(ir_eps_lamb)),
                          std.error(unlist(ir_eps_alpha)),
                          std.error(unlist(ir_eps_int))))

pdf('./Competition/Figures/figures/ir_mechanism_contribution.pdf', width = 4, height = 4)
ggplot(mech, aes(x = 1:4, y = average)) + 
        geom_point(size = 2) + 
        geom_errorbar(aes(ymin = average - 2*se, ymax = average + 2*se), width = 0.2) + 
        geom_hline(yintercept = 0, linetype = "dashed") + 
        ylab("Mean Contribution to Log LDGR") + xlab("Mechanism") +
        scale_x_continuous(breaks = 1:4, labels = c(expression(Delta^"0"), expression(Delta^lambda), expression(Delta^alpha),expression(Delta^paste(alpha, lambda)))) + 
        coord_cartesian(xlim = c(0.5,4.5)) +
        theme_classic() +
        theme(axis.text.x = element_text(size = 14))
dev.off()

### Variation-dependent average mechanism comparison ----

var_mech <- data.frame(mechanism = factor(c("RNL in Lambda","RNL in Alpha","E-C Covariance"),
                                          levels = c("RNL in Alpha","RNL in Lambda","E-C Covariance")),
                       average = c(mean(unlist(coexist_eps_lamb)),
                                   mean(unlist(coexist_eps_alpha)),
                                   mean(unlist(coexist_eps_int))),
                       se = c(std.error(unlist(coexist_eps_lamb)),
                              std.error(unlist(coexist_eps_alpha)),
                              std.error(unlist(coexist_eps_int))))

pdf('./Competition/Figures/VarDep_contribution.pdf', width = 3, height = 4)
ggplot(var_mech, aes(x = 1:3, y = average)) + 
        geom_point(size = 2) + 
        geom_errorbar(aes(ymin = average - 2*se, ymax = average + 2*se), width = 0.2) + 
        geom_hline(yintercept = 0, linetype = "dashed") + 
        ylab("Mean Contribution to Log LDGR") + xlab("Mechanism") +
        scale_x_continuous(breaks = 1:3, labels = c(expression(Delta^lambda), expression(Delta^alpha),expression(Delta^"cov(E,C)"))) + 
        coord_cartesian(xlim = c(0.5,3.5)) +
        theme_classic() +
        theme(axis.text.x = element_text(size = 14))
dev.off()


