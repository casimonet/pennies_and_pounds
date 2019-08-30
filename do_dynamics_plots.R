setwd('/path/to/simulations_output.RData')


mtexti <- function(text, side, off = 0.25,
                   srt = if(side == 2) 90  else
                     if(side == 4) 270 else 0, ...) {
  # dimensions of plotting region in user units
  usr <- par('usr')
  # dimensions of plotting region in inches
  pin <- par('pin')
  # user units per inch
  upi <- c(usr[2]-usr[1],
           usr[4]-usr[3]) / pin
  # default x and y positions
  xpos <- (usr[1] + usr[2])/2
  ypos <- (usr[3] + usr[4])/2
  if(1 == side)
    ypos <- usr[3] - upi[2] * off
  if(2 == side)
    xpos <- usr[1] - upi[1] * off
  if(3 == side)
    ypos <- usr[4] + upi[2] * off
  if(4 == side)
    xpos <- usr[2] + upi[1] * off
  text(x=xpos, y=ypos, text, xpd=NA, srt=srt, ...)
}


library(dplyr)
files<- list.files()


pdf('/path/to/output/dynamics_example.pdf', width = 9.5, height = 7.5)


for(j in 1:length(files)){
load(files[j])

for(i in 1:length(runs)){

  
  if(out_data[i,'prop_f_usage'] >= 0.5){
  # do ----
  t0<- data.frame(time = outs_burn_in[[i]][,'time'],
                  C_00 = rowSums(outs_burn_in[[i]][,c('C_00_u', 'C_00_f')]),
                  C_10 = rowSums(outs_burn_in[[i]][,c('C_10_u', 'C_10_f')]),
                  C_01 = rowSums(outs_burn_in[[i]][,c('C_01_u', 'C_01_f')]),
                  C_11 = rowSums(outs_burn_in[[i]][,c('C_11_u', 'C_11_f')]),
                  I_00 = rowSums(outs_burn_in[[i]][,c('I_00_l', 'I_00_f')]),
                  I_10 = rowSums(outs_burn_in[[i]][,c('I_10_l', 'I_10_f')]),
                  I_01 = rowSums(outs_burn_in[[i]][,c('I_01_l', 'I_01_f')]),
                  I_11 = rowSums(outs_burn_in[[i]][,c('I_11_l', 'I_11_f')])
  )
  
  
  t_adj_f<- data.frame(time = runs[[i]]$adj_f[,'time'],
                  C_00 = rowSums(runs[[i]]$adj_f[,c('C_00_u', 'C_00_f')]),
                  C_10 = rowSums(runs[[i]]$adj_f[,c('C_10_u', 'C_10_f')]),
                  C_01 = rowSums(runs[[i]]$adj_f[,c('C_01_u', 'C_01_f')]),
                  C_11 = rowSums(runs[[i]]$adj_f[,c('C_11_u', 'C_11_f')]),
                  I_00 = rowSums(runs[[i]]$adj_f[,c('I_00_l', 'I_00_f')]),
                  I_10 = rowSums(runs[[i]]$adj_f[,c('I_10_l', 'I_10_f')]),
                  I_01 = rowSums(runs[[i]]$adj_f[,c('I_01_l', 'I_01_f')]),
                  I_11 = rowSums(runs[[i]]$adj_f[,c('I_11_l', 'I_11_f')])
  )
  
  
  t_adj_l<- data.frame(time = runs[[i]]$adj_l[,'time'],
                       C_00 = rowSums(runs[[i]]$adj_l[,c('C_00_u', 'C_00_f')]),
                       C_10 = rowSums(runs[[i]]$adj_l[,c('C_10_u', 'C_10_f')]),
                       C_01 = rowSums(runs[[i]]$adj_l[,c('C_01_u', 'C_01_f')]),
                       C_11 = rowSums(runs[[i]]$adj_l[,c('C_11_u', 'C_11_f')]),
                       I_00 = rowSums(runs[[i]]$adj_l[,c('I_00_l', 'I_00_f')]),
                       I_10 = rowSums(runs[[i]]$adj_l[,c('I_10_l', 'I_10_f')]),
                       I_01 = rowSums(runs[[i]]$adj_l[,c('I_01_l', 'I_01_f')]),
                       I_11 = rowSums(runs[[i]]$adj_l[,c('I_11_l', 'I_11_f')])
  )
  
  
  t_diag_f<- data.frame(time = runs[[i]]$diag_f[,'time'],
                       C_00 = rowSums(runs[[i]]$diag_f[,c('C_00_u', 'C_00_f')]),
                       C_10 = rowSums(runs[[i]]$diag_f[,c('C_10_u', 'C_10_f')]),
                       C_01 = rowSums(runs[[i]]$diag_f[,c('C_01_u', 'C_01_f')]),
                       C_11 = rowSums(runs[[i]]$diag_f[,c('C_11_u', 'C_11_f')]),
                       I_00 = rowSums(runs[[i]]$diag_f[,c('I_00_l', 'I_00_f')]),
                       I_10 = rowSums(runs[[i]]$diag_f[,c('I_10_l', 'I_10_f')]),
                       I_01 = rowSums(runs[[i]]$diag_f[,c('I_01_l', 'I_01_f')]),
                       I_11 = rowSums(runs[[i]]$diag_f[,c('I_11_l', 'I_11_f')])
  )
  
  
  t_diag_l<- data.frame(time = runs[[i]]$diag_l[,'time'],
                       C_00 = rowSums(runs[[i]]$diag_l[,c('C_00_u', 'C_00_f')]),
                       C_10 = rowSums(runs[[i]]$diag_l[,c('C_10_u', 'C_10_f')]),
                       C_01 = rowSums(runs[[i]]$diag_l[,c('C_01_u', 'C_01_f')]),
                       C_11 = rowSums(runs[[i]]$diag_l[,c('C_11_u', 'C_11_f')]),
                       I_00 = rowSums(runs[[i]]$diag_l[,c('I_00_l', 'I_00_f')]),
                       I_10 = rowSums(runs[[i]]$diag_l[,c('I_10_l', 'I_10_f')]),
                       I_01 = rowSums(runs[[i]]$diag_l[,c('I_01_l', 'I_01_f')]),
                       I_11 = rowSums(runs[[i]]$diag_l[,c('I_11_l', 'I_11_f')])
  )
  
  t_new_f<- data.frame(time = runs[[i]]$new_f[,'time'],
                        C_000 = rowSums(runs[[i]]$new_f[,c('C_000_u', 'C_000_f')]),
                        C_100 = rowSums(runs[[i]]$new_f[,c('C_100_u', 'C_100_f')]),
                        C_010 = rowSums(runs[[i]]$new_f[,c('C_010_u', 'C_010_f')]),
                        C_001 = rowSums(runs[[i]]$new_f[,c('C_001_u', 'C_001_f')]),
                        C_110 = rowSums(runs[[i]]$new_f[,c('C_110_u', 'C_110_f')]),
                        C_011 = rowSums(runs[[i]]$new_f[,c('C_011_u', 'C_011_f')]),
                        C_101 = rowSums(runs[[i]]$new_f[,c('C_101_u', 'C_101_f')]),
                        C_111 = rowSums(runs[[i]]$new_f[,c('C_111_u', 'C_111_f')]),
                       I_000 = rowSums(runs[[i]]$new_f[,c('I_000_l', 'I_000_f', 'I_000_m')]),
                       I_100 = rowSums(runs[[i]]$new_f[,c('I_100_l', 'I_100_f', 'I_100_m')]),
                       I_010 = rowSums(runs[[i]]$new_f[,c('I_010_l', 'I_010_f', 'I_010_m')]),
                       I_001 = rowSums(runs[[i]]$new_f[,c('I_001_l', 'I_001_f', 'I_001_m')]),
                       I_110 = rowSums(runs[[i]]$new_f[,c('I_110_l', 'I_110_f', 'I_110_m')]),
                       I_011 = rowSums(runs[[i]]$new_f[,c('I_011_l', 'I_011_f', 'I_011_m')]),
                       I_101 = rowSums(runs[[i]]$new_f[,c('I_101_l', 'I_101_f', 'I_101_m')]),
                       I_111 = rowSums(runs[[i]]$new_f[,c('I_111_l', 'I_111_f', 'I_111_m')])
                       
  )
  
  t_new_l<- data.frame(time = runs[[i]]$new_l[,'time'],
                       C_000 = rowSums(runs[[i]]$new_l[,c('C_000_u', 'C_000_f')]),
                       C_100 = rowSums(runs[[i]]$new_l[,c('C_100_u', 'C_100_f')]),
                       C_010 = rowSums(runs[[i]]$new_l[,c('C_010_u', 'C_010_f')]),
                       C_001 = rowSums(runs[[i]]$new_l[,c('C_001_u', 'C_001_f')]),
                       C_110 = rowSums(runs[[i]]$new_l[,c('C_110_u', 'C_110_f')]),
                       C_011 = rowSums(runs[[i]]$new_l[,c('C_011_u', 'C_011_f')]),
                       C_101 = rowSums(runs[[i]]$new_l[,c('C_101_u', 'C_101_f')]),
                       C_111 = rowSums(runs[[i]]$new_l[,c('C_111_u', 'C_111_f')]),
                       I_000 = rowSums(runs[[i]]$new_l[,c('I_000_l', 'I_000_f', 'I_000_m')]),
                       I_100 = rowSums(runs[[i]]$new_l[,c('I_100_l', 'I_100_f', 'I_100_m')]),
                       I_010 = rowSums(runs[[i]]$new_l[,c('I_010_l', 'I_010_f', 'I_010_m')]),
                       I_001 = rowSums(runs[[i]]$new_l[,c('I_001_l', 'I_001_f', 'I_001_m')]),
                       I_110 = rowSums(runs[[i]]$new_l[,c('I_110_l', 'I_110_f', 'I_110_m')]),
                       I_011 = rowSums(runs[[i]]$new_l[,c('I_011_l', 'I_011_f', 'I_011_m')]),
                       I_101 = rowSums(runs[[i]]$new_l[,c('I_101_l', 'I_101_f', 'I_101_m')]),
                       I_111 = rowSums(runs[[i]]$new_l[,c('I_111_l', 'I_111_f', 'I_111_m')])
                       
  ) 
  
t_adj_f<- rbind(t0, t_adj_f)
t_adj_f$time<- seq(1:nrow(t_adj_f))

t_adj_l<- rbind(t0, t_adj_l)
t_adj_l$time<- seq(1:nrow(t_adj_l))

t_diag_f<- rbind(t0, t_diag_f)
t_diag_f$time<- seq(1:nrow(t_diag_f))

t_diag_l<- rbind(t0, t_diag_l)
t_diag_l$time<- seq(1:nrow(t_diag_l))


t0_new<- data.frame( time = t0$time, 
                     C_000 = t0$C_00,
            C_100 = t0$C_10,
            C_010 = rep(NA, nrow(t0)),
            C_001 = t0$C_01,
            C_110 = rep(NA, nrow(t0)),
            C_011 = rep(NA, nrow(t0)),
            C_101 = rep(NA, nrow(t0)),
            C_111 = t0$C_11,
            I_000 = t0$I_00,
            I_100 = t0$I_10,
            I_010 = rep(NA, nrow(t0)),
            I_001 = t0$I_01,
            I_110 = rep(NA, nrow(t0)),
            I_011 = rep(NA, nrow(t0)),
            I_101 = rep(NA, nrow(t0)),
            I_111 = t0$I_11)


t_new_f<- rbind(t0_new, t_new_f)
t_new_f$time<- seq(1:nrow(t_new_f))

t_new_l<- rbind(t0_new, t_new_l)
t_new_l$time<- seq(1:nrow(t_new_l))




m<- matrix(ncol = 12, nrow = 8, byrow = TRUE, data = c(
                                                        
                                                        0,1,1,0,5,5,0,9,9,0,13,13,
                                                        0,1,1,0,5,5,0,9,9,0,13,13,
                                                        
                                                        0,2,2,0,6,6,0,10,10,0,13,13,
                                                        0,2,2,0,6,6,0,10,10,0,13,13,
                                                        0,3,3,0,7,7,0,11,11,0,13,13,
                                                        0,3,3,0,7,7,0,11,11,0,13,13,
                                                        
                                                        0,4,4,0,8,8,0,12,12,0,13,13,
                                                        0,4,4,0,8,8,0,12,12,0,13,13))

  
layout(m)
par(mar = c(3, 0.5, 0.5, 0.5), oma = c(1,2.5,2.5,1))

#layout.show(13)


# ADJUVANT_FRONT_ILL
plot(t_adj_f$I_00~t_adj_f$time,  col = 'black', ylim = c(0, max(t_adj_f[,6:9])), type = 'l', lwd =3, ylab = '', xlab = '')
lines(t_adj_f$I_10~t_adj_f$time,  col = 'dodgerblue', lwd =3)
lines(t_adj_f$I_01~t_adj_f$time,  col = 'gold', lwd =3)
lines(t_adj_f$I_11~t_adj_f$time,  col = 'firebrick', lwd =3)
mtexti('Amount',   side = 2, font = 1, cex = 1.5, off = 0.4)
mtexti('ADJUVANT',   side = 3, font = 2, cex = 2, off = 0.2)
mtexti('SYMPTOMATIC CLASS',   side = 2, font = 1, cex = 1.5, off = 0.6)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)


# ADJUVANT_FRONT_CARRIER
plot(t_adj_f$C_00~t_adj_f$time,  col = 'black', ylim = c(0, max(t_adj_f[,2:5])), type = 'l', lwd =3, ylab = 'Amount', xlab = 'time')
lines(t_adj_f$C_10~t_adj_f$time,  col = 'dodgerblue', lwd =3)
lines(t_adj_f$C_01~t_adj_f$time,  col = 'gold', lwd =3)
lines(t_adj_f$C_11~t_adj_f$time,  col = 'firebrick', lwd =3)
mtexti('Amount',   side = 2, font = 1, cex = 1.5, off = 0.4)
mtexti('CARRIER CLASS',   side = 2, font = 1, cex = 1.5, off = 0.6)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)


# ADJUVANT_LAST_ILL
plot(t_adj_l$I_00~t_adj_l$time,  col = 'black', ylim = c(0, max(t_adj_l[,6:9])), type = 'l', lwd =3, ylab = 'Amount', xlab = '')
lines(t_adj_l$I_10~t_adj_l$time,  col = 'dodgerblue', lwd =3)
lines(t_adj_l$I_01~t_adj_l$time,  col = 'gold', lwd =3)
lines(t_adj_l$I_11~t_adj_l$time,  col = 'firebrick', lwd =3)
mtexti('Amount',   side = 2, font = 1, cex = 1.5, off = 0.4)
mtexti('SYMPTOMATIC CLASS',   side = 2, font = 1, cex = 1.5, off = 0.6)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)

# ADJUVANT_LAST_CARRIER
plot(t_adj_l$C_00~t_adj_l$time,  col = 'black', ylim = c(0, max(t_adj_l[,2:5])), type = 'l', lwd =3, ylab = 'Amount', xlab = 'time')
lines(t_adj_l$C_10~t_adj_l$time,  col = 'dodgerblue', lwd =3)
lines(t_adj_l$C_01~t_adj_l$time,  col = 'gold', lwd =3)
lines(t_adj_l$C_11~t_adj_l$time,  col = 'firebrick', lwd =3)
mtexti('Amount',   side = 2, font = 1, cex = 1.5, off = 0.4)
mtexti('Time',   side = 1, font = 1, cex = 1.5, off = 0.4)
mtexti('CARRIER CLASS',   side = 2, font = 1, cex = 1.5, off = 0.6)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)



# DIAGNOSTIC_FRONT_ILL
plot(t_diag_f$I_00~t_diag_f$time,  col = 'black', ylim = c(0, max(t_diag_f[,6:9])), type = 'l', lwd =3, ylab = '', xlab = '')
lines(t_diag_f$I_10~t_diag_f$time,  col = 'dodgerblue', lwd =3)
lines(t_diag_f$I_01~t_diag_f$time,  col = 'gold', lwd =3)
lines(t_diag_f$I_11~t_diag_f$time,  col = 'firebrick', lwd =3)
mtexti('DIAGNOSTIC',   side = 3, font = 2, cex = 2, off = 0.2)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)

# DIAGNOSTIC_FRONT_CARRIER
plot(t_diag_f$C_00~t_diag_f$time,  col = 'black', ylim = c(0, max(t_diag_f[,2:5])), type = 'l', lwd =3, ylab = '', xlab = 'time')
lines(t_diag_f$C_10~t_diag_f$time,  col = 'dodgerblue', lwd =3)
lines(t_diag_f$C_01~t_diag_f$time,  col = 'gold', lwd =3)
lines(t_diag_f$C_11~t_diag_f$time,  col = 'firebrick', lwd =3)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)


# DIAGNOSTIC_LAST_ILL
plot(t_diag_l$I_00~t_diag_l$time,  col = 'black', ylim = c(0, max(t_diag_l[,6:9])), type = 'l', lwd =3, ylab = '', xlab = '')
lines(t_diag_l$I_10~t_diag_l$time,  col = 'dodgerblue', lwd =3)
lines(t_diag_l$I_01~t_diag_l$time,  col = 'gold', lwd =3)
lines(t_diag_l$I_11~t_diag_l$time,  col = 'firebrick', lwd =3)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)

# DIAGNOSTIC_LAST_CARRIER
plot(t_diag_l$C_00~t_diag_l$time,  col = 'black', ylim = c(0, max(t_diag_l[,2:5])), type = 'l', lwd =3, ylab = '', xlab = 'time')
lines(t_diag_l$C_10~t_diag_l$time,  col = 'dodgerblue', lwd =3)
lines(t_diag_l$C_01~t_diag_l$time,  col = 'gold', lwd =3)
lines(t_diag_l$C_11~t_diag_l$time,  col = 'firebrick', lwd =3)
mtexti('Time',   side = 1, font = 1, cex = 1.5, off = 0.4)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)


# NEW_FRONT_ILL
plot(t_new_f$I_000~t_new_f$time,  col = 'black', ylim = c(0, max(c(max(t0[,6:9]),max(t_new_f[302:nrow(t_new_f),10:17])))), type = 'l', lwd =3, ylab = '', xlab = '')
lines(t_new_f$I_100~t_new_f$time,  col = 'dodgerblue', lwd =3)
lines(t_new_f$I_010~t_new_f$time,  col = 'darkorange', lwd =3)
lines(t_new_f$I_001~t_new_f$time,  col = 'gold', lwd =3)
lines(t_new_f$I_110~t_new_f$time,  col = 'tan1', lwd =3)
lines(t_new_f$I_011~t_new_f$time,  col = 'seagreen', lwd =3)
lines(t_new_f$I_101~t_new_f$time,  col = 'violet', lwd =3)
lines(t_new_f$I_111~t_new_f$time,  col = 'firebrick', lwd =3)
mtexti('NEW DRUG',   side = 3, font = 2, cex = 2, off = 0.2)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)

# NEW_FRONT_CARRIER
plot(t_new_f$C_000~t_new_f$time,  col = 'black', ylim = c(0, max(c(max(t0[,6:9]),max(t_new_f[302:nrow(t_new_f),2:9])))), type = 'l', lwd =3, ylab = '', xlab = 'time')
lines(t_new_f$C_100~t_new_f$time,  col = 'dodgerblue', lwd =3)
lines(t_new_f$C_010~t_new_f$time,  col = 'darkorange', lwd =3)
lines(t_new_f$C_001~t_new_f$time,  col = 'gold', lwd =3)
lines(t_new_f$C_110~t_new_f$time,  col = 'tan1', lwd =3)
lines(t_new_f$C_011~t_new_f$time,  col = 'seagreen', lwd =3)
lines(t_new_f$C_101~t_new_f$time,  col = 'violet', lwd =3)
lines(t_new_f$C_111~t_new_f$time,  col = 'firebrick', lwd =3)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)

# NEW_LAST_ILL
plot(t_new_l$I_000~t_new_l$time,  col = 'black', ylim = c(0, max(c(max(t0[,6:9]),max(t_new_l[302:nrow(t_new_l),10:17])))), type = 'l', lwd =3, ylab = '', xlab = '')
lines(t_new_l$I_100~t_new_l$time,  col = 'dodgerblue', lwd =3)
lines(t_new_l$I_010~t_new_l$time,  col = 'darkorange', lwd =3)
lines(t_new_l$I_001~t_new_l$time,  col = 'gold', lwd =3)
lines(t_new_l$I_110~t_new_l$time,  col = 'tan1', lwd =3)
lines(t_new_l$I_011~t_new_l$time,  col = 'seagreen', lwd =3)
lines(t_new_l$I_101~t_new_l$time,  col = 'violet', lwd =3)
lines(t_new_l$I_111~t_new_l$time,  col = 'firebrick', lwd =3)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)

# NEW_LAST_CARRIER
plot(t_new_l$C_000~t_new_l$time,  col = 'black', ylim = c(0, max(c(max(t0[,6:9]),max(t_new_l[302:nrow(t_new_l),2:9])))), type = 'l', lwd =3, ylab = '', xlab = 'time')
lines(t_new_l$C_100~t_new_l$time,  col = 'dodgerblue', lwd =3)
lines(t_new_l$C_010~t_new_l$time,  col = 'darkorange', lwd =3)
lines(t_new_l$C_001~t_new_l$time,  col = 'gold', lwd =3)
lines(t_new_l$C_110~t_new_l$time,  col = 'tan1', lwd =3)
lines(t_new_l$C_011~t_new_l$time,  col = 'seagreen', lwd =3)
lines(t_new_l$C_101~t_new_l$time,  col = 'violet', lwd =3)
lines(t_new_l$C_111~t_new_l$time,  col = 'firebrick', lwd =3)
mtexti('Time',   side = 1, font = 1, cex = 1.5, off = 0.4)
abline(v = 300, col = 'darkgrey', lty = 2, lwd = 1.5)

mtext('                                                                    FRONT LINE INTERVENTION', side = 2, outer = TRUE, font = 2)
mtext('LAST LINE INTERVENTION                                                                    ', side = 2, outer = TRUE, font = 2)

# legend ----

plot(0, xlim= c(0,10), ylim= c(0,10), xaxt = 'n', yaxt = 'n', bty = 'n', pch = '', ylab = '', xlab = '')

text(expression(beta['w']),x = 0.5, y = 10, cex = 2)
text(paste('=', round(pars[i,'bw'], 2)),x = 4.5, y = 10, cex = 2)

text(expression(beta['b']),x=0.5, y = 9.5, cex = 2)
text(paste('=', round(pars[i,'bb'], 2)),x = 4.5, y = 9.5, cex = 2)


#text(paste('c =', round(pars[i,'c'], 2)), x = 5, y = 8.0,  cex = 2)

text(expression(c),x = 0.5, y = 9.0, cex = 2)
text(paste('=', round(pars[i,'c'], 2)),x = 4.5, y = 9.0, cex = 2)

text(expression(delta),x = 0.5, y = 8.5, cex = 2)
text(paste('=', round(pars[i,'d'], 2)),x = 4.5, y = 8.5, cex = 2)

text(expression(gamma),x = 0.5, y = 8.0, cex = 2)
text(paste('=', round(pars[i,'g'], 2)),x = 4.5, y = 8.0, cex = 2)

text(expression(zeta),x = 0.5, y = 7.5, cex = 2)
text(paste('=', round(pars[i,'z'], 2)),x = 4.5, y = 7.5, cex = 2)

text(expression(alpha),x = 0.5, y = 7.0, cex = 2)
text(paste('=', round(pars[i,'a'], 2)),x = 4.5, y = 7.0, cex = 2)

text(expression(epsilon),x = 0.5, y = 6.5, cex = 2)
text(paste('=', round(pars[i,'e'], 2)),x = 4.5, y = 6.5, cex = 2)


#text(paste('p =', round(pars[i,'p'], 2)), x = 5, y = 5, cex = 2)

text(expression(p), x = 0.5, y = 6.0, cex = 2)
text(paste('=', round(pars[i,'p'], 2)),x = 4.5, y = 6.0, cex = 2)


text(expression(mu), x = 0.5, y = 5.5, cex = 2)
text(paste('=', round(pars[i,'mu'], 2)), x = 4, y = 5.5, cex = 2)

text(expression(f), x = 0.5, y = 5, cex = 2)
text(paste('=', round(pars[i,'f'], 2)), x = 5, y = 5, cex = 2)


#text(paste('f =', round(pars[i,'f'], 2)), x = 5, y = 4, cex = 2)
text(paste('seed =', pars[i,'my_seed']), x = 3, y = 4.5, cex = 2)
text(paste('draw =', round(pars[i,'draw'], 0)), x = 3.5, y = 4, cex = 2)


lines(x = c(0,1,2), y = c(3.5,3.5,3.5), lwd = 4, col = 'black')
text('00 or 000',x = 6, y = 3.5, cex = 2)

lines(x = c(0,1,2), y = c(3,3,3), lwd = 4, col = 'dodgerblue')
text('10 or 100',x = 6, y = 3, cex = 2)

lines(x = c(0,1,2), y = c(2.5,2.5,2.5), lwd = 4, col = 'gold')
text('01 or 001',x = 6, y = 2.5, cex = 2)

lines(x = c(0,1,2), y = c(2,2,2), lwd = 4, col = 'firebrick')
text('11 or 111',x = 6, y = 2, cex = 2)


lines(x = c(0,1,2), y = c(1.5,1.5,1.5), lwd = 4, col = 'darkorange')
text('010', x = 4, y = 1.5, cex = 2)

lines(x = c(0,1,2), y = c(1,1,1), lwd = 4, col = 'tan1')
text('110', x = 4, y = 1, cex = 2)

lines(x = c(0,1,2), y = c(0.5,0.5,0.5), lwd = 4, col = 'seagreen')
text('011', x = 4, y = 0.5, cex = 2)

lines(x = c(0,1,2), y = c(0,0,0), lwd = 4, col = 'violet')
text('101', x = 4, y = 0, cex = 2)

}

}

  print(j)
}


dev.off()











