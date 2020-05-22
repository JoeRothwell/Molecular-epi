library(survival)
base <- Cncr_Caco_Clrt ~ Qe_Energy + L_School + Smoke_Int + Smoke_Stat + #Smoke_Intensity + 
  Height_C + strata(Match_Caseset)
# Dairy product intake = QGE05

# Biocrates small, large, fatty acids (large fasted subset did not improve OR)
# CRC A: score and signature
library(survival)
fit1 <- clogit(update(base, ~. + comp1), data = crc1.ph)
fit2 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc1.ph)

fit1dp <- clogit(update(base, ~. + comp1 + Qge05), data = crc1.ph)
fit2dp <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = crc1.ph)

# CRC B Score, signature and by sex
fit3 <- clogit(update(base, ~. + comp1), data = crc2.ph)
fit4 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc2.ph)

fit3dp <- clogit(update(base, ~. + comp1 + Qge05), data = crc2.ph)
fit4dp <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = crc2.ph)


# CRC A Fatty acids score, signature and by sex
fit5 <- clogit(update(base, ~. + comp2), data = crc3.ph)
fit6 <- clogit(update(base, ~. + Wcrf_C_Cal), data = crc3.ph)

fit5dp <- clogit(update(base, ~. + comp2 + Qge05), data = crc3.ph)
fit6dp <- clogit(update(base, ~. + Wcrf_C_Cal + Qge05), data = crc3.ph)
