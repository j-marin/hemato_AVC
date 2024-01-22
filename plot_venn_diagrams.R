# 15/01/2024
# antoine.bridier-nahmias@inserm.fr

library(tidyverse)
library(eulerr)

df <- 
  readxl::read_xlsx('SupplementaryTables_Hemato_decembre2023.xlsx', sheet = 5)
names(df) <- df[2,]

df <- 
  df[-(1:2),]



############ Patient A
df |> 
  filter(Patient == 'A') |> 
  select(Site, 'SNP ID') |>
  distinct() |> 
  as.data.frame() -> patient_A

blood_A <- patient_A[patient_A$Site == 'blood', 'SNP ID']
stool_A <- patient_A[patient_A$Site == 'stool', 'SNP ID']
urine_A <- patient_A[patient_A$Site == 'urine', 'SNP ID']

blood_stool_urine_A <- length(intersect(intersect(blood_A, stool_A), urine_A))
blood_stool_A <- length(intersect(blood_A, stool_A)) - blood_stool_urine_A
blood_urine_A <- length(intersect(blood_A, urine_A)) - blood_stool_urine_A
stool_urine_A <- length(intersect(stool_A, urine_A)) - blood_stool_urine_A

blood_A <- length(blood_A) - blood_stool_urine_A - blood_stool_A - blood_urine_A
stool_A <- length(stool_A) - blood_stool_urine_A - blood_stool_A - stool_urine_A
urine_A <- length(urine_A) - blood_stool_urine_A - blood_urine_A - stool_urine_A

eul_A <-
  euler(combinations = c('blood' = blood_A, 
                         'stool' = stool_A,
                         'urine' = urine_A,
                         'blood&stool' = blood_stool_A,
                         'blood&urine' = blood_urine_A,
                         'stool&urine' = stool_urine_A,
                         'blood&stool&urine' = blood_stool_urine_A), 
        shape = 'ellipse')

plot(eul_A, quantities = TRUE, 
     fills = list(fill = c('firebrick1', 'blue', 'gold1'), 
                  alpha = c(0.5, 0.5, 0.5)), main = 'Patient A')




############ Patient C
df |> 
  filter(Patient == 'C') |> 
  select(Site, 'SNP ID') |>
  distinct() |> 
  as.data.frame() -> patient_C

blood_C <- patient_C[patient_C$Site == 'blood', 'SNP ID']
stool_C <- patient_C[patient_C$Site == 'stool', 'SNP ID']
urine_C <- patient_C[patient_C$Site == 'urine', 'SNP ID']

blood_stool_urine_C <- length(intersect(intersect(blood_C, stool_C), urine_C))
blood_stool_C <- length(intersect(blood_C, stool_C)) - blood_stool_urine_C
blood_urine_C <- length(intersect(blood_C, urine_C)) - blood_stool_urine_C
stool_urine_C <- length(intersect(stool_C, urine_C)) - blood_stool_urine_C

blood_C <- length(blood_C) - blood_stool_urine_C - blood_stool_C - blood_urine_C
stool_C <- length(stool_C) - blood_stool_urine_C - blood_stool_C - stool_urine_C
urine_C <- length(urine_C) - blood_stool_urine_C - blood_urine_C - stool_urine_C

eul_C <-
  euler(combinations = c('blood' = blood_C, 
                         'stool' = stool_C,
                         'urine' = urine_C,
                         'blood&stool' = blood_stool_C,
                         'blood&urine' = blood_urine_C,
                         'stool&urine' = stool_urine_C,
                         'blood&stool&urine' = blood_stool_urine_C), 
        shape = 'ellipse')

plot(eul_C, quantities = TRUE, 
     fills = list(fill = c('firebrick1', 'blue', 'gold1'), 
                  alpha = c(0.5, 0.5, 0.5)), main = 'Patient C')


############ Patient B

df |> 
  filter(Patient == 'B') |> 
  select(Site, 'SNP ID') |>
  distinct() |> 
  as.data.frame() -> patient_B

blood_B <- patient_B[patient_B$Site == 'blood', 'SNP ID']
stool_B <- patient_B[patient_B$Site == 'stool', 'SNP ID']
urine_B <- patient_B[patient_B$Site == 'urine', 'SNP ID']

blood_stool_urine_B <- length(intersect(intersect(blood_B, stool_B), urine_B))
blood_stool_B <- length(intersect(blood_B, stool_B)) - blood_stool_urine_B
blood_urine_B <- length(intersect(blood_B, urine_B)) - blood_stool_urine_B
stool_urine_B <- length(intersect(stool_B, urine_B)) - blood_stool_urine_B

blood_B <- length(blood_B) - blood_stool_urine_B - blood_stool_B - blood_urine_B
stool_B <- length(stool_B) - blood_stool_urine_B - blood_stool_B - stool_urine_B
urine_B <- length(urine_B) - blood_stool_urine_B - blood_urine_B - stool_urine_B

eul_B <-
  euler(combinations = c('blood' = blood_B, 
                         'stool' = stool_B,
                         'blood&stool' = blood_stool_B),
        shape = 'ellipse')

plot(eul_B, quantities = TRUE, 
     fills = list(fill = c('firebrick1', 'blue'), 
                  alpha = c(0.5, 0.5)), main = 'Patient B')



############ Plot Saving
dev.off()
pdf(file = 'Patient_A.pdf')
plot(eul_A, quantities = TRUE, 
     fills = list(fill = c('firebrick1', 'blue', 'gold1'), 
                  alpha = c(0.5, 0.5, 0.5)), main = 'Patient A', 
     labels = list(fontsize = 15))
dev.off()

pdf(file = 'Patient_C.pdf')
plot(eul_C, quantities = TRUE, 
     fills = list(fill = c('firebrick1', 'blue', 'gold1'), 
                  alpha = c(0.5, 0.5, 0.5)), main = 'Patient C',
     labels = list(fontsize = 15))
dev.off()

pdf(file = 'Patient_B.pdf')
plot(eul_B, quantities = TRUE, 
     fills = list(fill = c('firebrick1', 'blue'), 
                  alpha = c(0.5, 0.5)), main = 'Patient B',
     labels = list(fontsize = 15))
dev.off()
