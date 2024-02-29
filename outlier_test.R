library(tidyverse)
library(outliers)

data <- read_csv(file = "/Users/emilyfitzmeyer/Desktop/outlier_test.csv")
data <- slice(data, 1:16)

boxplot(data$data, horizontal = TRUE)
data_v <- data$data_1

summary(data_v)
iqr <- IQR(data_v)

tmin <- 0.27292-(3*iqr)
tmax <- 0.91007+(3*iqr)

data_v[which(data_v < tmin | data_v > tmax)]

chisq.out.test(data_v)

