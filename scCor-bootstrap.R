
count_files <- list.files("Data", pattern = "count_norm", full = TRUE)
count_files[1] |> basename()

count_12dpi <- readRDS(count_files[1])

class(count_12dpi)
dim(count_12dpi)

count_4dpi <- readRDS(count_files[2])

class(count_4dpi)
dim(count_4dpi)

library(tidyverse)
# rearrange matrix to list gene of interest in first column

temp_12dpi <-
count_12dpi |> as.data.frame() |>
  select(`nbis-gene-2-utr`, everything()) |>
  as.matrix()

temp_4dpi <-
count_4dpi |> as.data.frame() |>
  select(`nbis-gene-2-utr`, everything()) |>
  as.matrix()


library(scLink)
# use sc_link to get cor
# cannot run ncores > 1 on Windows

# check to make sure that the output has the correct column

temp_cor_12dpi <- sclink_cor(temp_12dpi, ncores = 1)
temp_cor_12dpi[, 1]

temp_cor_4dpi <- sclink_cor(temp_4dpi, ncores = 1)
temp_cor_4dpi[, 1]


# function to run sclink_cor and extract the first column
# I would increase the number of cores, but it's not allowed in Windows

get_cor <- function(data, indices) {
  d <- data[indices, ]
  x <- sclink_cor(expr = d, ncores = 1)
  return(x[, 1])
}


# this was my calculation to see how long this might take on the system

set.seed(52511)
system.time(
try_it <- get_cor(temp_12dpi, indices = 1:nrow(temp))
)

# I set these up to run overnight on my dinky computer
18.67 * 1000 / 60 / 60 # 5 hours of time for 1000 copies


# set start and end times just as a matter of curiosity

start_12 <- Sys.time()
set.seed(52511)
cor_boot_12 <- boot::boot(temp_12dpi, get_cor, R = 1000)
end_12 <- Sys.time()
save.image() 
# I save the results in case something goes wrong between the time the
#  code stops running and I return to access my computer


end_12 - start12

# Time difference of 4.794911 hours (5 hours to run)


start_4 <- Sys.time()
set.seed(52511)
cor_boot_4 <- boot::boot(temp_4dpi, get_cor, R = 1000)
end_4 <- Sys.time()
save.image()

end_4 - start4
# Time difference of 7.570419 hours (almost 8 hours for the bigger file)

# look at the structure of the output
glimpse(cor_boot_12) 


# calculate how many are positive or negative
p12_positive <- apply(cor_boot_12$t, 2, function(x) sum(x > 0, na.rm = TRUE))
p12_negative <- apply(cor_boot_12$t, 2, function(x) sum(x < 0, na.rm = TRUE))


p12 <- tibble(
  gene = names(cor_boot_12$t0),
  cor = cor_boot_12$t0,
  greater_than_zero = p12_positive,
  less_than_zero = p12_negative,
# there were warnings because of zero standard deviation in some correlations
# this happens when all of the values are the same (in this case, all are zero)
  all_zero = apply(cor_boot_12$t, 2, function(x) sum(is.na(x)))) |>
# divide by 1000 because there are 1000 bootstrap runs
  mutate(p_value = ifelse(cor > 0, less_than_zero, greater_than_zero),
         p_value = p_value / 1000)

# adjust for false discovery rate
p12$q <- p.adjust(p12$p_value, "fdr")


# I included the code inside the tibble for the 4dpi data

p4 <- tibble(
  gene = names(cor_boot_4$t0),
  cor = cor_boot_4$t0,
  greater_than_zero = apply(cor_boot_4$t, 2, function(x) sum(x > 0, na.rm = TRUE)),
  less_than_zero = apply(cor_boot_4$t, 2, function(x) sum(x < 0, na.rm = TRUE)),
  all_zero = apply(cor_boot_4$t, 2, function(x) sum(is.na(x)))) |>
  mutate(p_value = ifelse(cor > 0, less_than_zero, greater_than_zero),
         p_value = p_value / 1000)

# adjust for false discovery rate
p4$q <- p.adjust(p4$p_value, "fdr")



# visualize distribution of correlations


test <-
cor_boot_4$t |> as_tibble() 
names(test) <- names(cor_boot_4$t0)

mypath <- getwd()

library(trelliscopejs)

formatP <- function(x) {
  ifelse(x < 0.0001, "<0.0001",
    ifelse(x < 0.001, sprintf("%0.4f", x), 
      ifelse(x < 0.01, sprintf("%0.3f", x), sprintf("%0.2f", x))))
}

colors <- viridis::viridis(3)
names(colors) <-  c("c", "b", "a")

test |> 
pivot_longer(everything(),names_to = "gene", values_to = "cor") |>
  replace_na(list(cor = 0)) |>
  mutate(zeros = case_when(
    cor < 0 ~ "a",
    cor == 0 ~ "b",
    cor > 0 ~ "c")) |>
left_join(p24 |> select(gene, q), by = "gene") |>
ggplot(aes(x = cor)) +
  geom_histogram(binwidth = 0.05, aes(fill = zeros), color = "darkgrey") +
  theme_bw() +
  guides(fill = "none", color = "none") +
  labs(x = "Bootstrap correlations") +
  scale_fill_manual(values = colors) +
  geom_vline(xintercept = 0, color = "firebrick", linetype = "dotted") +
  geom_text(aes(label = paste0("p = ", formatP(q))), y = 100, x = -0.5)  +
  facet_trelliscope(~ gene,
    name = "4dpi Correlations",
    path = emily)


test2 <- cor_boot_12$t |> as_tibble() 
names(test2) <- names(cor_boot_12$t0)

test2 |> 
pivot_longer(everything(),names_to = "gene", values_to = "cor") |>
  replace_na(list(cor = 0)) |>
  mutate(zeros = case_when(
    cor < 0 ~ "a",
    cor == 0 ~ "b",
    cor > 0 ~ "c")) |>
left_join(p12 |> select(gene, q), by = "gene") |>
ggplot(aes(x = cor)) +
  geom_histogram(binwidth = 0.05, aes(fill = zeros), color = "darkgrey") +
  theme_bw() +
  guides(fill = "none", color = "none") +
  labs(x = "Bootstrap correlations") +
  scale_fill_manual(values = colors) +
  geom_vline(xintercept = 0, color = "firebrick", linetype = "dotted") +
  geom_text(aes(label = paste0("p = ", formatP(q))), y = 100, x = -0.5)  +
  facet_trelliscope(~ gene,
    name = "12dpi Correlations",
    path = emily)



################# save stuff for Emily

emily <- "Results for Fitzmeyer"

saveRDS(cor_boot_12, file = paste0(emily, "/Full_bootstrap_12dpi.rds"))
saveRDS(cor_boot_4, file = paste0(emily, "/Full_bootstrap_4dpi.rds"))
write_csv(p12, file = paste0(emily, "/P-value_calculations_12dpi.csv"))
write_csv(p4, file = paste0(emily, "/P-value_calculations_4dpi.csv"))
