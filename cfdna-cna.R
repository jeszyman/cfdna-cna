

library(tidyverse)

load("/mnt/ris/aadel/mpnst/data_model/data_model.RData")

ls()

names(libraries_full)

class(libraries_full$collect_date)

libraries_full$collect_date = as.Date(libraries_full$collect_date)

as.numeric(libraries_full$collect_date[[1]]- libraries_full$collect_date[[2]])

test =
  libraries_full %>% arrange(collect_date) %>% group_by(participant_id, isolation_type) %>%
  mutate(collect_day = as.numeric(collect_date - first(collect_date))) %>%
  mutate(collect_day = replace_na(collect_day, 0))

tf = read.table("/tmp/tf.tsv", header = F, sep = '\t')
colnames(tf) = c("libnfrag", "tf", "ploidy")
tf$library_id = substr(tf$libnfrag, 1, 6)

tf2 =
  tf %>% filter(grepl("sub20m_frag90", libnfrag))


test2=tf2 %>% left_join(test, by = "library_id")

write.csv(file ="/tmp/test.csv", test2)
test %>% select(participant_id, collect_day) %>% arrange(participant_id) %>% print(n = Inf)



test$collect_day

  case_when(collect_date == first(collect_date) ~ 0,
                                 collect_date > first(collect_date) ~ collect_date - first(collect_date)))




) %>% select(library_id, participant_id, collect_day)
