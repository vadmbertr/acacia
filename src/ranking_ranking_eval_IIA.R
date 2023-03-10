library(ggplot2)

step3 = df_median_centerscale_transfo %>%
  group_by(candidate,cat_score,name_score) %>%
  summarise(score3=mean(trendval, na.rm=T))
step4 = step3 %>%
  group_by(candidate,cat_score) %>%
  summarise(score4=geomMean(score3))
step5 = step4 %>%
  group_by(candidate) %>%
  summarise(score5=geomMean(score4))

step4b = df_median_centerscale_transfo %>%
  group_by(candidate,cat_score,dataset) %>%
  summarise(score4b=mean(trendval, na.rm=T))
step5b = step4b %>%
  group_by(candidate,dataset) %>%
  summarise(score5b=geomMean(score4b))
step3b = step5b %>%
  group_by(candidate) %>%
  summarise(score3b=geomMean(score5b))

#----
# Compare ranking in 12345 vs 12453
#----
plot(step3b$score3b, step5$score5)
df = data.frame(score=c(step5$score5, step3b$score3b),
                ranking=c(rep(c("345","453"),each=nrow(step5))),
                candidate=c(step5$candidate, step3b$candidate))
ggplot(df, aes(x=ranking, y=score, color=candidate, group=candidate)) +
  geom_point() +
  geom_line(linewidth=2, alpha=.5) +
  ggrepel::geom_label_repel(data=df[df$ranking=="345",], aes(label=candidate)) +
  theme_abyss()
ggsave("ranking_ranking_eval_IIA/comp_ranking_order.pdf",
       width = 6, height = 4)

#----
# Eval step 3
#----
compa_step3_bef = df_median_centerscale_transfo[,c("name_score","cat_score","candidate","dataset","trendval")] %>%
  rename(score=trendval)
compa_step3_aft = step3[,c("name_score","cat_score","candidate","score3")] %>%
  rename(score=score3)
compa_step3_aft$dataset = "Aggregated"

compa_step3 <-  bind_rows(compa_step3_bef,compa_step3_aft)
compa_step3$datatype <- sapply(compa_step3$dataset, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step3$datatype <- factor(compa_step3$datatype, levels=c("Single","Aggregated"))

for (i in unique(compa_step3$name_score)) {
print(ggplot(compa_step3[compa_step3$name_score==i,], aes(x=datatype, y=score, color=dataset)) +
    geom_point() +
    facet_wrap(~candidate) +
    theme_modern())
ggsave(paste0("ranking_ranking_eval_IIA/step3_",i,".pdf"),
       width = 10, height = 8)
}

compa_step3_bef_order = compa_step3_bef %>%
  group_by(name_score,dataset) %>%
  mutate(rank=order(score)) %>%
  select(name_score,candidate,dataset,rank)
compa_step3_aft_order = compa_step3_aft %>%
  group_by(name_score) %>%
  mutate(rank=order(score)) %>%
  select(name_score,candidate,rank)
compa_step3_aft_order$dataset = "Aggregated"

compa_step3_order <-  bind_rows(compa_step3_bef_order,compa_step3_aft_order)
compa_step3_order$datatype <- sapply(compa_step3_order$dataset, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step3_order$datatype <- factor(compa_step3_order$datatype, levels=c("Single","Aggregated"))

order_change = sapply(unique(compa_step3_order$name_score), function(x)
  sapply(unique(compa_step3_order$candidate), function(y) {
    df = compa_step3_order[compa_step3_order$name_score==x & compa_step3_order$candidate==y,]
    df_bef = df$rank[df$datatype=="Single"]
    df_aft = df$rank[df$datatype=="Aggregated"]
    df_aft <= max(df_bef) & df_aft >= min(df_bef)
  }))
mean(order_change)
compa_step3_order$OK = NA
for (x in unique(compa_step3_order$name_score)) {
  for (y in unique(compa_step3_order$candidate)) {
    compa_step3_order$OK[compa_step3_order$candidate==y & compa_step3_order$name_score==x] = order_change[y,x]
  }
}
for (i in unique(compa_step3_order$name_score)) {
  print(ggplot(compa_step3_order[compa_step3_order$name_score==i,], aes(x=datatype, y=rank, color=dataset)) +
          geom_point() +
          facet_wrap(~candidate) +
          title(i) +
          ggrepel::geom_label_repel(data=compa_step3_order[compa_step3_order$name_score==i & compa_step3_order$datatype=="Aggregated",], aes(label=OK)) +
          theme_modern())
  ggsave(paste0("ranking_ranking_eval_IIA/step3_order_",i,".pdf"),
         width = 10, height = 8)
}

#----
# Eval step 4
#----
compa_step4_bef = step3[,c("name_score","cat_score","candidate","score3")] %>%
  rename(score=score3)
compa_step4_aft = step4[,c("cat_score","candidate","score4")] %>%
  rename(score=score4)
compa_step4_aft$name_score = "Aggregated"

compa_step4 <-  bind_rows(compa_step4_bef,compa_step4_aft)
compa_step4$scoretype <- sapply(compa_step4$name_score, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step4$scoretype <- factor(compa_step4$scoretype, levels=c("Single","Aggregated"))

for (i in unique(compa_step4$cat_score)) {
  print(ggplot(compa_step4[compa_step4$cat_score==i,], aes(x=scoretype, y=score, color=name_score)) +
          geom_point() +
          facet_wrap(~candidate) +
          theme_modern())
  ggsave(paste0("ranking_ranking_eval_IIA/step4_",i,".pdf"),
         width = 10, height = 8)
}

compa_step4_bef_order = compa_step4_bef %>%
  group_by(name_score,cat_score) %>%
  mutate(rank=order(score)) %>%
  select(name_score,cat_score,candidate,rank)
compa_step4_aft_order = compa_step4_aft %>%
  group_by(cat_score) %>%
  mutate(rank=order(score)) %>%
  select(cat_score,candidate,rank)
compa_step4_aft_order$name_score = "Aggregated"

compa_step4_order <-  bind_rows(compa_step4_bef_order,compa_step4_aft_order)
compa_step4_order$scoretype <- sapply(compa_step4_order$name_score, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step4_order$scoretype <- factor(compa_step4_order$scoretype, levels=c("Single","Aggregated"))

order_change = sapply(unique(compa_step4_order$cat_score), function(x)
  sapply(unique(compa_step4_order$candidate), function(y) {
    df = compa_step4_order[compa_step4_order$cat_score==x & compa_step4_order$candidate==y,]
    df_bef = df$rank[df$scoretype=="Single"]
    df_aft = df$rank[df$scoretype=="Aggregated"]
    df_aft <= max(df_bef) & df_aft >= min(df_bef)
  }))
mean(order_change)
compa_step4_order$OK = NA
for (x in unique(compa_step4_order$cat_score)) {
  for (y in unique(compa_step4_order$candidate)) {
    compa_step4_order$OK[compa_step4_order$candidate==y & compa_step4_order$cat_score==x] = order_change[y,x]
  }
}
for (i in unique(compa_step4_order$cat_score)) {
  print(ggplot(compa_step4_order[compa_step4_order$cat_score==i,], aes(x=scoretype, y=rank, color=name_score)) +
          geom_point() +
          facet_wrap(~candidate) +
          title(i) +
          ggrepel::geom_label_repel(data=compa_step4_order[compa_step4_order$cat_score==i & compa_step4_order$scoretype=="Aggregated",], aes(label=OK)) +
          theme_modern())
  ggsave(paste0("ranking_ranking_eval_IIA/step4_order_",i,".pdf"),
         width = 10, height = 8)
}

#----
# Eval step 5
#----
compa_step5_bef = step4[,c("cat_score","candidate","score4")] %>%
  rename(score=score4)
compa_step5_aft = step5[,c("candidate","score5")] %>%
  rename(score=score5)
compa_step5_aft$cat_score = "Aggregated"

compa_step5 <-  bind_rows(compa_step5_bef,compa_step5_aft)
compa_step5$cattype <- sapply(compa_step5$cat_score, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step5$cattype <- factor(compa_step5$cattype, levels=c("Single","Aggregated"))

ggplot(compa_step5, aes(x=cattype, y=score, color=cat_score)) +
          geom_point() +
          facet_wrap(~candidate) +
          theme_modern()
ggsave("ranking_ranking_eval_IIA/step5.pdf",
         width = 10, height = 8)

compa_step5_bef_order = compa_step5_bef %>%
  group_by(cat_score) %>%
  mutate(rank=order(score)) %>%
  select(cat_score,candidate,rank)
compa_step5_aft_order = compa_step5_aft %>%
  mutate(rank=order(score)) %>%
  select(candidate,rank)
compa_step5_aft_order$cat_score = "Aggregated"

compa_step5_order <-  bind_rows(compa_step5_bef_order,compa_step5_aft_order)
compa_step5_order$cattype <- sapply(compa_step5_order$cat_score, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step5_order$cattype <- factor(compa_step5_order$cattype, levels=c("Single","Aggregated"))

order_change = sapply(unique(compa_step5_order$candidate), function(y) {
    df = compa_step5_order[compa_step5_order$candidate==y,]
    df_bef = df$rank[df$cattype=="Single"]
    df_aft = df$rank[df$cattype=="Aggregated"]
    df_aft <= max(df_bef) & df_aft >= min(df_bef)
})
mean(order_change)
compa_step5_order$OK = NA
for (y in unique(compa_step5_order$candidate)) {
  compa_step5_order$OK[compa_step5_order$candidate==y] = order_change[y]
}
ggplot(compa_step5_order, aes(x=cattype, y=rank, color=cat_score)) +
          geom_point() +
          facet_wrap(~candidate) +
          title(i) +
          ggrepel::geom_label_repel(data=compa_step5_order[compa_step5_order$cattype=="Aggregated",], aes(label=OK)) +
          theme_modern()
ggsave("ranking_ranking_eval_IIA/step5_order.pdf",
         width = 10, height = 8)

#----
# Eval step 4b
#----
compa_step4b_bef = df_median_centerscale_transfo[,c("name_score","cat_score","candidate","dataset","trendval")] %>%
  rename(score=trendval)
compa_step4b_aft = step4b[,c("cat_score","candidate","dataset","score4b")] %>%
  rename(score=score4b)
compa_step4b_aft$name_score = "Aggregated"

compa_step4b <-  bind_rows(compa_step4b_bef,compa_step4b_aft)
compa_step4b$scoretype <- sapply(compa_step4b$name_score, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step4b$scoretype <- factor(compa_step4b$scoretype, levels=c("Single","Aggregated"))

for (i in unique(compa_step4b$cat_score)) {
  print(ggplot(compa_step4b[compa_step4b$cat_score==i,], aes(x=scoretype, y=score, color=name_score)) +
          geom_point(aes(shape=dataset)) +
          facet_wrap(~candidate) +
          theme_modern())
  ggsave(paste0("ranking_ranking_eval_IIA/step4b_",i,".pdf"),
         width = 10, height = 8)
}

compa_step4b_bef_order = compa_step4b_bef %>%
  group_by(name_score,cat_score,dataset) %>%
  mutate(rank=order(score)) %>%
  select(name_score,cat_score,candidate,dataset,rank)
compa_step4b_aft_order = compa_step4b_aft %>%
  group_by(cat_score,dataset) %>%
  mutate(rank=order(score)) %>%
  select(cat_score,candidate,dataset,rank)
compa_step4b_aft_order$name_score = "Aggregated"

compa_step4b_order <-  bind_rows(compa_step4b_bef_order,compa_step4b_aft_order)
compa_step4b_order$scoretype <- sapply(compa_step4b_order$name_score, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step4b_order$scoretype <- factor(compa_step4b_order$scoretype, levels=c("Single","Aggregated"))

order_change = lapply(unique(compa_step4b_order$cat_score), function(x)
  sapply(unique(compa_step4b_order$dataset), function(z)
    sapply(unique(compa_step4b_order$candidate), function(y) {
      df = compa_step4b_order[compa_step4b_order$cat_score==x & compa_step4b_order$dataset==z & compa_step4b_order$candidate==y,]
      df_bef = df$rank[df$scoretype=="Single"]
      df_aft = df$rank[df$scoretype=="Aggregated"]
      df_aft <= max(df_bef) & df_aft >= min(df_bef)
      })))
names(order_change) <- unique(compa_step4b_order$cat_score)
mean(sapply(order_change,mean))
compa_step4b_order$OK = NA
for (x in unique(compa_step4b_order$cat_score)) {
  for (y in unique(compa_step4b_order$candidate)) {
    for (z in unique(compa_step4b_order$dataset)) {
      compa_step4b_order$OK[compa_step4b_order$dataset==z & compa_step4b_order$candidate==y & compa_step4b_order$cat_score==x] = order_change[[x]][y,z]
    }
  }
}
for (i in unique(compa_step4b_order$cat_score)) {
  print(ggplot(compa_step4b_order[compa_step4b_order$cat_score==i,], aes(x=scoretype, y=rank, color=name_score)) +
          geom_point(aes(shape=dataset)) +
          facet_wrap(~candidate) +
          title(i) +
          ggrepel::geom_label_repel(data=compa_step4b_order[compa_step4b_order$cat_score==i & compa_step4b_order$scoretype=="Aggregated",], aes(label=OK)) +
          theme_modern())
  ggsave(paste0("ranking_ranking_eval_IIA/step4b_order_",i,".pdf"),
         width = 10, height = 8)
}

#----
# Eval step 5b
#----
compa_step5b_bef = step4b[,c("dataset","cat_score","candidate","score4b")] %>%
  rename(score=score4b)
compa_step5b_aft = step5b[,c("dataset","candidate","score5b")] %>%
  rename(score=score5b)
compa_step5b_aft$cat_score = "Aggregated"

compa_step5b <-  bind_rows(compa_step5b_bef,compa_step5b_aft)
compa_step5b$cattype <- sapply(compa_step5b$cat_score, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step5b$cattype <- factor(compa_step5b$cattype, levels=c("Single","Aggregated"))

for (i in unique(compa_step5b$dataset)) {
  print(ggplot(compa_step5b[compa_step5b$dataset==i,], aes(x=cattype, y=score, color=cat_score)) +
          geom_point() +
          facet_wrap(~candidate) +
          theme_modern())
  ggsave(paste0("ranking_ranking_eval_IIA/step5b_",i,".pdf"),
         width = 10, height = 8)
}

compa_step5b_bef_order = compa_step5b_bef %>%
  group_by(dataset,cat_score) %>%
  mutate(rank=order(score)) %>%
  select(dataset,cat_score,candidate,rank)
compa_step5b_aft_order = compa_step5b_aft %>%
  group_by(dataset) %>%
  mutate(rank=order(score)) %>%
  select(dataset,candidate,rank)
compa_step5b_aft_order$cat_score = "Aggregated"

compa_step5b_order <-  bind_rows(compa_step5b_bef_order,compa_step5b_aft_order)
compa_step5b_order$cattype <- sapply(compa_step5b_order$cat_score, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step5b_order$cattype <- factor(compa_step5b_order$cattype, levels=c("Single","Aggregated"))

order_change = sapply(unique(compa_step5b_order$dataset), function(x)
  sapply(unique(compa_step5b_order$candidate), function(y) {
    df = compa_step5b_order[compa_step5b_order$dataset==x & compa_step5b_order$candidate==y,]
    df_bef = df$rank[df$cattype=="Single"]
    df_aft = df$rank[df$cattype=="Aggregated"]
    df_aft <= max(df_bef) & df_aft >= min(df_bef)
  }))
mean(order_change)
compa_step5b_order$OK = NA
for (x in unique(compa_step5b_order$dataset)) {
  for (y in unique(compa_step5b_order$candidate)) {
    compa_step5b_order$OK[compa_step5b_order$candidate==y & compa_step5b_order$dataset==x] = order_change[y,x]
  }
}
for (i in unique(compa_step5b_order$dataset)) {
  print(ggplot(compa_step5b_order[compa_step5b_order$dataset==i,], aes(x=cattype, y=rank, color=cat_score)) +
          geom_point() +
          facet_wrap(~candidate) +
          title(i) +
          ggrepel::geom_label_repel(data=compa_step5b_order[compa_step5b_order$cat_score==i & compa_step5b_order$cattype=="Aggregated",], aes(label=OK)) +
          theme_modern())
  ggsave(paste0("ranking_ranking_eval_IIA/step5b_order_",i,".pdf"),
         width = 10, height = 8)
}

#----
# Eval step 3b
#----
compa_step3b_bef = step5b[,c("dataset","candidate","score5b")] %>%
  rename(score=score5b)
compa_step3b_aft = step3b[,c("candidate","score3b")] %>%
  rename(score=score3b)
compa_step3b_aft$dataset = "Aggregated"

compa_step3b <-  bind_rows(compa_step3b_bef,compa_step3b_aft)
compa_step3b$datatype <- sapply(compa_step3b$dataset, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step3b$datatype <- factor(compa_step3b$datatype, levels=c("Single","Aggregated"))

ggplot(compa_step3b, aes(x=datatype, y=score, color=dataset)) +
  geom_point() +
  facet_wrap(~candidate) +
  theme_modern()
ggsave("ranking_ranking_eval_IIA/step3b.pdf",
       width = 10, height = 8)

compa_step3b_bef_order = compa_step3b_bef %>%
  group_by(dataset) %>%
  mutate(rank=order(score)) %>%
  select(dataset,candidate,rank)
compa_step3b_aft_order = compa_step3b_aft %>%
  mutate(rank=order(score)) %>%
  select(candidate,rank)
compa_step3b_aft_order$dataset = "Aggregated"

compa_step3b_order <-  bind_rows(compa_step3b_bef_order,compa_step3b_aft_order)
compa_step3b_order$datatype <- sapply(compa_step3b_order$dataset, function(x)
  ifelse(x=="Aggregated","Aggregated","Single"))
compa_step3b_order$datatype <- factor(compa_step3b_order$datatype, levels=c("Single","Aggregated"))

order_change = sapply(unique(compa_step3b_order$candidate), function(y) {
  df = compa_step3b_order[compa_step3b_order$candidate==y,]
  df_bef = df$rank[df$datatype=="Single"]
  df_aft = df$rank[df$datatype=="Aggregated"]
  df_aft <= max(df_bef) & df_aft >= min(df_bef)
})
mean(order_change)
compa_step3b_order$OK = NA
for (y in unique(compa_step3b_order$candidate)) {
  compa_step3b_order$OK[compa_step3b_order$candidate==y] = order_change[y]
}
ggplot(compa_step3b_order, aes(x=datatype, y=rank, color=dataset)) +
  geom_point() +
  facet_wrap(~candidate) +
  ggrepel::geom_label_repel(data=compa_step3b_order[compa_step3b_order$datatype=="Aggregated",], aes(label=OK)) +
  theme_modern()
ggsave("ranking_ranking_eval_IIA/step3b_order.pdf",
       width = 10, height = 8)

ggplot(compa_step3b, aes(x=dataset, y=score, color=candidate)) +
  geom_point() +
  geom_line(aes(group=candidate)) +
  ggrepel::geom_label_repel(data=compa_step3b[compa_step3b$dataset=="Aggregated",], aes(label=candidate)) +
  theme_modern()
ggsave("ranking_ranking_eval_IIA/step3b_order2.pdf",
       width = 10, height = 8)

