# Analysis for conflict resolution project


# Load libraries ----------------------------------------------------------
library(leadershipdata)
library(tidyverse)
library(ggmosaic)
library(viridis)
library(glmnet)
library(lme4)
library(emmeans)
library(car)
library(hagenutils)
library(rstan)
library(rstanarm)
library(modelbased)
library(ggeffects)
library(stargazer)


# Create functions --------------------------------------------------------

merge_dfs <- function(
  leader_text2, 
  all_ids, 
  leader_cult,
  documents,
  threshold = 4,
  vars = c(
    "c_culture_code",
    "region", 
    "subsistence", 
    "c_cultural_complexity", 
    "settlement_fixity", 
    "pop_density", 
    "com_size",
    "ehraf_pages",
    "number_leader_records"
  )
){
  leader_text2 %>%
    select_if(
      ~ is.character(.) | (is.numeric(.) && sum(., na.rm = T) >= threshold) # remove low evidence
    ) %>% 
    left_join(all_ids, by = 'cs_textrec_ID') %>%
    left_join(
      leader_cult[vars], by = c('d_culture' = 'c_culture_code')
    ) %>% 
    left_join(
      documents[c('d_ID', 'd_publication_date', 'female_coauthor')], # One doc missing pub date
      by = c('doc_ID' = 'd_ID')
    ) %>% 
    rename(
      pub_date = d_publication_date
    ) %>% 
    mutate(
      pagesZ = (ehraf_pages - mean(ehraf_pages, na.rm=T))/sd(ehraf_pages, na.rm=T),
      pagesSD = ehraf_pages/sd(ehraf_pages, na.rm=T),
      pub_date = as.numeric(pub_date),
      pub_dateZ = (pub_date - mean(pub_date, na.rm=T))/sd(pub_date, na.rm = T)
    )
}

model_words <- function(d, dtm, var, lam = 'lambda.min', exponentiate = T, title=''){
  
  df <- left_join(d[c('cs_textrec_ID', var)], dtm, by = 'cs_textrec_ID')
  y <- df[[2]]
  x <- as.matrix(df[-c(1:2)])
  
  m_cv <- cv.glmnet(x, y, family = 'binomial', alpha = 1, standardize = F)
  plot(m_cv)
  
  if (lam == 'mid'){
    lmda <- m_cv$lambda.min + (m_cv$lambda.1se - m_cv$lambda.min)/2
  } else if (lam == 'min'){
    lmda <- m_cv$lambda.min
  } else if(lam == '1se'){
    lmda <- m_cv$lambda.1se
  } else {
    lmda = lam
  }
  
  c.min <- coef(m_cv, s = lmda)
  
  coefs <- c()
  for (i in 1:length(c.min)){
    if (c.min[i] != 0){
      coefs <- c(coefs, c.min[i])
      names(coefs)[length(coefs)] <- rownames(c.min)[i]
    }
  }
  
  if (exponentiate){
    coefs <- exp(coefs)
    xintrcpt <- 1
  } else {
    xintrcpt <- 0
  }
  
  coefs <- sort(coefs[-1]) # delete intercept
  
  df <-
    tibble(
      Word = factor(names(coefs), levels = names(coefs)),
      Coefficient = coefs,
      Sign = ifelse(coefs > xintrcpt, 'Increase', 'Decrease')
    )
  
  plot <- 
    ggplot(df, aes(Coefficient, Word, colour = Sign, shape=Sign)) + 
    geom_point(size=3) + 
    geom_vline(xintercept = xintrcpt, linetype = 'dotted') +
    hagenutils::scale_color_binary() +
    guides(colour=F, shape=F) +
    labs(title = title, x = '', y = '') +
    theme_minimal(15)
  
  return(plot)
}

# Elasticnet models of one dimension by other dimensions

elastic_dimensions <- function(d, outcomevar, predictorvars, alpha = 1, lambda = 'lambda.min', threshold = 0){
  predvars <- variable_names(d, predictorvars)
  predvars <- predvars[predvars != outcomevar]
  
  y <- d[[outcomevar]]
  x <- d[predvars]
  x <- x[colSums(x)>threshold]
  x <- as.matrix(x)
  
  m <- glmnet::cv.glmnet(x, y, family = 'binomial', alpha = alpha)
  plot(m)
  coefs <- coef(m, s = m[[lambda]])[-1,1] # delete intercept
  names(coefs) <- var_names[names(coefs)] # var_names from leadershipdata
  ggdotchart(exp(coefs[coefs != 0]), threshold = 1) +
    geom_vline(xintercept = 1, linetype = 'dotted') +
    guides(colour=F, shape=F) +
    scale_x_log10()
}

# variable names encode their type
variable_names <- function(df, type){
  dfnames <- names(df)
  dfnames <- dfnames[!str_detect(dfnames, "functions_Context")]
  
  # build regex
  pattern <- str_c(type, collapse = '|')
  
  thevars <- dfnames[str_detect(dfnames, pattern)]
  names(thevars) <- var_names[thevars]
  return(thevars)
}

# logistic model EMM plot
all_emms <- function(m, specs, upperlimit, title){
  
  theplots = list()
  for (spec in specs){
    p <- 
      hagenutils::ggemmeans(emmeans(m, spec, type = 'response')) +
      scale_x_continuous(limits = c(0, upperlimit)) +
      labs(title = '', x = '', y = '')
    theplots <- c(theplots, list(p))
  }
  return(theplots)
}

# Fucntion to create BC ratio
benefit_cost_ratio <- function(d){
  m <-
    glmer(
      Evidence ~
        Benefit_cost * Variable +
        (1|d_culture/author_ID),
      family = binomial,
      data = d,
      nAGQ = 0
    )
  
  em_cb <- summary(emmeans(m, pairwise ~ Benefit_cost, type = 'response'))
  em_cb_var <- confint(emmeans(m, pairwise ~ Benefit_cost | Variable, type = 'response'))
  
  list(
    benefit_cost_OR = em_cb$contrasts$odds.ratio,
    benefit_cost_var_OR = em_cb_var$contrasts
  )
}



# Make new varnames -------------------------------------------------------

var_names2 <- c(
  "functions_BestowMate" = "Bestow mates",
  "functions_PoliticalAppointments" = "Political appointments",
  "functions_ConstructionInfrastructure"       = "Construction/infastructure",
  "functions_ControlEconomics"                 = "Control economics",
  "functions_CouncilMember"                    = "Council member",
  "functions_GroupDetermination"              = "Group determination/cohesiveness",
  "functions_Hospitality"                       = "Hospitality",
  "functions_MilitaryCommand"                  = "Military command",
  "functions_NewSettlement"                    =  "Movement/migration",
  "functions_ProsocialInvestment"              = "Prosocial investment",
  "functions_ProvideCounsel"                   = "Provide counsel/direction",
  "functions_Punishment"                        = "Punishment",
  "functions_ServeLeader"                      = "Serve a leader",
  "functions_StrategicPlanning" = "Strategic planning",
  "functions_OrganizeCooperation"  = "Organize cooperation",
  "functions_ResolveConflcit"      = "Resolve conflict",
  "functions_ControlCalendar"     = "Control calendar",
  "functions_ControlImmigration"  = "Control immigration",
  "functions_DistributeResources" = "Distribute resources",
  "functions_GroupRepresentative" = "Group representative",
  "functions_Medicinal"           = "Medicinal functions",
  "functions_MoralAuthority"     = "Moral authority",
  "functions_Policymaking"        =  "Policy making",
  "functions_Protection"          = "Protection",
  "functions_ProvideSubsistence" = "Provide subsistence",
  "functions_Ritual"              = "Ritual functions",
  "functions_SocialFunctions" = "Misc. social functions",
  "qualities_ArtisticPerformance"     = "Artistic performance",
  "qualities_Generous"                 = "Generosity",
  "qualities_Age"                      = "Age",
  "qualities_Attractive"              = "Attractive",
  "qualities_CoerciveAuthority"      = "Coercive authority",
  "qualities_CulturallyProgressive"  = "Culturally progressive",
  "qualities_FavorablePersonality"   = "Favorable personality",
  "qualities_Honest"                  = "Honesty",
  "qualities_IngroupMember"          = "Ingroup member",
  "qualities_Killer"                  = "Killer",
  "qualities_ManyChildren"           = "Many children",
  "qualities_PhysicallyStrong"       = "Physically formidable",
  "qualities_Prosocial"               = "Prosocial",
  "qualities_StrategicPlanner"       = "Strategic planner",
  "qualities_DrugConsumption"     = "Drug consumption",
  "qualities_HighStatus"            = "High status",
  "qualities_Aggressive"             = "Aggressiveness",
  "qualities_Bravery"               = "Bravery",
  "qualities_Confident"             = "Confidence",
  "qualities_Decisive"              = "Decisiveness/decision-making",
  "qualities_Feared"                = "Feared",
  "qualities_Humble"                = "Humility",
  "qualities_Innovative"            = "Innovative",
  "qualities_KnowlageableIntellect" = "Knowledgeable/intelligent",
  "qualities_OratorySkill"         = "Oratory skill",
  "qualities_Polygynous"            = "Polygynous",
  "qualities_SocialContacts"       = "Social contacts",
  "qualities_Supernatural"    = "Supernatural",
  "qualities_ExpAccomplished"       = "Experienced/accomplished",
  "qualities_Wealthy"                = "Wealthy",
  "qualities_Ambition"               = "Ambitious",
  "qualities_Charisma"               = "Charisma",
  "qualities_CulturallyConservative" = "Culturally conservative",
  "qualities_Fairness"               = "Fairness",
  "qualities_HighQualitySpouse"    = "High-quality spouse",
  "qualities_Industriousness"        = "Industriousness",
  "qualities_InterpersonalSkills"   = "Interpersonal skills",
  "qualities_Loyalty"                = "Loyalty",
  "qualities_PhysicalHealth"        = "Physical health",
  "qualities_ProperBehavior"        = "Proper behavior",
  "qualities_StrategicNepotism"     = "Strategic nepotism",
  "qualities_Xenophobic"   = "Xenophobia",
  "qualities_AntiHonest" = "Dishonest",
  "qualities_AntiFairness" = "Unfair",
  "qualities_AntiDrugConsumption" = "No drug consumption",
  "qualities_AntiCoerciveAuthority" = "No coercive authority",
  leader.benefits_Fitness = 'Leader benefit: Inclusive fitness',
  leader.benefits_Mating = 'Leader benefit: Mating',
  leader.benefits_Other = 'Leader benefit: Misc. non-material',
  leader.benefits_RiskHarmConflict = 'Leader benefit: Reduced risk of harm',
  leader.benefits_ResourceFood = 'Leader benefit: Food',
  leader.benefits_ResourceOther = 'Leader benefit: Material resources',
  leader.benefits_SocialServices = 'Leader benefit: Social services',
  leader.benefits_SocialStatusReputation = 'Leader benefit: Increased social status',
  leader.benefits_Territory = 'Leader benefit: Territory',
  follower.benefits_Fitness = 'Follower benefit: Inclusive fitness',
  follower.benefits_Mating = 'Follower benefit: Mating',
  follower.benefits_Other = 'Follower benefit: Misc. non-material',
  follower.benefits_RiskHarmConflict = 'Follower benefit: Reduced risk of harm',
  follower.benefits_ResourceFood = 'Follower benefit: Food',
  follower.benefits_ResourceOther = 'Follower benefit: Material',
  follower.benefits_SocialServices = 'Follower benefit: Social services',
  follower.benefits_SocialStatusReputation = 'Follower benefit: Increased social status',
  follower.benefits_Territory = 'Follower benefit: Territory',
  
  leader.costs_Fitness = 'Leader cost: Inclusive fitness',
  leader.costs_RiskHarmConflict = 'Leader cost: Increased risk of harm',
  leader.costs_Other = 'Leader cost: Misc. non-material',
  leader.costs_ResourceFood = 'Leader cost: Food',
  leader.costs_ResourceOther = 'Leader cost: Loss of material resources',
  leader.costs_SocialStatusReputation = 'Leader cost: Reduced social status',
  leader.costs_Territory = 'Leader cost: Loss of territory',
  leader.costs_Mating = 'Leader cost: Mating cost',
  leader.costs_SocialServices = 'Leader cost: Loss of social services',
  follower.costs_Fitness = 'Follower cost: Inclusive fitness',
  follower.costs_RiskHarmConflict = 'Follower cost: Increased risk of harm',
  follower.costs_Mating = 'Follower cost: Mating',
  follower.costs_Other = 'Follower cost: Misc. non-material',
  follower.costs_ResourceFood = 'Follower cost: Food',
  follower.costs_ResourceOther = 'Follower cost: Loss of material resources',
  follower.costs_SocialServices = 'Follower cost: Loss of social services',
  follower.costs_SocialStatusReputation = 'Follower cost: Reduced social status',
  follower.costs_Territory = 'Follower cost: Loss of territory'
)

# Prep data ---------------------------------------------------------------

# df that links the ids of texts, docs, and cultures -----------------------
text_doc_auth_cult_ID <- function(df_text, df_doc){
  df_text %>% 
    dplyr::select(cs_textrec_ID, doc_ID, author_ID) %>% 
    left_join(df_doc[c('d_ID', 'd_culture')], by = c("doc_ID" = "d_ID"))
}

all_ids = text_doc_auth_cult_ID(leader_text_original, documents)
all_data = merge_dfs(leader_text2, all_ids, leader_cult, documents, threshold = 1)

all_data_conflict = all_data[all_data$functions_ResolveConflcit==1,]

# Moasic plot of Conflict resolution context by group context  -----------------------------------------------------------

df_groups_conf <- 
  all_data %>% 
  dplyr::select(
    group.structure2,
    demo_sex,
    subsistence,
    functions_ResolveConflcit,
    functions_Context
  ) %>% 
  dplyr::mutate(
    sex = factor(demo_sex, levels = c('male', 'female')),
    group = factor(
      group.structure2,
      levels = c(
        'residential subgroup',
        'kin group',
        'economic group',
        'religious group',
        'military group',
        'political group (community)',
        'political group (supracommunity)'
      ),
      plyr::revalue(c(
        'residential subgroup'='residential\nsubgroup',
        'kin group'='kin\ngroup',
        'economic group'='economic\ngroup',
        'religious group'='\n\nreligious\ngroup',
        'military group'='military\ngroup',
        'political group (community)'='political group\n(community)',
        'political group (supracommunity)'='political group\n(supracommunity)'
      ))
    ),
    subsistence = factor(
      subsistence,
      levels = c("hunter gatherers",
                 "pastoralists",
                 "mixed",
                 "horticulturalists",
                 "agriculturalists"
      )
    ),
    context = factor(
      functions_Context,
      levels = c(
        "within-group",
        "between-group",
        "both"
      )
    )
    
  )

  
## REMOVE 4 NA's of function contexts for now
df_groups_conf <- df_groups_conf[!is.na(df_groups_conf$context),]

## Remove all Unkown from this data frame

plot_sub_context <-
  ggplot(df_groups_conf) +
  geom_mosaic(aes(x = product(context, subsistence), fill = context)) +
  scale_fill_viridis(discrete = T) +
  labs(x="", y="", fill = "Context of conflict") +
  guides(fill = guide_legend(reverse = T)) +
  theme_minimal(20) 
plot_sub_context
ggsave("subsistence.pdf", plot_sub_context, width = 17)


plot_group_context <-
  ggplot(df_groups_conf) +
  geom_mosaic(aes(x = product(context, group), fill = context)) +
  scale_fill_viridis(discrete = T) +
  labs(x="", y="", fill = "Context of conflict") +
  guides(fill = guide_legend(reverse = T)) +
  theme_minimal(20)
plot_group_context
ggsave("context.pdf", plot_group_context, width = 15)



# # Heatmap/clust of conflict data ------------------------------------------------
# 
# all_data_conflict_heat <- all_data_conflict[,c(3:23,25:113,115,120)]
# 
# all_data_conflict_heat_numeric <- all_data_conflict_heat[,c(1:22,24:45,47:110)]
# 
# all_data_conflict_heat_numeric <- all_data_conflict_heat_numeric[!rowSums(all_data_conflict_heat_numeric)==0,]

#conflict_clust <- pvclust::pvclust(all_data_conflict_heat_numeric, method.hclust = "ward.D", method.dist = "euclidean", 
                                  # nboot = 1000, parallel = TRUE)
# plot(conflict_clust)

# Text analysis -----------------------------------------------------------

conflict_words_plot = model_words(all_data, leader_dtm, 'functions_ResolveConflcit', lam = "1se", title = '')
conflict_words_plot

# elasticnet --------------------------------------------------------------

plot_elastic_conflict = elastic_dimensions(all_data, 'functions_ResolveConflcit', c('functions', 'qualities'), alpha = 1, lambda = 'lambda.1se')
plot_elastic_conflict <- plot_elastic_conflict +
  scale_color_viridis_d()
plot_elastic_conflict

elastic_fit_conflict_cb = glmnet::cv.glmnet(as.matrix(all_data[c(3:20,50:67)]),
                                             all_data$functions_ResolveConflcit,
                                             family = 'binomial',
                                             alpha = 1)

coefs_cb <- coef(elastic_fit_conflict_cb, s = elastic_fit_conflict_cb[["lambda.min"]])[-1,1] # delete intercept

names(coefs_cb) <- var_names2[names(coefs_cb)] # var_names from leadershipdata

plot_elastic_conflict_cb <- ggdotchart(exp(exp(coefs_cb[coefs_cb != 0])), threshold = 1) +
  geom_vline(xintercept = 1, linetype = 'dotted') +
  guides(colour=F, shape=F) +
  scale_color_viridis_d() +
  scale_x_log10()
plot_elastic_conflict_cb

# Logistic model ----------------------------------------------------------

mm_conflict_group_int <-
  glmer(
    functions_ResolveConflcit ~
      (1|d_culture/author_ID),
    family = binomial,
    data = all_data,
    nAGQ = 0
  )

mm_conflict_group3 <-
  glmer(
    functions_ResolveConflcit ~
      group.structure2 +
      subsistence +
      region +
      (1|d_culture/author_ID),
    family = binomial,
    data = all_data,
    nAGQ = 0
  )


Anova(mm_conflict_group3)

car::Anova(mm_conflict_group3, type=c(2))

mm_conflict_group <-
  glmer(
    functions_ResolveConflcit ~
      group.structure2 +
      subsistence +
      (1|d_culture/author_ID),
    family = binomial,
    data = all_data,
    nAGQ = 0
  )

mm_conflict_group_1 <-
  glmer(
    functions_ResolveConflcit ~
      group.structure2 +
      (1|d_culture/author_ID),
    family = binomial,
    data = all_data,
    nAGQ = 0
  )

aic_INT <- AIC(mm_conflict_group_int)
aic_3pred <- AIC(mm_conflict_group3)
aic_FIN <- AIC(mm_conflict_group)
aic_ONE <- AIC(mm_conflict_group_1)



AIC(mm_conflict_group_int, mm_conflict_group3, mm_conflict_group, mm_conflict_group_1)

Anova(mm_conflict_group3)

theplots <- all_emms(mm_conflict_group3, c('region', 'subsistence', 'group.structure2'), 0.95, 'Predictors of conflict resolution')

modelplot_reg <- theplots[[1]]
modelplot_sub <- theplots[[2]]
modelplot_grp <- theplots[[3]]

## Adjust p-values
mm_conflict_group3_ADJP <- p.adjust(as.matrix(Anova(mm_conflict_group3)[3]), method = "BH")
mm_conflict_group_ADJP <- p.adjust(as.matrix(Anova(mm_conflict_group)[3]), method = "BH")

means<-emmeans(mm_conflict_group3, "group.structure2")
emmeans::contrast(means)
pairs(means)
pwpm(means)
plot(means, comparisons=TRUE)
contrasts_plot <- pwpp(means) +
  labs(y="", x = "\nTukey-adjusted P value") +
  scale_colour_viridis(discrete = T, direction = -1, option = "D") +
  theme_minimal(20)

# Logistic rstanarm model -------------------------------------------------

# m_conflict_group_stan <-
#   stan_glmer(
#     functions_ResolveConflcit ~
#       group.structure2 +
#       subsistence +
#       (1|d_culture/author_ID),
#     data = all_data,
#     family = binomial(link = "logit"),
#     prior = normal(0,1),
#     prior_intercept = normal(0,1),
#     
#     )
# 
# means <- estimate_means(m_conflict_group_stan)
# 
# marginal_eff <- ggeffect(m_conflict_group_stan)
# marginal_eff2 <- rstantools::posterior_predict(m_conflict_group_stan)
# 
# theplots_stan <- all_emms(m_conflict_group_stan, c('subsistence', 'group.structure2'), 0.8, 'Predictors of conflict resolution')
# 
# 
# pplot<-plot(m_conflict_group_stan, pars = c("group.structure2kin group",
#                                             "group.structure2military group",
#                                             "group.structure2political group (community)",
#                                             "group.structure2political group (supracommunity)",
#                                             "group.structure2religious group",
#                                             "group.structure2residential subgroup",
#                                             "group.economic group",
#                                             "subsistencehorticulturalists",
#                                             "subsistencehunter gatherers",
#                                             "subsistencemixed",
#                                             "subsistencepastoralists",
#                                             "subsistenceagriculturalits"
#                                             ),
#             "areas", prob = 0.90, prob_outer = .95)
# pplot+ geom_vline(xintercept = 0)


# Conflict evidence data --------------------------------------------------

conflict_data <- all_data[all_data$functions_ResolveConflcit == 1,] %>% 
  filter(functions_Context != "unkown") %>% 
  filter(functions_Context != "both")
 
conflict_data$functions_Context <- factor(conflict_data$functions_Context)

m_conflict_context <- glmer(
  functions_Context ~
    group.structure2 +
    subsistence +
    region +
    #demo_sex +
    (1|d_culture/author_ID),
  family = binomial,
  data = conflict_data,
  nAGQ = 0)


# Cost benefit ratio ------------------------------------------------------

df_costs_benefits =
  all_data_conflict %>% 
  select(d_culture, author_ID, contains('benefits'), contains('costs')) %>% 
  gather(key = Variable, value = Evidence, -d_culture, -author_ID) %>% 
  separate(Variable, into = c('Type', 'Variable'), sep = '\\.') %>% 
  separate(Variable, into = c('Benefit_cost', 'Variable'), sep = '_')
  #mutate(Variable = bc_dict[Variable]

# Leaders only
df_leader_costs_benefits =
  df_costs_benefits %>% 
  filter(Type == 'leader') %>% 
  select(-Type)

# Followers only
df_follower_costs_benefits =
  df_costs_benefits %>% 
  filter(Type == 'follower') %>% 
  select(-Type)

leader_bc_ratio = benefit_cost_ratio(df_leader_costs_benefits)
follower_bc_ratio = benefit_cost_ratio(df_follower_costs_benefits)

lvls <- leader_bc_ratio$benefit_cost_var_OR$Variable[order(leader_bc_ratio$benefit_cost_var_OR$odds.ratio)]

bc_OR <- 
  bind_rows(
    Leader = leader_bc_ratio$benefit_cost_var_OR, 
    Follower = follower_bc_ratio$benefit_cost_var_OR,
    .id = 'Type'
  ) %>% 
  mutate(
    Variable = factor(Variable, levels = lvls)
  )

plot_bc_OR <-
  ggplot(bc_OR, aes(odds.ratio, Variable, xmin = asymp.LCL, xmax = asymp.UCL, colour=Type)) + 
  geom_errorbarh(height = 0, lwd = 2.5, alpha = 0.7, position = position_dodge(0.7)) + 
  geom_point(position = position_dodge(0.7)) + 
  geom_vline(xintercept = 1, linetype = 'dotted') +
  hagenutils::scale_colour_binary() +
  scale_x_log10() +
  guides(colour = guide_legend(reverse=T)) +
  labs(title = 'Relative evidence of benefits vs. costs', x = '\nOdds ratio', y = '') +
  theme_minimal(15)
plot_bc_OR



# Playing with model sytax ------------------------------------------------

LMM <-
  lmer(
    functions_ResolveConflcit ~
      group.structure2 +
      subsistence +
      region +
      (1|d_culture/author_ID),
    data = all_data
  )

# Save Image --------------------------------------------------------------

save.image(file = "conflict.RData")
