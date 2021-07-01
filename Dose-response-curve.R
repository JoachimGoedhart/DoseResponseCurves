# Fitting a dose-response curve


library(tidyverse)
library(broom)
library(drc)
library(modelr)

#Get a list with all csv files
filelist = list.files(pattern="*.csv$")

#read all csv files and put in df_input_list
df_input_list <- lapply(filelist, read.csv)

#get the filenames, remove extension and use as "id"
names(df_input_list) <- gsub(filelist, pattern="\\..*", replacement="")

#Merge all the dataframes and use the filenames as id    
df <- bind_rows(df_input_list, .id = "id")

#Seperate replicate and date
df <- df %>% separate(Experiment, c('data', 'replicate'), sep=" ")

#Set zero values to a 'low' value, to enable a log10 scale on x-axis
df <- df %>%
  mutate(Concen = ifelse((Concen == 0 & id == 'His nac'),
                  yes = 0.019,
                  no = Concen),
         Concen = ifelse((Concen == 0 & id == 'S1P nac'),
                         yes = 4,
                         no = Concen),
         Concen = ifelse((Concen == 0 & id == 'UK nac'),
                         yes = 0.04,
                         no = Concen)
         ) %>% filter(Condition!='ymptx')


df_summary <- df %>% group_by(Condition, Concen, id, replicate,) %>% summarise(mean_ERK=mean(ERKnac), mean_Akt=mean(Aktnac))


# https://stackoverflow.com/questions/36780357/plotting-dose-response-curves-with-ggplot2-and-drc/37043751

ggplot(data = df, aes(x = Concen, y = ERKnac)) + 
  geom_jitter(alpha=0.1, width=0.1,size=0.4) + 
  geom_point(data=df_summary, aes(x = Concen, y = mean_ERK, fill=replicate), size=8, shape=21, color='black', alpha=0.5) +
  geom_smooth(data=df_summary, aes(x = Concen, y = mean_ERK),color='black',method = drm, method.args = list(fct = L.4()), se = FALSE) +
  scale_x_log10()+ylim(-1,12)+
  facet_grid(Condition~id, scales="free_x") +theme_light()


##### Plots the Akt data if uncommented
# ggplot(data = df, aes(x = Concen, y = Aktnac)) + 
#   geom_jitter(alpha=0.1) + 
#   geom_point(data=df_summary, aes(x = Concen, y = mean_Akt, fill=replicate), size=8, shape=21, color='black', alpha=0.5) +
#   geom_smooth(data=df_summary, aes(x = Concen, y = mean_Akt),color='black',method = drm, method.args = list(fct = L.4()), se = FALSE) +
#   scale_x_log10()+ylim(-1,4)+
#   facet_grid(Condition~id, scales="free_x") +theme_light()


#### Scatter plot Erk vs. Akt
# ggplot(data = df, aes(x = Aktnac, y = ERKnac))+ 
#   # geom_point(alpha=0.1)+
#   geom_hex() + scale_fill_viridis_c() +
#   # geom_bin2d() +
#   facet_grid(id~Condition)+xlim(-2,5)+ylim(-2,12)

# https://gist.github.com/angelovangel/c69d78a27c4360df3057b40f4a705837
# define drm function to use with map
drm.func <- function(x) {
  drm(mean_ERK ~ Concen, 
      fct = LL.4(names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), 
      data = x)
}

drm.ci <- function(x) {confint(drm.func(x))}

coefs.fun <- function(x) {coef(x) %>% tidy}

coefs <- df_summary %>% group_by(id, Condition) %>% nest() %>%
  mutate(drmod = map(data, drm.func), co = map(drmod, coefs.fun), ci = map(data, drm.ci)
         )


# summary of coefficients
coefs %>% unnest(co) %>% spread(names,x)

# 95%CI
coefs$ci


