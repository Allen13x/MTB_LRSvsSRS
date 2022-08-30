library(rstatix)
library(tidyverse)


# Import quast data -------------------------------------------------------



read_delim('Assembly/quast/report.tsv',delim='\t')->quast

read_delim('Assembly/quast/genome_data.tsv',delim=' ')[-1,]->quast_gen

quast %>% 
  gather(Sample,val,-Assembly) %>% 
  separate(Sample,c('type','sample'),sep='_') %>% 
  mutate(type=factor(type,levels=c('SRS','HybA','Hybrid','LRS'))
  ) %>% 
  droplevels() %>% 
  group_by(type) %>% 
  filter(Assembly%in%c('Genome fraction (%)',
                       'Largest alignment',
                       'NG50','NA50',
                       '# misassembled contigs',
                       '# contigs')) %>% 
  mutate(val=as.numeric(val)) %>% ungroup()->quast_metric


quast_gen %>% 
  separate(assembly,c('type','sample'),sep='_') %>% 
  mutate(genes=100*as.numeric(genes)/3954) %>% 
  rename('partial'='partial...6') %>% 
  select(type,sample,genes,gaps,partial) %>% 
  gather(Assembly,val,-type,-sample) %>% 
  mutate(val=as.numeric(val))->quest_gmetric



# Statistics --------------------------------------------------------------

quast_metric %>% 
  rbind(quest_gmetric) %>% 
  group_by(Assembly) %>% 
  wilcox_test(val~type) %>% 
  add_xy_position(scales='free_y',step.increase = 0.2) %>% 
  droplevels() %>% ungroup()->quast_stat



# Visualization -----------------------------------------------------------



quast_lab<-c(
  'partial' = 'N. Partial genes',
  'gaps' = 'N. Gaps',
  'genes'='Genes Fraction (%)',
  'Genome fraction (%)' = 'Genome Fraction (%)',
  'Largest alignment'='Largest Alignement',
  'NG50' = 'NG50','NA50'='NA50',
  '# misassembled contigs' = 'N. Misassembled Contigs',
  '# contigs' = 'N. Contigs'
)


quast_metric %>% 
  rbind(quest_gmetric) %>% 
  ggplot(aes(x=type,y=val,fill=type))+
  geom_boxplot(outlier.shape=NA,width=0.5)+
  geom_jitter(shape=21,fill='white',height=0,width=0.3)+
  theme_classic()+
  scale_fill_manual(values=c('#EC9B28',
                             'olivedrab3',
                             '#8A0845','#0083AC'))+
  facet_wrap(~Assembly,scale='free_y',labeller = labeller(Assembly = quast_lab))+
  ylab('')+xlab('')+
  theme(axis.text.x = element_text(angle=45,vjust=.5),
        text=element_text(size=12))+
  stat_pvalue_manual(quast_stat %>% droplevels(),label='p.adj.signif',hide.ns=T)

dir.create('FIGS')

ggsave('FIGS/assembly_stat.png',dpi=300)

