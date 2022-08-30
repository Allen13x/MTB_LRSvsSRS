library(rstatix)
library(tidyverse)


# Coverage 1k ----------------------------------------------------------------

sample_list<- read_delim('sample_list',delim='\t')
cov=data.frame()
for (i in sample_list){
  for (j in c('SRS','LRS','Hybrid')){
    cov=rbind(cov, read_delim(paste('Coverage/',j,'_',i,'_1k.thresholds.bed.gz',sep=''),delim='\t')%>%
                mutate(samp=i,r=j))
    
  }
}


cov%>%mutate(l=end-start)%>%
  separate(region,c('seq','n'),sep='_') %>% 
  mutate(n=str_pad(n,3,'left',pad = 0)) %>% 
  unite(seq:n,col = 'region',sep='_') %>% 
  mutate_at(c("1X","2X","4X","6X","8X","12X","16X","20X","40X"),~.x/l*100)%>%
  group_by(r,samp)%>%summarise_at(c("1X","2X","4X","6X","8X","12X","16X","20X","40X"),mean)%>%
  select(Tech=r,samp,coverage8=`8X`) %>% 
  ungroup() ->cov_8
  anova_test(coverage8~Tech) %>% pull(p)

cov_8 %>% 
  t_test(coverage8~Tech) %>% adjust_pvalue() %>% 
  add_significance(p.col='p.adj') %>% 
  add_y_position(step.increase =0.5, scales = 'free_y') %>% 
  add_x_position()->bar_pval



cov_8 %>% 
  mutate(Tech=factor(Tech,levels=c('SRS','Hybrid','LRS'))) %>% 
  ggplot(aes(x=Tech,y=coverage8)) +
  geom_errorbar(stat='summary',aes(col=Tech),width=0.1,lwd=1)+
  geom_bar(aes(fill=Tech,col=Tech),stat = "summary",fun='mean',width=0.5,lwd=1.3)+
  
  scale_y_continuous(limits=c(95,102),oob = scales::rescale_none)+
  theme_classic()+
  stat_pvalue_manual(bar_pval %>% droplevels(),hide.ns=T,tip.length = 0.0,step.increase = 0.01)+
  scale_color_manual(values=c('#EC9B28',
                              '#8A0845','#0083AC'),
                     label=c('SRS',
                             'Hybrid','LRS'),name='')+
  scale_fill_manual(values=alpha(c('#EC9B28',
                                   '#8A0845','#0083AC'),0.6),
                    label=c('SRS',
                            'Hybrid','LRS'),name='')+
  ylab('Breadth coverage (%)')+xlab('')+
  theme(axis.text.x = element_text(angle=45,vjust=.5),
        text=element_text(size=12))


dir.create('FIGS')


ggsave('FIGS/Coverage_barplot.png',dpi=300)
