
# set up the environment

import pickle
from matplotlib_venn import venn2
from onto_vae.ontobj import * 

project_dir = os.getcwd()

# import results of paired Wilcoxon testing 
res = pd.read_csv(project_dir + '/03_analysis/04_IFN_response/results/Genes_upregulating_IFN_activity_Wilcoxon_results.csv')

# import CD4T stimulated genes
diff_genes = pd.read_csv(project_dir + '/03_analysis/04_IFN_response/results/CD4T_IFN-ß_stim_up.txt', header=None)
diff_genes = diff_genes.iloc[:,0].tolist()



# make Fig4e (Hockeystick plot)

res_plot = res.iloc[0:717,:]   # 717 is the rank from GSEA
res_plot['log10p'] = -np.log10(res_plot['w_pval'])
res_plot['loglog10p'] = np.log2(res_plot['log10p'])
res_plot = res_plot.sort_values('loglog10p', ascending=False)
res_plot['rank'] = np.arange(1,res_plot.shape[0]+1)
res_plot['diff_rank'] = [diff_genes.index(g) + 1 if g in diff_genes else 'none' for g in res_plot.gene.tolist()]

# create color map for differential rank
color_map = dict(zip(list(np.arange(1,len(diff_genes))[::-1]), sns.color_palette("Reds", n_colors=len(diff_genes))))
color_map['none'] = 'black'

res_plot['grad_color'] = res_plot['diff_rank'].map(color_map)


g = sns.lineplot(data=res_plot,
                 x=res_plot['rank'],
                 y=res_plot['loglog10p'],
                 color='black')
g.set(ylabel='log2(-log10(p-value))')


g = sns.scatterplot(data=res_plot[res_plot.gene.isin(diff_genes)],
                    x='rank',
                    y='loglog10p',
                    hue='diff_rank',
                    palette=color_map)
plt.legend([],[], frameon=False)

for i, g in enumerate(res_plot.gene.to_list()):
    if g in diff_genes:
        plt.text(i + 1, 
                res_plot['loglog10p'].iloc[i], 
                g, 
                color=res_plot['grad_color'].iloc[i]) 

plt.savefig(project_dir + '/03_analysis/04_IFN_response/figures/06-01_ranking_influence_of_in_silico_stimulated_genes_on_IFN_signalling_leading_edge_Fig3c.pdf')
plt.clf()


# import ontobj 
with open(project_dir + '/ontobj/GO_symbol_ontobj.pickle', 'rb') as f:
    go = pickle.load(f) 

# get descendant genes dict
desc_genes = go.desc_genes['1000_30']

ifn_genes = desc_genes['GO:0060337'] # ID for type I interferon signaling pathway


# make Venn diagram

venn2([set(diff_genes), set(res_plot.gene.tolist())], ('diff genes', 'in silico genes'))
plt.savefig(project_dir + '/03_analysis/04_IFN_response/figures/06-02_overlap_diff_genes_in_silico_stimulated_genes_leading_edge.pdf')
plt.clf()