import os
from os.path import join as opj
from multiprocessing import Pool
import pickle 
from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import scipy

def get_all_counts(remapped_diff_list, input_path):  # remapped_diff_list : [{input}/{patient}/{exp}/pipeline_output/remapped_diff_gene.pkl, ...]
    all_counts = {}

    strong_evidences_dict = {}
    retrocopy_clues_dict = {}
    alteration_clues_dict = {}
    homology_check_dict = {}

    for remapped_diff_path in remapped_diff_list:
        exp_path = os.path.dirname(os.path.dirname(remapped_diff_path))
        exp_name = os.path.basename(os.path.dirname(remapped_diff_path))

        print(f'starting {exp_name}...', flush=True)
        reads_stats = pickle.load(open(f'{exp_path}/pipeline_output/reads_stats.pkl', 'rb'))
        remapped_same_gene = pickle.load(open(f'{exp_path}/pipeline_output/remapped_same_gene.pkl', 'rb'))
        remapped_diff_gene = pickle.load(open(remapped_diff_path, 'rb'))
        all_counts[exp_name] = {
        'good_reads': reads_stats['good_reads'],
        'disc_reads': reads_stats['disc_reads'],
        'remapped_same_gene': len(remapped_same_gene),
        'remapped_diff_gene': len(remapped_diff_gene)}


        all_read_infos = {drp[0]['name']:drp for drp in remapped_diff_gene}
        read_situation_dict = get_read_situation_dict(remapped_diff_gene)
        all_drp_by_situation = get_all_drp_by_situation(read_situation_dict)
        strong_evidences, retrocopy_clues, alteration_clues, homology_check = clues_organisation(all_drp_by_situation, all_read_infos)
        
        strong_evidences_dict.update({exp_name:strong_evidences})
        
        retrocopy_clues_dict.update({exp_name:retrocopy_clues})
        
        alteration_clues_dict.update({exp_name:alteration_clues})
        
        homology_check_dict.update({exp_name:homology_check})
            


    return all_counts, strong_evidences_dict, retrocopy_clues_dict, alteration_clues_dict, homology_check_dict


def get_violin_plots(data_N, data_T, out_name):
    plt.figure(figsize=(3, 4))

    # Fond blanc
    sns.set(style="white")

    fig, ax = plt.subplots()

    # Violin plot 
    ax = sns.violinplot(data=[data_N, data_T], palette=['#C8fcc8', 'gold'], cut=0, inner='quartiles')

    # Swarplot
    sns.swarmplot(data=[data_N, data_T], palette=['grey', 'grey'], s=4)

    # Contour du violin (vert foncé)
    for vp in ax.collections:
        vp.set_edgecolor('#076842')

    ax.set_xticklabels(['% DRPs normal tissue', '% DRPs tumoral tissue'])
    ax.yaxis.set_tick_params(color='black', labelcolor='black')

    # Désactiver le cadre autour de la figure
    plt.box(False)

    # Extraire les valeurs de la médiane, du 1er quartile et du 3ème quartile avec NumPy
    statistics_N = np.percentile(data_N, [25, 50, 75])

    q1_N = statistics_N[0]  # 1er quartile
    median_N = statistics_N[1]  # Médiane
    q3_N = statistics_N[2]  # 3ème quartile

    print("Médiane_N:", median_N)
    print("1er quartile_N:", q1_N)
    print("3ème quartile_N:", q3_N)

    statistics_T = np.percentile(data_T, [25, 50, 75])

    q1_T = statistics_T[0]  # 1er quartile
    median_T = statistics_T[1]  # Médiane
    q3_T = statistics_T[2]  # 3ème quartile

    print("Médiane_N:", median_T)
    print("1er quartile_N:", q1_T)
    print("3ème quartile_N:", q3_T)

    stat_test = scipy.stats.mannwhitneyu(data_N, data_T)
    print('MannWhitney results', stat_test)

    plt.title('Discordant Read Pair Proportion', fontsize=15)
    plt.savefig(out_name)

# NEW REMAPPED NOTEBOOK

def best_score_diff_gene(drp):
    r1, r2 = drp
    score_dict_1 = {'score_orig' : r1['score_orig'], 'score_mate' : r1['score_mate'], 'score_exon' : r1['score_exon'], 'score_exon_mate' : r1['score_exon_mate']}
    score_dict_2 = {'score_orig' : r2['score_orig'], 'score_mate' : r2['score_mate'], 'score_exon' : r2['score_exon'], 'score_exon_mate' : r2['score_exon_mate']}

    best_scores_1 = [score_type  for score_type in score_dict_1 if score_dict_1[score_type] == max(score_dict_1.values())]
    best_scores_2 = [score_type  for score_type in score_dict_2 if score_dict_2[score_type] == max(score_dict_2.values())]

    return best_scores_1, best_scores_2, score_dict_1, score_dict_2

def get_read_situation_dict(remapped_diff_gene):
    read_situation_dict = {}
    for drp in remapped_diff_gene:
        best_scores_1, best_scores_2, score_dict_1, score_dict_2 = best_score_diff_gene(drp)
        read_situation_dict.update({drp[0]['name']:(best_scores_1, best_scores_2)})
    return read_situation_dict

def get_all_drp_by_situation(read_situation_dict):
    all_drp_by_situation = {}
    for drp, (situation_1, situation_2) in read_situation_dict.items():
        
        real_situation_1 = get_real_situation(situation_1)
        real_situation_2 = get_real_situation(situation_2)
        
        if (real_situation_1, real_situation_2) in all_drp_by_situation:
            all_drp_by_situation[(real_situation_1, real_situation_2)].append(drp)
        elif (real_situation_1, real_situation_2) not in all_drp_by_situation:
            all_drp_by_situation.update({(real_situation_1, real_situation_2):[drp]})
    return all_drp_by_situation

def get_real_situation(situation):
    if situation == ['score_orig'] or situation == ['score_orig', 'score_exon']:
        real_situation = 'Original gene'
    elif situation == ['score_mate'] or situation == ['score_mate', 'score_exon_mate']:
        real_situation = 'Mate gene'
    elif situation == ['score_exon']:
        real_situation = 'Gene orig Retrocopy'
    elif situation == ['score_exon_mate']:
        real_situation = 'Gene mate Retrocopy'
    else:
        real_situation = 'Homology'
    return real_situation

def clues_organisation(all_drp_by_situation, all_read_infos):
    retrocopy_clues = {}
    strong_evidences = {}
    alteration_clues = {}
    homology_check = {}
    for situation_key in all_drp_by_situation.keys():
        if situation_key == ('Gene orig Retrocopy', 'Gene mate Retrocopy') or situation_key == ('Gene mate Retrocopy', 'Gene orig Retrocopy'):
            retrocopied_gene_ind = situation_key.index('Gene orig Retrocopy') # This index is to see if its read A or B the retrocopied gene
            for drp_name in all_drp_by_situation[situation_key]:
                retrocopied_gene = all_read_infos[drp_name][retrocopied_gene_ind]['gene_orig']
                if retrocopied_gene in strong_evidences:
                    strong_evidences[retrocopied_gene].append(drp_name)
                elif retrocopied_gene not in strong_evidences:
                    strong_evidences[retrocopied_gene] = [drp_name]

        elif 'Homology' in situation_key:
            for drp_name in all_drp_by_situation[situation_key]:
                gene_1 = all_read_infos[drp_name][0]['gene_orig']
                gene_2 = all_read_infos[drp_name][0]['gene_mate']
                if (gene_1, gene_2) in homology_check:
                    homology_check[(gene_1, gene_2)].append(drp_name)
                elif (gene_1, gene_2) not in homology_check:
                    homology_check[(gene_1, gene_2)] = [drp_name]
        
        elif 'Gene orig Retrocopy' in situation_key:
            retrocopied_gene_ind = situation_key.index('Gene orig Retrocopy')
            for drp_name in all_drp_by_situation[situation_key]:
                retrocopied_gene = all_read_infos[drp_name][retrocopied_gene_ind]['gene_orig']
                host_gene = all_read_infos[drp_name][retrocopied_gene_ind]['gene_mate']
                if (retrocopied_gene, host_gene) in retrocopy_clues:
                    retrocopy_clues[(retrocopied_gene, host_gene)].append(drp_name)
                elif (retrocopied_gene, host_gene) not in retrocopy_clues:
                    retrocopy_clues[(retrocopied_gene, host_gene)] = [drp_name]
    
        
        elif 'Gene mate Retrocopy' in situation_key:
            retrocopied_gene_ind = situation_key.index('Gene mate Retrocopy')
            for drp_name in all_drp_by_situation[situation_key]:
                retrocopied_gene = all_read_infos[drp_name][retrocopied_gene_ind]['gene_mate']
                host_gene = all_read_infos[drp_name][retrocopied_gene_ind]['gene_orig']
                if (retrocopied_gene, host_gene) in retrocopy_clues:
                    retrocopy_clues[(retrocopied_gene, host_gene)].append(drp_name)
                elif (retrocopied_gene, host_gene) not in retrocopy_clues:
                    retrocopy_clues[(retrocopied_gene, host_gene)] = [drp_name]

        else : 
            if 'Original gene' not in situation_key:
                #print(situation_key)
                loc_index = situation_key.index('Mate gene')
                for drp_name in all_drp_by_situation[situation_key]:
                    gene_1 = all_read_infos[drp_name][loc_index]['gene_mate']
                    gene_2 = all_read_infos[drp_name][loc_index]['gene_orig']
                    if (gene_1, gene_2) in alteration_clues:
                        alteration_clues[(gene_1, gene_2)].append(drp_name)
                    elif (gene_1, gene_2) not in alteration_clues:
                        alteration_clues[(gene_1, gene_2)] = [drp_name]
            else:
                #print(situation_key)
                loc_index = situation_key.index('Original gene')
                for drp_name in all_drp_by_situation[situation_key]:
                    gene_1 = all_read_infos[drp_name][loc_index]['gene_orig']
                    gene_2 = all_read_infos[drp_name][loc_index]['gene_mate']
                    if (gene_1, gene_2) in alteration_clues:
                        alteration_clues[(gene_1, gene_2)].append(drp_name)
                    elif (gene_1, gene_2) not in alteration_clues:
                        alteration_clues[(gene_1, gene_2)] = [drp_name]

    
    return strong_evidences, retrocopy_clues, alteration_clues, homology_check

def get_evidences_fig(strong_evidences_compilation, out_name):
    strong_ev_fig = {}
    clues_fig = {}
    alt_fig = {}
    for ex, gene_evidences in strong_evidences_compilation.items():
        strong_count = 0
        clues_count = 0
        alt_count = 0
        for gene, evidences in gene_evidences.items():
            strong_count += len(evidences['strong'])
            for gene_duo, clues_list in evidences['clues']:
                clues_count += len(clues_list)
            for gene_duo, alt_list in evidences['alteration']:
                alt_count += len(alt_list)

        strong_ev_fig.update({ex:strong_count})
        clues_fig.update({ex:clues_count})
        alt_fig.update({ex:alt_count})

    strong_ev_fig_sorted = dict(sorted(strong_ev_fig.items(), key=lambda kv: kv[0]))
    clues_fig_sorted = dict(sorted(clues_fig.items(), key=lambda kv: kv[0]))
    alt_fig_sorted = dict(sorted(alt_fig.items(), key=lambda kv: kv[0]))

    fig, ax = plt.subplots(nrows=1, ncols=3)

    exome = [exp_name.split('_')[0] for exp_name in strong_ev_fig_sorted.keys()]
    y_pos = np.arange(len(exome))

    ax[0].barh(y_pos, list(strong_ev_fig_sorted.values()), align='center')
    ax[0].set_yticks(y_pos, labels=exome)
    ax[0].set_xlabel('Number of strong retrocopy evidences', fontsize = 15)

    ax[1].barh(y_pos, list(clues_fig_sorted.values()), align='center')
    ax[1].set_yticks(y_pos, labels=exome)
    ax[1].set_xlabel('Number of retrocopy clues', fontsize = 15)

    ax[2].barh(y_pos, list(alt_fig_sorted.values()), align='center')
    ax[2].set_yticks(y_pos, labels=exome)
    ax[2].set_xlabel('Number of alteration clues', fontsize = 15)

    fig.set_size_inches(15, 10.5)
    plt.savefig(out_name)
#NOTEBOOK END 

def get_TDG_evidences_fig(strong_evidences_compilation, out_name):
    strong_ev_fig = {}
    clues_fig = {}
    alt_fig = {}
    for ex, gene_evidences in strong_evidences_compilation.items():
        strong_count = 0
        clues_count = 0
        alt_count = 0
        for gene, evidences in gene_evidences.items():
            if gene == 'ENSG00000139372':
                strong_count += len(evidences['strong'])
                for gene_duo, clues_list in evidences['clues']:
                    clues_count += len(clues_list)
                for gene_duo, alt_list in evidences['alteration']:
                    alt_count += len(alt_list)

        strong_ev_fig.update({ex:strong_count})
        clues_fig.update({ex:clues_count})
        alt_fig.update({ex:alt_count})

    strong_ev_fig_sorted = dict(sorted(strong_ev_fig.items(), key=lambda kv: kv[0]))
    clues_fig_sorted = dict(sorted(clues_fig.items(), key=lambda kv: kv[0]))
    alt_fig_sorted = dict(sorted(alt_fig.items(), key=lambda kv: kv[0]))

    fig, ax = plt.subplots(nrows=1, ncols=3)

    exome = [exp_name.split('_')[0] for exp_name in strong_ev_fig_sorted.keys()]
    y_pos = np.arange(len(exome))

    ax[0].barh(y_pos, list(strong_ev_fig_sorted.values()), align='center')
    ax[0].set_yticks(y_pos, labels=exome)
    ax[0].set_xlabel('Number of strong retrocopy evidences', fontsize = 15)

    ax[1].barh(y_pos, list(clues_fig_sorted.values()), align='center')
    ax[1].set_yticks(y_pos, labels=exome)
    ax[1].set_xlabel('Number of retrocopy clues', fontsize = 15)

    ax[2].barh(y_pos, list(alt_fig_sorted.values()), align='center')
    ax[2].set_yticks(y_pos, labels=exome)
    ax[2].set_xlabel('Number of alteration clues', fontsize = 15)

    fig.set_size_inches(15, 10.5)
    plt.savefig(out_name)

def get_noTDG_evidences_fig(strong_evidences_compilation, out_name):
    strong_ev_fig = {}
    clues_fig = {}
    alt_fig = {}
    for ex, gene_evidences in strong_evidences_compilation.items():
        strong_count = 0
        clues_count = 0
        alt_count = 0
        for gene, evidences in gene_evidences.items():
            if gene != 'ENSG00000139372':
                strong_count += len(evidences['strong'])
                for gene_duo, clues_list in evidences['clues']:
                    clues_count += len(clues_list)
                for gene_duo, alt_list in evidences['alteration']:
                    alt_count += len(alt_list)

        strong_ev_fig.update({ex:strong_count})
        clues_fig.update({ex:clues_count})
        alt_fig.update({ex:alt_count})

    strong_ev_fig_sorted = dict(sorted(strong_ev_fig.items(), key=lambda kv: kv[0]))
    clues_fig_sorted = dict(sorted(clues_fig.items(), key=lambda kv: kv[0]))
    alt_fig_sorted = dict(sorted(alt_fig.items(), key=lambda kv: kv[0]))

    fig, ax = plt.subplots(nrows=1, ncols=3)

    exome = [exp_name.split('_')[0] for exp_name in strong_ev_fig_sorted.keys()]
    y_pos = np.arange(len(exome))

    ax[0].barh(y_pos, list(strong_ev_fig_sorted.values()), align='center')
    ax[0].set_yticks(y_pos, labels=exome)
    ax[0].set_xlabel('Number of strong retrocopy evidences', fontsize = 15)

    ax[1].barh(y_pos, list(clues_fig_sorted.values()), align='center')
    ax[1].set_yticks(y_pos, labels=exome)
    ax[1].set_xlabel('Number of retrocopy clues', fontsize = 15)

    ax[2].barh(y_pos, list(alt_fig_sorted.values()), align='center')
    ax[2].set_yticks(y_pos, labels=exome)
    ax[2].set_xlabel('Number of alteration clues', fontsize = 15)

    fig.set_size_inches(15, 10.5)
    plt.savefig(out_name)

def run_all(remapped_diff_list, input_path):

    print('starting ...', flush=True)
    # Get all_counts
    all_counts, strong_evidences_dict, retrocopy_clues_dict, alteration_clues_dict, homology_check_dict = get_all_counts(remapped_diff_list, input_path)

    print('All dict done', flush=True)

    fig_path = f'{input_path}figures_WES/'
    if os.path.exists(fig_path) == False:
        os.mkdir(fig_path)

    pickle.dump(strong_evidences_dict, open(f'{fig_path}/strong_evidences_dict.pkl', 'wb'))
    pickle.dump(retrocopy_clues_dict, open(f'{fig_path}/retrocopy_clues_dict.pkl', 'wb'))
    pickle.dump(alteration_clues_dict, open(f'{fig_path}/alteration_clues_dict.pkl', 'wb'))
    pickle.dump(homology_check_dict, open(f'{fig_path}/homology_check_dict.pkl', 'wb'))

    print('All dict dump', flush=True)

    # Get exp_list for which we have N and T
    T_patient = []
    N_patient = []
    for exp_name in all_counts.keys():
        patient_id = exp_name.split('_')[0]
        if 'T_WES' in exp_name:
            T_patient.append(patient_id)
        if 'N_WES' in exp_name:
            N_patient.append(patient_id)
    both_patient = list(set(T_patient) & set(N_patient))

    print(f'Both_list : {both_patient}, {len(both_patient)}', flush=True)

    done_exp = []
    for exp_name in all_counts.keys():
        patient_id = exp_name.split('_')[0]
        if patient_id in both_patient:
            done_exp.append(exp_name)

    print(f'Done exp : {done_exp}, {len(done_exp)}', flush=True)

    # Get violin plots
    tumor_data = [count_dict['disc_reads']/(count_dict['disc_reads'] + count_dict['good_reads'])*100 for exp_name, count_dict in all_counts.items() if 'T_WES' in exp_name and exp_name in done_exp]
    normal_data = [count_dict['disc_reads']/(count_dict['disc_reads'] + count_dict['good_reads'])*100 for exp_name, count_dict in all_counts.items() if 'N_WES' in exp_name and exp_name in done_exp]

    print(f'tumor_data : {tumor_data}, {len(tumor_data)}', flush=True)
    print(f'normal_data : {normal_data}, {len(normal_data)}', flush=True)

    out_name = f'{fig_path}DRP%_violin.svg'
    get_violin_plots(normal_data, tumor_data, out_name)

    print('Violin plot done', flush=True)

    # Gets evidences compilation
    print('evidences_compilation starting...', flush=True)
    strong_evidences_compilation_T = {}
    strong_evidences_compilation_N = {}

    for exp_name in strong_evidences_dict:
        if 'T_WES' in exp_name:
            strong_evidences_compilation_T.update({exp_name:{}})
            for gene, drps in strong_evidences_dict[exp_name].items():
                strong_evidences_compilation_T[exp_name].update({gene:{'strong':drps, 'clues':[], 'alteration':[]}})
                for gene_comb, clues_drps in retrocopy_clues_dict[exp_name].items():
                    if gene in gene_comb:
                        strong_evidences_compilation_T[exp_name][gene]['clues'].append((gene_comb, clues_drps))
                for gene_comb, clues_drps in alteration_clues_dict[exp_name].items():
                    if gene in gene_comb:
                        strong_evidences_compilation_T[exp_name][gene]['alteration'].append((gene_comb, clues_drps))
        
        if 'N_WES' in exp_name:
            strong_evidences_compilation_N.update({exp_name:{}})
            for gene, drps in strong_evidences_dict[exp_name].items():
                strong_evidences_compilation_N[exp_name].update({gene:{'strong':drps, 'clues':[], 'alteration':[]}})
                for gene_comb, clues_drps in retrocopy_clues_dict[exp_name].items():
                    if gene in gene_comb:
                        strong_evidences_compilation_N[exp_name][gene]['clues'].append((gene_comb, clues_drps))
                for gene_comb, clues_drps in alteration_clues_dict[exp_name].items():
                    if gene in gene_comb:
                        strong_evidences_compilation_N[exp_name][gene]['alteration'].append((gene_comb, clues_drps))


    pickle.dump(strong_evidences_compilation_T, open(f'{fig_path}/strong_evidences_compilation_T.pkl', 'wb'))
    pickle.dump(strong_evidences_compilation_N, open(f'{fig_path}/strong_evidences_compilation_N.pkl', 'wb'))

    print('evidences_compilation done', flush=True)
    # Get evidence figures

    fig_name = f'{fig_path}tumor_evidences.svg'
    get_evidences_fig(strong_evidences_compilation_T, fig_name)

    fig_name = f'{fig_path}normal_evidences.svg'
    get_evidences_fig(strong_evidences_compilation_N, fig_name)

    print('Evidences fig done', flush=True)

    fig_name = f'{fig_path}tumor_evidences_TDG.svg'
    get_TDG_evidences_fig(strong_evidences_compilation_T, fig_name)

    fig_name = f'{fig_path}normal_evidences_TDG.svg'
    get_TDG_evidences_fig(strong_evidences_compilation_N, fig_name)

    fig_name = f'{fig_path}tumor_evidences_noTDG.svg'
    get_noTDG_evidences_fig(strong_evidences_compilation_T, fig_name)

    fig_name = f'{fig_path}normal_evidences_noTDG.svg'
    get_noTDG_evidences_fig(strong_evidences_compilation_N, fig_name)

    print('Evidences fig done', flush=True)


def main():
    remapped = snakemake.input.remapped_drps
    same_gene_path = snakemake.input.same_gene
    reads_stats_path = snakemake.input.reads_stats
    out_comp = snakemake.output.compilation

    # Normaliser en liste
    paths = [remapped] if isinstance(remapped, str) else remapped

    # Lire directement les fichiers fournis par Snakemake
    reads_stats = pickle.load(open(reads_stats_path, 'rb'))
    remapped_same_gene = pickle.load(open(same_gene_path, 'rb'))
    remapped_diff_gene = pickle.load(open(remapped, 'rb'))

    # En déduire exp_name (si utile pour noms ou figures)
    exp_path = os.path.dirname(os.path.dirname(remapped))
    exp_name = os.path.basename(os.path.dirname(remapped))

    # Remplacer cette fonction si elle reste utile
    all_counts = {
        exp_name: {
            'good_reads': reads_stats['good_reads'],
            'disc_reads': reads_stats['disc_reads'],
            'remapped_same_gene': len(remapped_same_gene),
            'remapped_diff_gene': len(remapped_diff_gene)
        }
    }

    all_read_infos = {drp[0]['name']: drp for drp in remapped_diff_gene}
    read_situation_dict = get_read_situation_dict(remapped_diff_gene)
    all_drp_by_situation = get_all_drp_by_situation(read_situation_dict)
    strong_ev, retrocopy_clues, alteration_clues, homol = clues_organisation(all_drp_by_situation, all_read_infos)

    # Compilation T uniquement (tu pourrais factoriser le reste comme avant)
    strong_evidences_compilation_T = {
        exp_name: {
            gene: {
                'strong': drps,
                'clues': [(k, v) for k, v in retrocopy_clues.items() if gene in k],
                'alteration': [(k, v) for k, v in alteration_clues.items() if gene in k]
            }
            for gene, drps in strong_ev.items()
        }
    }

    pickle.dump(strong_evidences_compilation_T, open(out_comp, "wb"))


if __name__ == "__main__":
    main()