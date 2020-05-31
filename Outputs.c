#include "Outputs.h"
#include "Constants.h"
#include "GlobalVariables.h"

void NetworksDOT (FILE *file_netw, network *networks , int *dominanceScore, int g, double *target1, double *target2, double *target3, double *Ff, const char *argvCOND){
    
    int i,j,k;
    int netw_count;
    char string_netw[500];
    netw_count=0;
    for (i=0;i<SIZEPOP;i++)
    {
        if (dominanceScore[i]==0)
        {
            //Storage of Pareto-optimal network's topology in .dot format

            sprintf(string_netw,"./%03d_%s_dr%.2f_mrr%.2f_a%d/DOT/Gen%d(Netw%d).dot",REPLICA,argvCOND, damage_rate, minimum_regeneration_rate, TIME_AGEING,g,netw_count);
            file_netw = fopen(string_netw,"w");
            
            fprintf(file_netw, "digraph G{\nrankdir=LR\nsplines=line\nnodesep=.05;\nsubgraph cluster_0{\n\ncolor=white;\nnode [style=solid,color=blue4, shape=circle];\ni1 i2 i3;\nlabel = \"Input\";\n}\n node [label=\"\"];\nsubgraph cluster_1{\ncolor=white;\nnode [style=solid,color=red2, shape=circle];");
            for (j=0;j<networks[i].hid1;j++)
                if (networks[i].vec_aux_truly_feas_hid1[j]==1)
                    fprintf(file_netw, "a%d2 ",j+1);
            fprintf(file_netw, ";\n");
            fprintf(file_netw, "label = \"Hidden neurons(layer 1)\";\n}\nsubgraph cluster_2{\ncolor=white;\nnode [style=solid,color=red2, shape=circle];\n");
            
            for (j=0;j<networks[i].hid2;j++)
                if (networks[i].vec_aux_truly_feas_hid2[j]==1)
                    fprintf(file_netw, "a%d3 ",j+1);
            fprintf(file_netw, ";\n");
            fprintf(file_netw, "label = \"Hidden neurons(layer 2)\";\n}\nsubgraph cluster_3{\ncolor=white;\nnode [style=solid,color=seagreen2, shape=circle];\nO;\nlabel=\"Output\";\n}\n");
            
            for (j=0;j<INPUT;j++)
                for(k=0;k<networks[i].hid1;k++)
                {
                    if (networks[i].vec_aux_truly_feas_hid1[k]==1)
                    {
                        if(networks[i].inhid1[j*networks[i].hid1+k] == 1)
                            fprintf(file_netw,"edge[color=black];\ni%d->a%d2;\n",j+1,k+1);
                        if(networks[i].inhid1[j*networks[i].hid1+k] == -1)
                            fprintf(file_netw,"edge[color=red];\ni%d->a%d2;\n",j+1,k+1);
                    }
                }
            for (j=0;j<networks[i].hid1;j++)
                if (networks[i].vec_aux_truly_feas_hid1[j]==1)
                {
                    for(k=0;k<networks[i].hid2;k++)
                    {
                        if (networks[i].vec_aux_truly_feas_hid2[k]==1)
                        {
                            if(networks[i].hid1hid2[j*networks[i].hid2+k] == 1)
                                fprintf(file_netw,"edge[color=black];\na%d2->a%d3;\n",j+1,k+1);
                            if(networks[i].hid1hid2[j*networks[i].hid2+k] == -1)
                                fprintf(file_netw,"edge[color=red];\na%d2->a%d3;\n",j+1,k+1);
                        }
                    }
                }
            
            for (j=0;j<networks[i].hid2;j++)
                if (networks[i].vec_aux_truly_feas_hid2[j]==1)
                {
                    for(k=0;k<OUTPUT;k++)
                    {
                        if(networks[i].hid2out[j] == 1)
                            fprintf(file_netw,"edge[color=black];\na%d3->O;\n",j+1);
                        if(networks[i].hid2out[j] == -1)
                            fprintf(file_netw,"edge[color=red];\na%d3->O;\n",j+1);
                    }
                    
                }
            
            fprintf(file_netw, "}\n\n /*Target1 = %lf, Target2 = %lf, Target3 = %lf, Ff = %lf*/\n",target1[i], target2[i],target3[i],Ff[i]);
            
            netw_count++;
            fclose(file_netw);
        }
        
    }
};

void TargetsData (FILE *pareto_data, FILE *nonpareto_data, int g, int *dominanceScore, double *target1, double *target2, double *target3, double *Ff, const char *argvCOND){
    
    char string_pareto[1000], string_nonpareto[1000];
    int i;
    sprintf(string_pareto,"./%03d_%s_dr%.2f_mrr%.2f_a%d/ParetoNonPareto/Gen%d.dat",REPLICA, argvCOND,damage_rate, minimum_regeneration_rate, TIME_AGEING,g);
    sprintf(string_nonpareto, "./%03d_%s_dr%.2f_mrr%.2f_a%d/ParetoNonPareto/Gen%d_nonpareto.dat",REPLICA, argvCOND,damage_rate, minimum_regeneration_rate, TIME_AGEING,g);
    
    pareto_data = fopen(string_pareto,"w");
    nonpareto_data = fopen(string_nonpareto,"w");
    for (i=0;i<SIZEPOP;i++)
    {
        if (dominanceScore[i]==0)
            fprintf(pareto_data, "%lf %lf %lf %lf\n", target1[i],target2[i],target3[i],Ff[i]);
        else
            fprintf(nonpareto_data, "%lf %lf %lf %lf\n", target1[i],target2[i],target3[i],Ff[i]);
    }

    fclose(pareto_data);
    fclose(nonpareto_data);
};

void ConnectivityNetworksStorage(FILE *connectivity, network *networks, int g, int *domscore, const char *argvCOND){
    int i,j,k,counter;
    char string[500];
    sprintf(string, "./%03d_%s_dr%.2f_mrr%.2f_a%d/NetworksData/Gen%d_ConnectivityStored.dat",REPLICA, argvCOND, damage_rate, minimum_regeneration_rate, TIME_AGEING,g);
    connectivity = fopen(string,"w");
    counter=0;
    for (i=0;i<SIZEPOP;i++)
    {
        if (domscore[i]==0)
        {
            fprintf(connectivity, "N%d\n",counter);
            fprintf(connectivity, "%d\n%d\n",networks[i].truly_feas_hid1,networks[i].truly_feas_hid2);
            for (j=0;j<INPUT;j++)
                for(k=0;k<networks[i].hid1;k++)
                    if (networks[i].vec_aux_truly_feas_hid1[k]==1)
                        fprintf(connectivity,"%d",networks[i].inhid1[j*networks[i].hid1+k]);
            fprintf(connectivity,"\n");
            
            for (j=0;j<networks[i].hid1;j++)
                if (networks[i].vec_aux_truly_feas_hid1[j]==1)
                    for(k=0;k<networks[i].hid2;k++)
                        if (networks[i].vec_aux_truly_feas_hid2[k]==1)
                            fprintf(connectivity,"%d",networks[i].hid1hid2[j*networks[i].hid2+k]);
            fprintf(connectivity,"\n");
            
            for (j=0;j<networks[i].hid2;j++)
                if (networks[i].vec_aux_truly_feas_hid2[j]==1)
                    for(k=0;k<OUTPUT;k++)
                        fprintf(connectivity,"%d", networks[i].hid2out[j*OUTPUT+k]);
            fprintf(connectivity,"\n");
            counter++;
        }
    }
    fclose(connectivity);

};

