#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include </home/choi/mt19937-64.c>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

int temp;
char file_path[500];
char file_path_nw[500];
char file_path_nw2[500];
char header_file[500];
FILE *fp;
FILE *fp2;
FILE *fpop;
FILE *fpnw;
FILE *fpnw2;

void decimalToBinary(unsigned long long n, int binary[],int num_genes_tot) {
	for (int i = num_genes_tot-1; i >= 0; i--) {
		binary[i] = n % 2;
		n = n / 2;
	}
}

void ReadNW(char*file_path,int*nl,int**ll) {
	int i,j;
	fp=fopen(file_path,"r"); 
	if(fp!=NULL){
		//printf("%s\n",file_path);
		while(fscanf(fp,"%d %d",&i,&j)!=EOF) { 
			//printf("%d %d\n",i,j);
			ll[i][nl[i]]=j;
			nl[i]++;

			ll[i]=(int*)realloc(ll[i],sizeof(int)*(nl[i]+1));
		}
		fclose(fp);
	}
	else printf("no nw files\n");
}

void Read_num_genes(char*file_path,int numbers[]) {
	int i;
	numbers[0]=0; //genes have truth table
	numbers[1]=0; //total num genes
	char key[100];
	fp=fopen(file_path,"r"); 
	if(fp!=NULL){
		//printf("%s\n",file_path);
		while(fscanf(fp,"%d\t%s",&i,key)!=EOF) { 
			if(i<1000000) { numbers[0]++; }
			numbers[1]++;
		}
		fclose(fp);
	}
	else printf("no species key file\n");

}
void Read_matrix_header(char*file_path,int num_genes_tot,char key[][500]) {
	int i;
	char key0[4000];
	char *token;

	fp=fopen(file_path,"r"); 
	if(fp!=NULL){
		printf("%s\n",file_path);

		fscanf(fp,"%s",key0);
		token=strtok(key0, ",");
		i=0;
		while (token != NULL) {
			strcpy(key[i], token);
			token = strtok(NULL, ","); // 다음 토큰으로 이동
			i++;
		}
		fclose(fp);
	}
	else printf("no matrix file, header\n");

}


void Truth_table_scan(char*file_path,int num_genes_tot,char gene[],double bias_input[3]) {
	int i;
	char file[500];
	char file2[500];
	char line[500];
	char *token;
	char key_tt[num_genes_tot][500];
	int num_input;
	double bias_data;

	sprintf(file,"%s/%s_nospace.csv",file_path,gene);
	fp=fopen(file,"r"); 
	if(fp!=NULL){
		fscanf(fp,"%s",line); 
		//printf("%s\n",line);

		token=strtok(line, ",");
		i=0;
		while (token != NULL) {
			strcpy(key_tt[i], token);
			token = strtok(NULL, ","); // 다음 토큰으로 이동
			i++;
		}
		num_input=i-1;

		int num_outlier=0;
		for(i=0;i<num_input;i++) {
			sprintf(file2,"%s/%s_nospace.csv",file_path,key_tt[i]);
			fp2=fopen(file2,"r"); 
			if(fp2==NULL) { 
				num_outlier++;
				//printf("nofile\n");
			}
			else { 
				fclose(fp2);  
				//printf("%s\n",file2); 
			}
		}
	

		char last_char;
		int len;
		int bias_count=0;
		for(i=0;i<pow(2,num_input);i++) {
			fscanf(fp,"%s",line); 

			len=strlen(line);
			last_char=line[len - 1]; 
			if(last_char=='1') bias_count+=1;
			//printf("%s\n",line);
		}
		//if(bias_count==0) printf("%s: %g\n",gene,bias_count/pow(2,num_input));
		bias_data=bias_count/pow(2,num_input);

		fclose(fp);
		bias_input[0]=bias_data;
		bias_input[1]=num_input;
		bias_input[2]=num_outlier;
	}
	else {
		bias_input[0]=0;
		bias_input[1]=0;
		bias_input[2]=0;
	}

}

int main(int argc, char *argv[]) {  
	char filename[400];
	char origin_filename[400];
	int nn; 
	sprintf(filename,"numnodes-filenames");
	fpop=fopen(filename,"r"); 

	sprintf(file_path_nw,"./bias");
	fpnw=fopen(file_path_nw,"w"); 
	fprintf(fpnw,"#filename bias avg_num_input avg_num_outlier num_genes\n");

	sprintf(file_path_nw2,"./bias_by_genes");
	fpnw2=fopen(file_path_nw2,"w"); 
	fprintf(fpnw2,"#filename gene_key bias num_input num_outlier\n");
	
	while (fscanf(fpop,"%d %s",&nn,origin_filename)!= EOF) {
		//printf("%s %d\n",origin_filename,nn);

		int numbers[2];
		numbers[0]=0; numbers[1]=0;
		sprintf(file_path,"./DATA_truth/%s_TT/SPECIES_KEY.csv",origin_filename);
		Read_num_genes(file_path,numbers);
		int num_genes_tot=numbers[1];

		if (num_genes_tot>0) {
			char key[num_genes_tot][500];
			sprintf(file_path,"./DATA_matrix/%s_MATRIX_pm.csv",origin_filename);
			Read_matrix_header(file_path,num_genes_tot,key);

			double bias_input[3];
			double avg_bias=0;
			double bias;
			double avg_num_input=0;
			double avg_num_outlier=0;
			for(int i=0;i<num_genes_tot;i++) {

				sprintf(file_path,"./DATA_truth/%s_TT",origin_filename);
				Truth_table_scan(file_path,num_genes_tot,key[i],bias_input);
				if(bias_input[0]>0.5) bias=1-bias_input[0];
				else bias=bias_input[0];

				fprintf(fpnw2,"%s %s %g %g %g\n",origin_filename,key[i],bias,bias_input[1],bias_input[2]); 
				avg_bias+=bias;
				avg_num_input+=bias_input[1];
				avg_num_outlier+=bias_input[2];
			}
	
			fprintf(fpnw,"%s %g %g %g %d\n",origin_filename,avg_bias/numbers[0],avg_num_input/numbers[0],avg_num_outlier/numbers[0],numbers[1]); 

		}
	}
	fclose(fpnw);
	fclose(fpnw2);
	fclose(fpop);

	return 0;
}


