#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <time.h>
#include </home/choi/mt19937-64.c>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

int temp;
char file_path[400];
char file_path_nw[400];
char header_file[400];
FILE *fp;
FILE *fp2;
FILE *fpnw;

void decimalToBinary(unsigned long long n, int binary[],int num_genes_tot,int input_check[]) {
	for (int i = num_genes_tot-1; i >= 0; i--) {
		if(input_check[i]==1) {
			binary[i] = n % 2;
			n = n / 2;
		}
		else binary[i]=0;
		//printf("%d %d %d\n",i,input_check[i],binary[i]);
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

void Read_matrix(char*file_path,int num_genes_tot,char key[][500],int pm[][num_genes_tot]) {
	int i,j;
	char key0[500];
	char *token;

	fp=fopen(file_path,"r"); 
	if(fp!=NULL){
		//printf("%s\n",file_path);

		fscanf(fp,"%s",key0);
		token=strtok(key0, ",");
		i=0;
		while (token != NULL) {
			strcpy(key[i], token);
			token = strtok(NULL, ","); // 다음 토큰으로 이동
			i++;
		}

		int sgn;
		for(i=0;i<num_genes_tot;i++) { 
			for(j=0;j<num_genes_tot-1;j++) { fscanf(fp,"%d,",&sgn); pm[i][j]=abs(sgn); }
			fscanf(fp,"%d",&sgn); pm[i][j]=abs(sgn); 
		}
		fclose(fp);
	}
	else printf("no matrix file\n");

}

void Read_matrix_header(char*origin_filename,int num_genes_tot,char key[][500],int input_check[]) {
	int i;
	char key0[500];
	char *token;
	char file_path1[400];
	char file_path2[400];
	sprintf(file_path1,"./DATA_matrix/%s_MATRIX_pm.csv",origin_filename);
	fp=fopen(file_path1,"r"); 
	if(fp!=NULL){
		//printf("%s\n",file_path);

		fscanf(fp,"%s",key0);
		token=strtok(key0, ",");
		i=0;
		while (token != NULL) {
			strcpy(key[i], token);

			sprintf(file_path2,"./DATA_truth/%s_TT/%s_nospace.csv",origin_filename,key[i]);
			fp2=fopen(file_path2,"r"); 
			if(fp2!=NULL) { 
				input_check[i]=1; fclose(fp2); 
				//printf("%d %s %s\n",i,key[i],file_path2);
			}
			else input_check[i]=0;

			token = strtok(NULL, ","); // 다음 토큰으로 이동
			i++;
		}
		fclose(fp);
	}
	else printf("no matrix file, header\n");

}


char Truth_table_ex(char*file_path,int num_genes_tot,unsigned long long state,int binary[],char key[][500],char gene[]) {
	int i,j;
	char file[500];
	char line[500];
	char *token;
	char key_tt[num_genes_tot][500];
	int num_input;

	sprintf(file,"%s/%s_nospace.csv",file_path,gene);
	fp=fopen(file,"r"); 
	if(fp!=NULL){

		fscanf(fp,"%s",line); 
		//printf("header: %s\n",line);
		//printf("--- key ------------  gene:%s\n",gene);
		token=strtok(line, ",");
		i=0;
		while (token != NULL) {
			strcpy(key_tt[i], token);
			//printf("%s\n",key_tt[i]);
			token = strtok(NULL, ","); // 다음 토큰으로 이동
			i++;
		}
		num_input=i-1;
		//printf("--- key has read ---\n");

		int binary_tt[num_input];
		for(i=0;i<num_input;i++) for(j=0;j<num_genes_tot;j++)  { 
			if(strcmp(key[j],key_tt[i])==0) binary_tt[i]=binary[j];
		}
		//for(i=0;i<num_input;i++) printf("%d",binary_tt[i]);

		int tt_deci=0;
		for(i=0;i<num_input;i++) { tt_deci+=binary_tt[num_input-1-i]*pow(2,i); }
		//printf(" = %d(demimal)\n",tt_deci);
	
		for(i=0;i<tt_deci+1;i++) fscanf(fp,"%s",line); 
		//printf("line: %s\n",line);

		int len = strlen(line);
		char last_char = line[len - 1]; 
		//printf("%c\n",last_char);
	
		fclose(fp);

		return last_char;

	}
	else { 
//		printf("no species key file\n"); 
		return '\0'; 
	}
		
}



int main(int argc, char *argv[]) {  
	char *origin_filename; origin_filename=argv[1];
	
	int numbers[2];
	sprintf(file_path,"./DATA_truth/%s_TT/SPECIES_KEY.csv",origin_filename);
	Read_num_genes(file_path,numbers);
	int num_genes_tot=numbers[1];
	int num_genes_input=numbers[0];

	//if (num_genes_tot<20 && num_genes_tot>0 && numbers[0]==numbers[1]) {
	
	if (numbers[0]<20 && num_genes_tot>0) {
		printf("%d %d\n",numbers[0],numbers[1]);
	
		sprintf(file_path_nw,"./STATE_Networks/%s.nw",origin_filename);
		fpnw=fopen(file_path_nw,"w"); 
		
		unsigned long long num_states = 1ULL << num_genes_input;
		//printf("%llu\n",num_states);

		int input_check[num_genes_tot];
		char key[num_genes_tot][500];
		Read_matrix_header(origin_filename,num_genes_tot,key,input_check);
		//print total key of genes 
		//for(int j=0;j<num_genes_tot;j++) { printf("%s ",key[j]); } printf("\n");
		//for(int i=0;i<num_genes_tot;i++) printf("%d %d\n",i,input_check[i]);

		unsigned long long state,new_state;
		char result;
		char new_binary[num_genes_tot+1]; new_binary[num_genes_tot] = '\0';
		char new_binary_input[num_genes_input+1]; new_binary[num_genes_input] = '\0';

		//time_t now; double seed=(time(&now)); 
		init_genrand64(0); 

		for(state=0;state<num_states;state++) {
			//printf("%llu\t",state);
			
			int binary[num_genes_tot]; 
			decimalToBinary(state,binary,num_genes_tot,input_check); 
			//for(int i=0;i<num_genes_tot;i++) { printf("%d",binary[i]); } printf("\t"); printf("\n");

			for(int i=0;i<num_genes_tot;i++) {
				sprintf(file_path,"./DATA_truth/%s_TT",origin_filename);
				result=Truth_table_ex(file_path,num_genes_tot,state,binary,key,key[i]);
				if(result!='\0') new_binary[i]=result;
				else new_binary[i] = binary[i] + '0'; // 이진 값을 문자로 변환하여 저장
			}
			new_binary[num_genes_tot] = '\0'; // null 문자 추가
			
			int ii=0;
			for(int i=0;i<num_genes_tot;i++) {
				if(input_check[i]==1) {
					new_binary_input[ii]=new_binary[i];
					ii++;
				}
			}
			//printf("%s\n",new_binary_input);

			new_state=(int)strtol(new_binary_input, NULL, 2);
			//printf("--> %llu\t",new_state);
			//printf("%s\n",new_binary);

			fprintf(fpnw,"%llu %llu\n",state,new_state);
		}


		fclose(fpnw);
	}


	return 0;
}


