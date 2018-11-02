void print_mat(char *init, char matrix[256][256]){
  int i,j;
  printf("%s {\n",init);

  printf("/* ");
  for(j=0;j<256;j++)
    printf("%c  ",(j>='A'&&j<='Z')?j:' ');
  printf(" */\n");
    
  for(i=0;i<256;i++){
    printf(" { ");
    for(j=0;j<256;j++)
      printf("%d%s ",matrix[i][j], (j==255)?"},":",");
    if(i>='A'&&i<='Z') printf("/* %c */",i);
    printf(" \n");
  }
  printf("};\n\n");
}

int main(){
  char matrix[256][256];
  char row[256];
  char tmp, *chars;
  int i,j,k;


  /*DNA*/
  for(i=0;i<256;i++)
    for(j=0;j<256;j++)
      matrix[i][j]=1;

  j='A';
  chars="AMRWVHDXN\0";
  for(k=0;k<strlen(chars);k++)
    matrix[chars[k]][j]=0;
  j='C';
  chars="CMSYVHBXN";
  for(k=0;k<strlen(chars);k++)
    matrix[chars[k]][j]=0;
  j='G';
  chars="GRSKVDBXN";
  for(k=0;k<strlen(chars);k++)
    matrix[chars[k]][j]=0;
  j='T';
  chars="TWYKHDBXN";
  for(k=0;k<strlen(chars);k++)
    matrix[chars[k]][j]=0;

  print_mat("static char IUPAC_dna_decode[256][256] =",matrix);
  
  for(i=0;i<256;i++){
    tmp=matrix[i]['A'];
    matrix[i]['A']=matrix[i]['T'];
    matrix[i]['T']=tmp;
    tmp=matrix[i]['C'];
    matrix[i]['C']=matrix[i]['G'];
    matrix[i]['G']=tmp;
  }

  print_mat("static char IUPAC_dna_compl_decode[256][256] =", matrix);


  /*Protein*/
  for(i=0;i<256;i++)
    for(j=0;j<256;j++)
      matrix[i][j]=1;

  chars="ACDEFGHIKLMNPQRSTVWY";
  for(k=0;k<strlen(chars);k++)
    matrix[chars[k]][chars[k]]=0;

  matrix['B']['D']=0;
  matrix['B']['N']=0;
  matrix['Z']['E']=0;
  matrix['Z']['Q']=0;

  print_mat("static char IUPAC_pro_decode[256][256] =", matrix);

  return 0;

}


