#include<stdio.h>
#include<stdlib.h>
main(){
  char c;
  printf("do you really want to make your laptop hot, type y or n\n");
  scanf("%c", &c);
  if(c == 'y'){
    while(1)
      ;
  }
  else{
    exit(0);
  }
}
